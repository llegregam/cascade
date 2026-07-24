#' Process compare peaks
#'
#' @export
#'
#' @include check_export_dir.R
#' @include compare_peaks.R
#' @include extract_ms_progress.R
#' @include extract_ms_peak.R
#' @include load_chromatograms.R
#' @include load_features.R
#' @include load_name.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#' @include transform_ms.R
#'
#' @param file File path
#' @param features Features path
#' @param type Type. "original", "baselined" or "improved"
#' @param detector Detector
#' @param headers Headers
#' @param export_dir Export directory
#' @param show_example Show example? Default to FALSE
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param min_area Min area
#' @param min_intensity Min intensity
#' @param resample Resample
#' @param shift Shift
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples NULL
process_compare_peaks <- function(
  file = NULL,
  features = NULL,
  type = "baselined",
  detector = "cad",
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  export_dir = "data/interim/peaks",
  show_example = FALSE,
  fourier_components = 0.01,
  frequency = 1,
  min_area = 0.005,
  min_intensity = 1E4,
  resample = 1,
  shift = 0.05,
  time_min = 0.5,
  time_max = 32.5
) {
  ## ---- 1. Load raw inputs ------------------------------------------------
  ## Unlike check_peaks_integration(), this function also needs the raw MS1
  ## scans: mode = "onDisk" keeps only an index in RAM and reads spectra on
  ## demand, which is what makes full LC-MS files tractable.
  message("loading MS data")
  ms_data <- file |>
    load_ms_data(show_example = show_example)

  ## Named list (bpi/pda/cad) of detector traces read from the mzML via mzR.
  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example, headers = headers)

  ## Sample name = basename(file); keys this sample's columns in the
  ## multi-sample MZmine export ("datafile:<name>:rt_range:min", ...).
  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  message("preparing features")
  ## Demo mode only: subsample so the (expensive) EIC extraction below stays
  ## fast. Seeded, so the example is reproducible.
  if (show_example) {
    message("selecting 10 random features for the example")
    set.seed(42)
    feature_table <- feature_table |>
      tidytable::slice_sample(n = 10)
  }
  ## Reshape + filter (DETECTED, >= min_intensity) + setkey(rt_1, rt_2).
  df_features <- feature_table |>
    prepare_features(min_intensity = min_intensity, name = name)

  ## ---- 2. Signal processing on the chosen detector ------------------------
  ## NB: `switch`/`list` shadow the base-R functions of the same name.
  ## Note this call omits the DSP knobs check_peaks_integration() exposes
  ## (k2/k4/sigma/baseline_method/...), so preprocess_chromatograms() runs
  ## with its own defaults here.
  message("preprocessing chromatograms")
  switch <- switch(
    detector,
    "bpi" = headers["bpi"],
    "cad" = headers["cad"],
    "pda" = headers["pda"]
  )
  list <- chromatograms_all[switch |> names()]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = list,
    # signal_name = signal_name,
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample
  )

  ## ---- 3. Detect peaks and match them to features -------------------------
  ## Returns 5 artifacts; the ones used below are:
  ##   list_df_features_with_peaks_long : one element per (sample, peak)
  ##   list_chromato_peaks              : detector peak as MSnbase Chromatogram
  ##   list_rtr / list_mzr              : the rt x mz box of each EIC to pull
  ##   df_features_without_peaks        : features that matched no peak
  message("preprocessing peaks")
  peaks_prelist <- preprocess_peaks(
    detector = detector,
    df_features = df_features,
    df_long = switch(
      type,
      "original" = chromatograms_list$chromatograms_original_long,
      "improved" = chromatograms_list$chromatograms_improved_long,
      "baselined" = chromatograms_list$chromatograms_baselined_long
    ) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    df_xy = switch(
      type,
      "original" = chromatograms_list$chromatograms_original[[1]],
      "improved" = chromatograms_list$chromatograms_improved[[1]],
      "baselined" = chromatograms_list$chromatograms_baselined[[1]]
    ),
    min_area = min_area,
    shift = shift,
    name = name
  )

  message("processing ", detector, " peaks")
  message("extracting ms chromatograms (longest step)")
  message(
    "count approx 1 minute per worker per 1000 features (increasing with features number)"
  )
  message("varies a lot depending on features distribution")
  ## ---- 4. Pull one MS1 EIC per matched peak (the bottleneck) --------------
  ## seq_along(...) yields integer indices 1..N; extract_ms_progress() uses
  ## each index to look up that peak's rt window (rts) and m/z windows (mzs).
  ## `nrows` is passed but unused by the inner map (vestigial).
  list_ms_chromatograms <- seq_along(
    peaks_prelist$list_df_features_with_peaks_long
  ) |>
    extract_ms_progress(
      ms_data = ms_data,
      rts = peaks_prelist$list_rtr,
      mzs = peaks_prelist$list_mzr,
      nrows = peaks_prelist$list_df_features_with_peaks_long |>
        purrr::map(.f = nrow)
    )

  ## Wrap each normalised EIC back into an MSnbase::Chromatogram S4 object,
  ## the input type compareChromatograms() requires.
  message("extracting ms peaks")
  list_ms_peaks <- list_ms_chromatograms |>
    purrr::map(
      .f = extract_ms_peak
    )

  ## ---- 5. Score detector peak vs MS1 EIC ---------------------------------
  ## Shape similarity via MSnbase::compareChromatograms(method = "closest").
  ## Both traces were rt-normalised to [0, 1] upstream, so this compares
  ## peak shape rather than absolute retention time.
  message("comparing peaks")
  list_comparison_score <- seq_along(list_ms_peaks) |>
    purrr::map(
      .f = compare_peaks,
      list_ms_peaks = list_ms_peaks,
      peaks_prelist = peaks_prelist
    )

  ## Collapse the nested list to one score per (sample, peak).
  message("summarizing comparison scores")
  n_scores <- sum(lengths(list_comparison_score))

  message("there are ", n_scores, " scores calculated")

  ## Map() is an element-wise zip: staple each score back onto its table.
  message("selecting features with peaks")
  list_df_features_with_scores <- Map(
    function(df, score) {
      df$comparison_score <- unlist(score, use.names = FALSE)
      df
    },
    peaks_prelist$list_df_features_with_peaks_long,
    list_comparison_score
  )

  ## Guard for "no feature matched any peak": `[0, ]` is an empty row slice
  ## that preserves the column schema (like df.iloc[0:0]), so the select()
  ## and export below still work.
  if (length(list_df_features_with_scores) == 0) {
    # Keep expected columns for downstream select/export when no peak match exists.
    df_features_with_peaks <- peaks_prelist$df_features_without_peaks[0, ]
    df_features_with_peaks$comparison_score <- numeric(0)
  } else {
    df_features_with_peaks <- list_df_features_with_scores |>
      tidytable::bind_rows()
  }

  ## ---- 6. Two tidy output tables ------------------------------------------
  ## select(new = old) renames while selecting. Both tables share the exact
  ## same 11 columns so they can be row-bound downstream
  ## (see prepare_comparison()); the only difference is that unmatched
  ## features carry comparison_score = NA.
  message("final aesthetics")
  df_features_with_peaks_scored <- df_features_with_peaks |>
    tidytable::select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = area,
      comparison_score
    ) |>
    tidytable::distinct()

  df_features_without_peaks_scored <-
    peaks_prelist$df_features_without_peaks |>
    tidytable::mutate(comparison_score = NA) |>
    tidytable::select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = area,
      comparison_score
    ) |>
    tidytable::distinct()

  ## ---- 7. Export ----------------------------------------------------------
  ## Hand-rolled os.makedirs(export_dir, exist_ok = TRUE).
  message("checking export directory")
  check_export_dir(export_dir)

  message("exporting")

  ## Filename: gsub("\\.[^.]+$", "") strips the extension off the sample name
  ## (re.sub(r"\.[^.]+$", "", name)), then the tag and detector are pasted on.
  ## No explicit return(): the products of this function are the two files.
  df_features_with_peaks_scored |>
    tidytable::fwrite(
      file = file.path(
        export_dir,
        name |>
          gsub(pattern = "\\.[^.]+$", replacement = "") |>
          paste("featuresInformed", detector, sep = "_") |>
          paste0(".tsv")
      ),
      sep = "\t"
    )

  df_features_without_peaks_scored |>
    tidytable::fwrite(
      file = file.path(
        export_dir,
        name |>
          gsub(pattern = "\\.[^.]+$", replacement = "") |>
          paste("featuresNotInformed", detector, sep = "_") |>
          paste0(".tsv")
      ),
      sep = "\t"
    )
}
