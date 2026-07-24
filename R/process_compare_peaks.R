#' Process compare peaks
#'
#' @export
#'
#' @include cascade_params.R
#' @include cascade_run.R
#' @include check_export_dir.R
#' @include compare_peaks.R
#' @include extract_ms_progress.R
#' @include extract_ms_peak.R
#' @include load_ms_data.R
#' @include transform_ms.R
#'
#' @description
#' Scores every feature that was matched to a detector peak by comparing the
#' shape of the detector peak with the shape of that feature's MS1 extracted ion
#' chromatogram, then writes two tables: features that a peak informs, and
#' features that no peak informs.
#'
#' @details
#' This function used to expose only 15 of the pipeline's parameters and ran the
#' rest at helper defaults, which meant a peak set tuned and validated in
#' [check_peaks_integration()] could not be reproduced on export. It now shares
#' the full parameter surface, and both functions run the same
#' [cascade_run()] front half.
#'
#' For a stronger guarantee still, build the run once and hand it to both — then
#' the exported peaks are the very objects you inspected, not merely peaks
#' computed from the same settings:
#'
#' ```r
#' run <- cascade_run(file, features, detector = "cad", load_ms = TRUE)
#' check_peaks_integration(run)
#' process_compare_peaks(run)
#' ```
#'
#' A `<name>_params_<detector>.tsv` sidecar is written next to the tables,
#' recording every resolved parameter for reproducibility.
#'
#' @param file File path, or a `cascade_run` object as returned by
#'   [cascade_run()].
#' @param features Features path
#' @param run An existing `cascade_run`. When supplied, every other tuning
#'   argument is ignored. It must have been built with `shapes = TRUE` and
#'   `load_ms = TRUE`.
#' @param chromatogram Which processed trace to use. One of "original",
#'   "improved", "baselined".
#' @param type Deprecated. Use `chromatogram`.
#' @param detector Detector type. One of "cad", "bpi", "pda".
#' @param headers Named vector mapping detector types to header names.
#' @param export_dir Export directory
#' @param show_example Show example? Default to FALSE
#' @param fourier_components Fraction of Fourier components to keep.
#' @param frequency Acquisition frequency in Hz.
#' @param min_area Relative area cutoff.
#' @param min_area_absolute Absolute area cutoff.
#' @param min_intensity Minimum intensity for feature filtering.
#' @param resample Resampling factor.
#' @param shift Detector-vs-MS time offset, minutes.
#' @param time_min Time min in minutes.
#' @param time_max Time max in minutes.
#' @param intensity_floor Small positive value for intensity floor.
#' @param k2 K2 parameter for signal sharpening.
#' @param k4 K4 parameter for signal sharpening.
#' @param sigma Sigma parameter for signal sharpening.
#' @param smoothing_width Smoothing width, in grid points.
#' @param smoothing_width_minutes Smoothing width in minutes; overrides
#'   `smoothing_width` when supplied.
#' @param baseline_method Method for baseline correction.
#' @param sd_max Maximum fitted peak width, in grid points.
#' @param sd_max_minutes Same, in minutes; overrides `sd_max` when supplied.
#' @param max_iter Maximum iterations for peak fitting.
#' @param min_peak_height Detection sensitivity, as a fraction of the trace
#'   maximum.
#' @param noise_threshold Deprecated and inert.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw".
#' @param intensity_threshold Minimum normalized intensity kept before shape
#'   comparison.
#' @param improve_signal Logical. Whether to apply signal improvement.
#'
#' @return Invisibly, a list with `informed`, `not_informed`, `files` and
#'   `params`. Previously this function ended on a write call and so returned
#'   `NULL`.
#'
#' @examples
#' \dontrun{
#' process_compare_peaks(show_example = TRUE)
#' }
process_compare_peaks <- function(
  file = NULL,
  features = NULL,
  run = NULL,
  chromatogram = cascade_defaults$chromatogram,
  type = NULL,
  detector = cascade_defaults$detector,
  headers = cascade_defaults$headers,
  export_dir = cascade_defaults$export_dir,
  show_example = FALSE,
  fourier_components = cascade_defaults$fourier_components,
  frequency = cascade_defaults$frequency,
  min_area = cascade_defaults$min_area,
  min_area_absolute = cascade_defaults$min_area_absolute,
  min_intensity = cascade_defaults$min_intensity,
  resample = cascade_defaults$resample,
  shift = cascade_defaults$shift,
  time_min = cascade_defaults$time_min,
  time_max = cascade_defaults$time_max,
  intensity_floor = cascade_defaults$intensity_floor,
  k2 = cascade_defaults$k2,
  k4 = cascade_defaults$k4,
  sigma = cascade_defaults$sigma,
  smoothing_width = cascade_defaults$smoothing_width,
  smoothing_width_minutes = cascade_defaults$smoothing_width_minutes,
  baseline_method = cascade_defaults$baseline_method,
  sd_max = cascade_defaults$sd_max,
  sd_max_minutes = cascade_defaults$sd_max_minutes,
  max_iter = cascade_defaults$max_iter,
  min_peak_height = cascade_defaults$min_peak_height,
  noise_threshold = NULL,
  fit = cascade_defaults$fit,
  intensity_threshold = cascade_defaults$intensity_threshold,
  improve_signal = cascade_defaults$improve_signal
) {
  ## `type` was this function's name for what check_peaks_integration() called
  ## `chromatogram`. One name now; the old one still works for a release.
  if (!is.null(type)) {
    deprecate_arg("type", "chromatogram", "process_compare_peaks")
    chromatogram <- type
  }

  ## ---- 1. Resolve the run -------------------------------------------------
  if (inherits(file, "cascade_run")) {
    run <- file
    file <- NULL
  }

  if (is.null(run)) {
    params <- cascade_params(
      detector = detector,
      chromatogram = chromatogram,
      headers = headers,
      time_min = time_min,
      time_max = time_max,
      shift = shift,
      frequency = frequency,
      resample = resample,
      improve_signal = improve_signal,
      fourier_components = fourier_components,
      intensity_floor = intensity_floor,
      smoothing_width = smoothing_width,
      smoothing_width_minutes = smoothing_width_minutes,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      baseline_method = baseline_method,
      fit = fit,
      sd_max = sd_max,
      sd_max_minutes = sd_max_minutes,
      max_iter = max_iter,
      min_peak_height = min_peak_height,
      min_area = min_area,
      min_area_absolute = min_area_absolute,
      min_intensity = min_intensity,
      intensity_threshold = intensity_threshold,
      noise_threshold = noise_threshold
    )

    run <- cascade_run(
      file = file,
      features = features,
      params = params,
      show_example = show_example,
      ## Keep the example fast: the EIC extraction below dominates runtime.
      subsample = if (show_example) 10L else NULL,
      shapes = TRUE,
      load_ms = TRUE
    )
  }

  ## A run built for QC (shapes = FALSE) lacks the objects this function needs.
  ## Say so plainly instead of failing on a NULL subscript later.
  if (!isTRUE(run$shapes)) {
    stop(
      "This `run` was built with `shapes = FALSE` and cannot be scored. ",
      "Rebuild it with `cascade_run(..., shapes = TRUE, load_ms = TRUE)`.",
      call. = FALSE
    )
  }
  if (is.null(run$ms_data)) {
    stop(
      "This `run` was built without MS1 data. ",
      "Rebuild it with `cascade_run(..., load_ms = TRUE)`.",
      call. = FALSE
    )
  }

  params <- run$params
  peaks_prelist <- run$peaks
  name <- run$name

  ## ---- 2. Pull one MS1 EIC per matched peak (the bottleneck) --------------
  message("processing ", params$detector, " peaks")
  message("extracting ms chromatograms (longest step)")
  message(
    "count approx 1 minute per worker per 1000 features (increasing with features number)"
  )
  message("varies a lot depending on features distribution")
  list_ms_chromatograms <- seq_along(
    peaks_prelist$list_df_features_with_peaks_long
  ) |>
    extract_ms_progress(
      ms_data = run$ms_data,
      rts = peaks_prelist$list_rtr,
      mzs = peaks_prelist$list_mzr
    )

  ## Wrap each normalised EIC back into an MSnbase::Chromatogram S4 object,
  ## the input type compareChromatograms() requires.
  message("extracting ms peaks")
  list_ms_peaks <- list_ms_chromatograms |>
    purrr::map(
      .f = extract_ms_peak
    )

  ## ---- 3. Score detector peak vs MS1 EIC ---------------------------------
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
  ## that preserves the column schema, so the select() below still works.
  if (length(list_df_features_with_scores) == 0) {
    df_features_with_peaks <- peaks_prelist$df_features_without_peaks[0, ]
    df_features_with_peaks$comparison_score <- numeric(0)
  } else {
    df_features_with_peaks <- list_df_features_with_scores |>
      tidytable::bind_rows()
  }

  ## ---- 4. Two tidy output tables ------------------------------------------
  ## Identical 11-column schema, so they can be row-bound downstream; the only
  ## difference is that unmatched features carry comparison_score = NA.
  message("final aesthetics")
  df_features_with_peaks_scored <- df_features_with_peaks |>
    select_output_columns()

  df_features_without_peaks_scored <-
    peaks_prelist$df_features_without_peaks |>
    tidytable::mutate(comparison_score = NA) |>
    select_output_columns()

  ## ---- 5. Export ----------------------------------------------------------
  message("checking export directory")
  check_export_dir(export_dir)

  message("exporting")
  stem <- name |>
    gsub(pattern = "\\.[^.]+$", replacement = "")

  path_informed <- file.path(
    export_dir,
    paste0(paste(stem, "featuresInformed", params$detector, sep = "_"), ".tsv")
  )
  path_not_informed <- file.path(
    export_dir,
    paste0(
      paste(stem, "featuresNotInformed", params$detector, sep = "_"),
      ".tsv"
    )
  )
  path_params <- file.path(
    export_dir,
    paste0(paste(stem, "params", params$detector, sep = "_"), ".tsv")
  )

  df_features_with_peaks_scored |>
    tidytable::fwrite(file = path_informed, sep = "\t")

  df_features_without_peaks_scored |>
    tidytable::fwrite(file = path_not_informed, sep = "\t")

  ## Provenance sidecar: the exact settings that produced these two tables.
  params |>
    params_as_table() |>
    tidytable::fwrite(file = path_params, sep = "\t")

  invisible(list(
    informed = df_features_with_peaks_scored,
    not_informed = df_features_without_peaks_scored,
    files = c(
      informed = path_informed,
      not_informed = path_not_informed,
      params = path_params
    ),
    params = params
  ))
}

#' Select and rename the exported columns
#'
#' @description
#' Both output tables must share one schema so they can be row-bound downstream.
#' Defining it once removes the risk of the two `select()` blocks drifting apart.
#'
#' @param df A scored feature table.
#'
#' @return A data frame with the 11 export columns.
#'
#' @keywords internal
#'
#' @examples NULL
select_output_columns <- function(df) {
  df |>
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
}
