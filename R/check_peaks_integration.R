#' Check chromatograms alignment
#'
#' @export
#'
#' @include load_chromatograms.R
#' @include load_features.R
#' @include load_name.R
#' @include plot_peak_detection.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#'
#' @param file File path
#' @param features Features path
#' @param detector Detector type (e.g., "cad", "bpi", "pda")
#' @param chromatogram Chromatogram type. One of "original", "improved", or
#'   "baselined". Default is "baselined".
#' @param headers Named vector mapping detector types to header names.
#' @param min_area Minimum area fraction for peak filtering. Default is 0.005.
#' @param min_intensity Minimum intensity for feature filtering. Default is
#'   1E4.
#' @param shift Time shift in minutes. Default is 0.05.
#' @param show_example Show example data? Default is FALSE.
#' @param fourier_components Fraction of Fourier components to keep. Default is
#'   0.01.
#' @param time_min Time min in minutes. Default is 0.5.
#' @param time_max Time max in minutes. Default is 32.5.
#' @param frequency Acquisition frequency in Hz. Default is 1.
#' @param resample Resampling factor. Default is 1.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection".
#' @param sd_max Maximum standard deviation for peak filtering. Default is 50.
#' @param max_iter Maximum iterations for peak fitting. Default is 1000.
#' @param noise_threshold Noise threshold for peak detection. Default is 0.001.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw". Default
#'   is "egh".
#' @param intensity_threshold Minimum normalized intensity threshold for
#'   filtering. Default is 0.1.
#' @param improve_signal Logical. Whether to apply signal improvement. Default
#'   is TRUE.
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples
#' \dontrun{
#' check_peaks_integration(show_example = TRUE)
#' }
check_peaks_integration <- function(
  file = NULL,
  features = NULL,
  detector = "cad",
  chromatogram = "baselined",
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  min_area = 0.005,
  min_intensity = 1E4,
  shift = 0.05,
  show_example = FALSE,
  fourier_components = 0.01,
  time_min = 0.5,
  time_max = 32.5,
  frequency = 1,
  resample = 1,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh",
  intensity_threshold = 0.1,
  improve_signal = TRUE
) {
  ## ---- 1. Load raw inputs ------------------------------------------------
  ## Named list (bpi/pda/cad) of detector traces read from the mzML via mzR.
  ## Only the chromatogram channels listed in `headers` are pulled.
  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example, headers = headers)

  ## Sample name = basename(file). Needed because the MZmine export is
  ## multi-sample: columns are keyed "datafile:<name>:rt_range:min" etc.
  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  ## Raw MZmine "comprehensive export" table (one row per feature).
  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  ## Reshape to cascade's schema, keep status == "DETECTED" and
  ## intensity_max >= min_intensity, then setkey(rt_1, rt_2) so the
  ## interval join in preprocess_peaks() can use it directly.
  message("preparing features")
  df_features <- feature_table |>
    prepare_features(min_intensity = min_intensity, name = name)

  ## ---- 2. Signal processing on the chosen detector ------------------------
  ## NB: `switch` and `list` shadow the base-R functions of the same name.
  ## Legacy style, harmless in this scope. `switch` holds a *named* vector
  ## (name = "cad", value = "UV#1_CAD_1_0"); names() recovers the short key
  ## used to index into chromatograms_all.
  message("Preprocessing chromatograms")
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
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width,
    baseline_method = baseline_method,
    improve_signal = improve_signal
  )
  ## Returns 6 objects: {original, improved, baselined} x {wide, _long}.
  ## "wide"  = list of (rtime, intensity) data.frames -> fed to the NLS fitter
  ## "_long" = row-bound with an `id` column + degenerate rt_1/rt_2 interval
  ##           columns -> fed to the data.table interval joins.

  ## ---- 3. Detect peaks and match them to features -------------------------
  ## `chromatogram` picks which of the three trace versions to work on.
  ## df_long is height-normalised to [0, 1]; df_xy is the same trace in wide
  ## form ([[1]] unwraps the single-sample list).
  peaks <-
    preprocess_peaks(
      df_features = df_features,
      df_long = switch(
        chromatogram,
        "original" = chromatograms_list$chromatograms_original_long,
        "improved" = chromatograms_list$chromatograms_improved_long,
        "baselined" = chromatograms_list$chromatograms_baselined_long
      ) |>
        tidytable::mutate(intensity = intensity / max(intensity)),
      df_xy = switch(
        chromatogram,
        "original" = chromatograms_list$chromatograms_original[[1]],
        "improved" = chromatograms_list$chromatograms_improved[[1]],
        "baselined" = chromatograms_list$chromatograms_baselined[[1]]
      ),
      min_area = min_area,
      shift = shift,
      name = name,
      sd_max = sd_max,
      max_iter = max_iter,
      noise_threshold = noise_threshold,
      fit = fit,
      intensity_threshold = intensity_threshold
    )

  ## ---- 4. Build the plotting tables ---------------------------------------
  ## One long, height-normalised table for the signal line.
  ## Quirk: the "original" branch also keeps every 10th point only, to thin
  ## the un-resampled raw trace. improved/baselined are already resampled.
  chromatogram_normalized <- switch(
    chromatogram,
    "original" = chromatograms_list$chromatograms_original_long |>
      tidytable::bind_rows() |>
      tidytable::filter(row_number() %% 10 == 1) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "improved" = chromatograms_list$chromatograms_improved_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "baselined" = chromatograms_list$chromatograms_baselined_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity))
  )

  ## Detected-peak markers, put on the same [0, 1] axis as the signal.
  peaks_normalized <- peaks$list_df_features_with_peaks_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(
      intensity = intensity_max / max(intensity_max),
      peak_max = peak_max / max(peak_max)
    )

  ## Linear interpolator over the trace (~ scipy.interpolate.interp1d).
  ## plot_peak_detection() calls it to read the trace height at a peak's
  ## start/end retention time, so the boundary markers sit on the curve.
  approx_f <- stats::approxfun(
    x = chromatogram_normalized |>
      tidytable::pull(rtime),
    y = chromatogram_normalized |>
      tidytable::pull(intensity)
  )

  ## ---- 5. Draw ------------------------------------------------------------
  ## Last expression = return value: an interactive plotly widget.
  chromatogram_normalized |>
    plot_peak_detection(df2 = peaks_normalized, fun = approx_f)
}
