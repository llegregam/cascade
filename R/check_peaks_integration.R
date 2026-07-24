#' Check peaks integration
#'
#' @export
#'
#' @include cascade_params.R
#' @include cascade_run.R
#' @include plot_peak_detection.R
#'
#' @description
#' Visual QC of peak detection: draws the conditioned detector trace with the
#' detected apex and boundary markers on top. Returns an interactive plotly
#' figure and writes nothing.
#'
#' The detected peak table is attached to the figure as the `peaks` attribute, so
#' you can inspect what was found without re-running anything:
#' `attr(p, "peaks")`.
#'
#' @details
#' Every parameter is validated by [cascade_params()] before any work happens, so
#' a mistyped `chromatogram` or `detector` is reported by name rather than
#' failing several calls deeper.
#'
#' To guarantee that what you look at here is exactly what
#' [process_compare_peaks()] exports, build the run once and pass it to both:
#'
#' ```r
#' run <- cascade_run(file, features, detector = "cad", sigma = 0.08)
#' check_peaks_integration(run)
#' process_compare_peaks(run)
#' ```
#'
#' @param file File path, or a `cascade_run` object as returned by
#'   [cascade_run()].
#' @param features Features path
#' @param run An existing `cascade_run`. When supplied, every other tuning
#'   argument is ignored because the run already carries its own resolved
#'   parameters.
#' @param detector Detector type. One of "cad", "bpi", "pda".
#' @param chromatogram Chromatogram type. One of "original", "improved", or
#'   "baselined".
#' @param headers Named vector mapping detector types to header names.
#' @param min_area Relative area cutoff: peak integral divided by the sum of all
#'   peak integrals in the sample.
#' @param min_area_absolute Absolute area cutoff, applied in addition to
#'   `min_area`.
#' @param min_intensity Minimum intensity for feature filtering.
#' @param shift Detector-vs-MS time offset, minutes.
#' @param show_example Show example data?
#' @param fourier_components Fraction of Fourier components to keep.
#' @param time_min Time min in minutes.
#' @param time_max Time max in minutes.
#' @param frequency Acquisition frequency in Hz.
#' @param resample Resampling factor.
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
#' @param intensity_threshold Minimum normalized intensity threshold used by the
#'   shape comparison. It does not affect this figure.
#' @param improve_signal Logical. Whether to apply signal improvement.
#'
#' @return A plotly figure, with the peak table attached as the `peaks`
#'   attribute.
#'
#' @examples
#' \dontrun{
#' check_peaks_integration(show_example = TRUE)
#' }
check_peaks_integration <- function(
  file = NULL,
  features = NULL,
  run = NULL,
  detector = cascade_defaults$detector,
  chromatogram = cascade_defaults$chromatogram,
  headers = cascade_defaults$headers,
  min_area = cascade_defaults$min_area,
  min_area_absolute = cascade_defaults$min_area_absolute,
  min_intensity = cascade_defaults$min_intensity,
  shift = cascade_defaults$shift,
  show_example = FALSE,
  fourier_components = cascade_defaults$fourier_components,
  time_min = cascade_defaults$time_min,
  time_max = cascade_defaults$time_max,
  frequency = cascade_defaults$frequency,
  resample = cascade_defaults$resample,
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
  ## ---- 1. Resolve the run -------------------------------------------------
  ## A cascade_run may arrive either as `run =` or, for convenience, in `file`
  ## (a run object is never a valid path, so this is unambiguous).
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

    ## shapes = FALSE: the per-peak normalised chromatograms are only consumed
    ## by the MS comparison, and this figure never reads them.
    run <- cascade_run(
      file = file,
      features = features,
      params = params,
      show_example = show_example,
      shapes = FALSE,
      load_ms = FALSE
    )
  }

  ## ---- 2. Build the plotting tables ---------------------------------------
  ## One long, height-normalised table for the signal line. The "original"
  ## branch is additionally thinned to every 10th point, because unlike the
  ## other two it has not been through the resampling step.
  chromatogram_normalized <- pick_chromatogram(
    run$chromatograms,
    run$params$chromatogram,
    long = TRUE
  ) |>
    tidytable::bind_rows()

  if (run$params$chromatogram == "original") {
    chromatogram_normalized <- chromatogram_normalized |>
      tidytable::filter(tidytable::row_number() %% 10 == 1)
  }

  chromatogram_normalized <- chromatogram_normalized |>
    tidytable::mutate(intensity = intensity / max(intensity))

  ## Detected-peak markers, put on the same [0, 1] axis as the signal.
  peaks_normalized <- run$peaks$list_df_features_with_peaks_long |>
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

  ## ---- 3. Draw ------------------------------------------------------------
  figure <- chromatogram_normalized |>
    plot_peak_detection(df2 = peaks_normalized, fun = approx_f)

  ## Attach the peak table so a QC run is inspectable, not just lookable-at.
  attr(figure, "peaks") <- peaks_normalized
  attr(figure, "params") <- run$params
  figure
}
