#' Improve signals progress
#'
#' @include cascade_defaults.R
#' @include improve_signal.R
#'
#' @param xs List of dataframes with 'rtime' and 'intensity' columns
#' @param fourier_components Fraction of Fourier components to keep. Default is
#'   0.01.
#' @param frequency Acquisition frequency in Hz. Default is 2.
#' @param resample Resampling factor. Default is 1.
#' @param time_min Time min in minutes. Default is 0.
#' @param time_max Time max in minutes. Default is Inf.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#'
#' @return A list of data frames with improved signals
#'
#' @examples NULL
improve_signals_progress <- function(
  xs,
  fourier_components = cascade_defaults$fourier_components,
  frequency = cascade_defaults$frequency,
  resample = cascade_defaults$resample,
  time_min = cascade_defaults$time_min,
  time_max = cascade_defaults$time_max,
  intensity_floor = cascade_defaults$intensity_floor,
  k2 = cascade_defaults$k2,
  k4 = cascade_defaults$k4,
  sigma = cascade_defaults$sigma,
  smoothing_width = cascade_defaults$smoothing_width
) {
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(
        x,
        fourier_components,
        frequency,
        resample,
        time_min,
        time_max,
        intensity_floor,
        k2,
        k4,
        sigma,
        smoothing_width
      ) {
        improve_signal(
          df = x |>
            tidytable::select(rtime, intensity),
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample,
          time_min = time_min,
          time_max = time_max,
          intensity_floor = intensity_floor,
          k2 = k2,
          k4 = k4,
          sigma = sigma,
          smoothing_width = smoothing_width
        )
      },
      fourier_components = fourier_components,
      frequency = frequency,
      resample = resample,
      time_min = time_min,
      time_max = time_max,
      intensity_floor = intensity_floor,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      smoothing_width = smoothing_width
    )
}
