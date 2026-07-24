#' Signal sharpening
#'
#' @include cascade_defaults.R
#' @include second_der.R
#'
#' @param time time
#' @param intensity intensity
#' @param k2 Divisor of the second-derivative term. Lower values increase the
#'   sharpening effect from the second derivative.
#' @param k4 Divisor of the fourth-derivative term. Lower values increase the
#'   sharpening effect from the fourth derivative.
#' @param sigma Overall sharpening gain. Higher values increase the effect.
#' @param smoothing_width Smoothing width for the running mean filter, in grid
#'   points. Higher values provide more smoothing but reduce resolution.
#'   Formerly spelled `Smoothing_width`.
#' @param baseline_adjust Baseline adjustment value. Formerly spelled
#'   `Baseline_adjust`.
#'
#' @return A sharpened signal
#'
#' @examples NULL
signal_sharpening <- function(
  time,
  intensity,
  k2 = cascade_defaults$k2,
  k4 = cascade_defaults$k4,
  sigma = cascade_defaults$sigma,
  smoothing_width = cascade_defaults$smoothing_width,
  baseline_adjust = 0
) {
  smooth_1 <- caTools::runmean(
    x = intensity,
    k = smoothing_width,
    align = "center"
  ) +
    baseline_adjust

  smooth_2 <- caTools::runmean(
    x = smooth_1,
    k = smoothing_width,
    align = "center"
  )

  deriv_2 <- second_der(
    x = time,
    y = smooth_2
  )

  smooth_3 <- caTools::runmean(
    x = deriv_2,
    k = smoothing_width,
    align = "center"
  )

  deriv_4 <- second_der(
    x = time[3:length(time)],
    y = smooth_3
  )

  smooth_4 <- caTools::runmean(
    x = deriv_4,
    k = smoothing_width,
    align = "center"
  )

  sharpened <- smooth_1[5:length(smooth_1)] -
    (sigma / k2 * smooth_3[3:length(smooth_3)]) +
    (sigma / k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  # sharpened <- sharpened / max(sharpened)

  return(sharpened)
}
