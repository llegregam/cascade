#' Normalize chromato
#'
#' @include cascade_defaults.R
#'
#' @param x X
#' @param df_xy Df X Y
#' @param intensity_threshold Minimum normalized intensity threshold for
#'   filtering. Set to 0 to keep all points.
#'
#' @return A normalized chromato
#'
#' @examples NULL
normalize_chromato <- function(
  x,
  df_xy,
  intensity_threshold = cascade_defaults$intensity_threshold
) {
  df_xy |>
    tidytable::filter(
      rtime >= x$rt_min[1] &
        rtime <= x$rt_max[1]
    ) |>
    tidytable::mutate(
      intensity = (intensity - min(intensity)) /
        (max(intensity) -
          min(intensity))
    ) |>
    tidytable::filter(intensity >= intensity_threshold)
}
