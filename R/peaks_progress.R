#' Peaks progress
#'
#' @include cascade_defaults.R
#' @include get_peaks.R
#'
#' @param df_xy Df X Y
#' @param sd_max Maximum fitted peak width, in grid points.
#' @param max_iter Maximum iterations for peak fitting.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw".
#' @param min_peak_height Detection sensitivity, as a fraction of the trace
#'   maximum. `0` keeps every candidate the derivative test finds.
#'
#' @return A list of peaks
#'
#' @examples NULL
peaks_progress <- function(
  df_xy,
  sd_max = cascade_defaults$sd_max,
  max_iter = cascade_defaults$max_iter,
  fit = cascade_defaults$fit,
  min_peak_height = cascade_defaults$min_peak_height
) {
  ## get_peaks() expects a matrix with retention times as row names and one
  ## column per wavelength; "666" is a placeholder for the single channel.
  matrix_666 <- df_xy |>
    tidytable::filter(rtime >= 0) |>
    tidytable::select(rtime, intensity) |>
    tidytable::rename(`666` = intensity) |>
    tibble::column_to_rownames("rtime") |>
    as.matrix()

  ## `min_peak_height` is a fraction of the trace maximum; find_peaks() wants an
  ## absolute amplitude, so convert here. This is the detection-sensitivity knob
  ## that was previously unreachable from any exported function.
  amp_thresh <- min_peak_height * max(matrix_666, na.rm = TRUE)

  list("666" = matrix_666) |>
    get_peaks(
      lambdas = c("666"),
      fit = fit,
      sd_max = sd_max,
      max_iter = max_iter,
      amp_thresh = amp_thresh
    ) |>
    purrr::pluck("666") |>
    purrr::pluck("666") |>
    tidytable::mutate(
      peak_id = tidytable::row_number(),
      peak_max = height,
      rt_apex = rt,
      rt_min = start,
      rt_max = end
    )
}
