#' Join peaks
#'
#' @include cascade_defaults.R
#'
#' @description
#' Attaches every trace sample inside a peak's RT window to that peak, integrates
#' it, and drops peaks that are too small.
#'
#' @param chromatograms Chromatograms
#' @param peaks Peaks
#' @param min_area **Relative** cutoff: a peak is kept when its integral divided
#'   by the sum of all peak integrals in the sample is at least this. Because it
#'   is relative, a single dominant peak pushes every other peak's fraction down.
#' @param min_area_absolute **Absolute** cutoff on the same integral, applied in
#'   addition to `min_area`. `0` disables it.
#'
#' @return A dataframe with joined peaks
#'
#' @examples NULL
join_peaks <- function(
  chromatograms,
  peaks,
  min_area = cascade_defaults$min_area,
  min_area_absolute = cascade_defaults$min_area_absolute
) {
  data.table::setkey(peaks, rt_min, rt_max)
  data.table::setkey(chromatograms, rt_1, rt_2)

  data.table::foverlaps(peaks, chromatograms) |>
    # tidytable::filter(id == i.id) |>
    tidytable::group_by(peak_id, id) |>
    ## Discrete integration: sum of trace intensities inside the peak bounds.
    tidytable::mutate(integral = sum(intensity)) |>
    tidytable::ungroup() |>
    tidytable::distinct(
      peak_id,
      id,
      peak_max,
      rt_apex,
      rt_min,
      rt_max,
      integral
    ) |>
    tidytable::group_by(id) |>
    tidytable::filter(
      integral / sum(integral) >= min_area &
        integral >= min_area_absolute
    ) |>
    tidytable::data.table()
}
