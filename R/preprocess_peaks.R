#' Preprocess peaks
#'
#' @include cascade_defaults.R
#' @include join_peaks.R
#' @include normalize_chromato.R
#' @include peaks_progress.R
#' @include prepare_mz.R
#' @include prepare_peaks.R
#' @include prepare_rt.R
#'
#' @param detector Detector
#' @param df_features DF features
#' @param df_long DF long
#' @param df_xy DF X Y
#' @param name Name
#' @param shift Detector-vs-MS time offset, minutes.
#' @param min_area Relative area cutoff: peak integral divided by the sum of all
#'   peak integrals in the sample.
#' @param min_area_absolute Absolute area cutoff, applied in addition to
#'   `min_area`.
#' @param sd_max Maximum fitted peak width, in grid points.
#' @param max_iter Maximum iterations for peak fitting.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw".
#' @param min_peak_height Detection sensitivity, as a fraction of the trace
#'   maximum.
#' @param intensity_threshold Minimum normalized intensity threshold for
#'   filtering in normalize_chromato.
#' @param shapes Compute the per-peak normalised shapes used for MS comparison?
#'   `FALSE` skips `normalize_chromato()`, `prepare_peaks()`, `prepare_rt()` and
#'   `prepare_mz()`, which a purely visual QC pass never reads.
#'
#' @return A list of lists and dataframe with preprocessed peaks
#'
#' @examples NULL
preprocess_peaks <- function(
  detector = cascade_defaults$detector,
  df_features,
  df_long,
  df_xy,
  name,
  shift = cascade_defaults$shift,
  min_area = cascade_defaults$min_area,
  min_area_absolute = cascade_defaults$min_area_absolute,
  sd_max = cascade_defaults$sd_max,
  max_iter = cascade_defaults$max_iter,
  fit = cascade_defaults$fit,
  min_peak_height = cascade_defaults$min_peak_height,
  intensity_threshold = cascade_defaults$intensity_threshold,
  shapes = TRUE
) {
  message("preprocessing ", detector, " peaks")
  ## data.table call outside of future because buggy else
  peaks <- peaks_progress(
    df_xy = df_xy,
    sd_max = sd_max,
    max_iter = max_iter,
    fit = fit,
    min_peak_height = min_peak_height
  )

  ## data.table call outside of future because buggy else
  peaks_long <- tidytable::bind_rows(peaks, .id = "id") |>
    tidytable::data.table()

  message("joining peaks")
  df_peaks <-
    join_peaks(
      chromatograms = df_long,
      peaks = peaks_long,
      min_area = min_area,
      min_area_absolute = min_area_absolute
    )

  data.table::setkey(df_peaks, rt_min, rt_max)

  message("joining within given rt tolerance")
  df_features_peaks <-
    data.table::foverlaps(df_features, df_peaks)

  df_features_with_peaks <- df_features_peaks |>
    tidytable::select(-rt_1, -rt_2) |>
    tidytable::filter(!is.na(peak_id)) |>
    tidytable::distinct()

  message("selecting features outside peaks")
  df_features_without_peaks <- df_features_peaks |>
    tidytable::filter(is.na(peak_id)) |>
    tidytable::distinct()

  message("splitting by file")
  list_df_features_with_peaks <- df_features_with_peaks |>
    tidytable::group_split(id)

  names(list_df_features_with_peaks) <-
    unique(df_features_with_peaks$id)

  message("splitting by peak")
  list_df_features_with_peaks_per_peak <- list_df_features_with_peaks |>
    purrr::map(
      .f = function(x) {
        x <- x |>
          tidytable::group_split(peak_id)
        return(x)
      }
    )

  list_df_features_with_peaks_long <-
    list_df_features_with_peaks_per_peak |>
    purrr::flatten()

  ## The four artifacts below are only consumed by the MS comparison path.
  ## check_peaks_integration() plots none of them, so it asks for shapes = FALSE
  ## rather than paying for work it discards.
  if (shapes) {
    message("normalizing chromato")
    list_chromato_with_peak <- list_df_features_with_peaks_long |>
      purrr::map(
        .f = normalize_chromato,
        df_xy = df_xy,
        intensity_threshold = intensity_threshold
      )

    message("preparing peaks chromato")
    list_chromato_peaks <- list_chromato_with_peak |>
      purrr::map(
        .f = prepare_peaks
      )

    message("preparing rt")
    list_rtr <- list_df_features_with_peaks_long |>
      purrr::map(
        .f = prepare_rt,
        shift = shift
      )

    message("preparing mz")
    list_mzr <- list_df_features_with_peaks_long |>
      purrr::map(
        .f = prepare_mz
      )
  } else {
    list_chromato_peaks <- NULL
    list_rtr <- NULL
    list_mzr <- NULL
  }

  returned_list <- list(
    list_df_features_with_peaks_long,
    list_chromato_peaks,
    list_rtr,
    list_mzr,
    df_features_without_peaks
  )
  names(returned_list) <- c(
    "list_df_features_with_peaks_long",
    "list_chromato_peaks",
    "list_rtr",
    "list_mzr",
    "df_features_without_peaks"
  )
  return(returned_list)
}
