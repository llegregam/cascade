#' Baseline chromatogram
#'
#' @include cascade_defaults.R
#'
#' @param df Dataframe
#' @param method Baseline correction method. See
#'   \code{\link[baseline]{baseline}} for available methods including: "als",
#'   "fillPeaks", "irls", "lowpass", "medianWindow", "modpolyfit", "peakDetection",
#'   "rfbaseline", "rollingBall", "shirley", "TAP".
#' @param ... Additional arguments passed to \code{\link[baseline]{baseline}}.
#'
#' @return A dataframe with baselined chromatogram
#'
#' @examples NULL
baseline_chromatogram <- function(
  df,
  method = cascade_defaults$baseline_method,
  ...
) {
  df2 <- df

  intensity <- df$intensity

  intensity[is.na(intensity)] <- 0

  intensity_baseline <- baseline::baseline(
    spectra = t(intensity),
    method = method,
    ...
  )

  intensity_new <- t(intensity_baseline@corrected) |>
    tidytable::data.table()

  df2$intensity <- intensity_new$V1

  return(df2)
}
