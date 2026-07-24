#' Download the latest LOTUS release
#'
#' @description
#' Thin guarded wrapper around `tima::get_last_version_from_zenodo()`.
#'
#' `tima` sits in `Suggests` rather than `Imports` because it is needed only by
#' the annotation-side functions (`generate_ids()`, `generate_tables()`,
#' `generate_pseudochromatograms()`). The peak-integration and alignment path —
#' `cascade_run()`, `check_peaks_integration()`, `process_compare_peaks()`,
#' `check_chromatograms_alignment()` — does not touch it. Keeping it optional
#' means the package loads and the integration workflow runs without a working
#' `tima` install.
#'
#' @param destination Where to write the downloaded archive.
#' @param doi Zenodo DOI of the LOTUS deposit.
#' @param pattern File name pattern to pull from the deposit.
#'
#' @return The value of `tima::get_last_version_from_zenodo()`, invisibly.
#'
#' @keywords internal
#'
#' @examples NULL
get_lotus <- function(
  destination = "data/source/libraries/lotus.csv.gz",
  doi = "10.5281/zenodo.5794106",
  pattern = "frozen_metadata.csv.gz"
) {
  if (!requireNamespace("tima", quietly = TRUE)) {
    stop(
      "This function needs the `tima` package, which is optional.\n",
      "Install it with:\n",
      "  remotes::install_github(\"taxonomicallyinformedannotation/tima\")\n",
      "The peak-integration and alignment workflow does not require it.",
      call. = FALSE
    )
  }
  message("Getting last LOTUS version")
  invisible(
    tima::get_last_version_from_zenodo(
      doi = doi,
      pattern = pattern,
      destination
    )
  )
}
