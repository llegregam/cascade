#' Default parameter values for the whole package
#'
#' @export
#'
#' @description
#' Single source of truth for every parameter shared by more than one function.
#'
#' Historically each function carried its own literal defaults, which drifted:
#' `frequency` was `1` in the entry points but `2` in `preprocess_chromatograms()`
#' and `improve_signal()`, `time_max` was `32.5` at the top and `Inf` underneath,
#' and `get_peaks()` used `max.iter = 100` where everything above it used `1000`.
#' The effective default therefore depended on which function you happened to call.
#'
#' Every signature in the package now reads its default from this list, so a value
#' is defined exactly once. Do not add a literal default to a function signature
#' for anything that appears here.
#'
#' @details
#' Units are stated explicitly because they are mixed:
#'
#' \describe{
#'   \item{minutes}{`time_min`, `time_max`, `shift`, `sd_max_minutes`,
#'     `smoothing_width_minutes`}
#'   \item{Hz}{`frequency`}
#'   \item{grid points}{`sd_max`, `smoothing_width` — these rescale silently when
#'     `frequency` or `resample` change, which is why the `_minutes` variants exist}
#'   \item{fraction of a total}{`fourier_components`, `min_area`,
#'     `intensity_threshold`, `min_peak_height`}
#'   \item{absolute detector counts}{`min_intensity`, `min_area_absolute`}
#' }
#'
#' @format A named list.
#'
#' @examples
#' \dontrun{
#' ## every knob and its current value
#' str(cascade_defaults)
#' }
cascade_defaults <- list(
  ## ---- I/O and channel selection -----------------------------------------
  ## Maps the short detector key to the vendor's chromatogramId in the mzML.
  ## Previously repeated as a literal in four different signatures.
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  detector = "cad",
  ## Which of the three processed trace versions to detect peaks on.
  ## Was named `chromatogram` in check_peaks_integration() and `type` in
  ## process_compare_peaks(); `chromatogram` is now canonical.
  chromatogram = "baselined",
  export_dir = "data/interim/peaks",

  ## ---- (A) time axis / acquisition grid ----------------------------------
  time_min = 0.5,
  time_max = 32.5,
  ## Detector-vs-MS dead-volume offset, minutes. Instrument specific.
  shift = 0.05,
  ## Detector acquisition rate, Hz. Grid step = 1 / (frequency * 60 * resample).
  frequency = 1,
  resample = 1,

  ## ---- (B) signal conditioning -------------------------------------------
  improve_signal = TRUE,
  ## Fraction of Fourier coefficients kept by the low-pass. Lower = smoother.
  fourier_components = 0.01,
  intensity_floor = 0.001,
  ## Running-mean window, in GRID POINTS (see smoothing_width_minutes).
  smoothing_width = 8,
  smoothing_width_minutes = NULL,
  ## Sharpening: f - (sigma/k2) f'' + (sigma/k4) f''''
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  baseline_method = "peakDetection",

  ## ---- (C) peak detection and fitting ------------------------------------
  fit = "egh",
  ## Max fitted peak width, in GRID POINTS (see sd_max_minutes).
  sd_max = 50,
  sd_max_minutes = NULL,
  max_iter = 1000,
  ## Detection sensitivity, as a fraction of the trace maximum. 0 = detect
  ## everything the derivative test finds. Previously unreachable: this is
  ## `amp_thresh` inside get_peaks()'s find_peaks(), which no public function
  ## exposed.
  min_peak_height = 0,

  ## ---- (D)/(E) filtering and matching ------------------------------------
  ## RELATIVE area cutoff: peak integral / sum of all peak integrals.
  min_area = 0.005,
  ## ABSOLUTE area cutoff, applied in addition to the relative one.
  min_area_absolute = 0,
  min_intensity = 1E4,
  intensity_threshold = 0.1
)

#' Valid values for the enumerated parameters
#'
#' @description
#' Used by [cascade_params()] for `match.arg()`-style validation. Kept next to
#' [cascade_defaults] so the allowed set and the default cannot drift apart.
#'
#' @format A named list of character vectors.
#'
#' @keywords internal
#'
#' @examples NULL
cascade_choices <- list(
  detector = c("cad", "bpi", "pda"),
  chromatogram = c("baselined", "improved", "original"),
  fit = c("egh", "gaussian", "raw"),
  ## As accepted by baseline::baseline().
  baseline_method = c(
    "als",
    "fillPeaks",
    "irls",
    "lowpass",
    "medianWindow",
    "modpolyfit",
    "peakDetection",
    "rfbaseline",
    "rollingBall",
    "shirley",
    "TAP"
  )
)

#' Warn about a deprecated argument
#'
#' @description
#' Minimal stand-in for `lifecycle::deprecate_warn()`. `lifecycle` sits in
#' Suggests rather than Imports, so relying on it here would mean guarding every
#' call with `requireNamespace()`.
#'
#' @param old Name of the deprecated argument.
#' @param new Name of the replacement argument.
#' @param what Name of the function the argument belongs to.
#'
#' @return `NULL`, invisibly. Called for the warning side effect.
#'
#' @keywords internal
#'
#' @examples NULL
deprecate_arg <- function(old, new, what) {
  warning(
    sprintf(
      "The `%s` argument of `%s()` is deprecated; use `%s` instead. The value you supplied was used.",
      old,
      what,
      new
    ),
    call. = FALSE
  )
  invisible(NULL)
}

#' Grid step in minutes
#'
#' @description
#' The uniform resampling step used by [improve_signal()], and therefore the
#' conversion factor between "grid points" and minutes for every parameter
#' expressed in points (`sd_max`, `smoothing_width`).
#'
#' @param frequency Acquisition frequency in Hz.
#' @param resample Resampling factor.
#'
#' @return A single numeric: minutes per grid point.
#'
#' @keywords internal
#'
#' @examples NULL
grid_step <- function(frequency, resample) {
  1 / (frequency * 60 * resample)
}
