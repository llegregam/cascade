#' Resolve and validate pipeline parameters
#'
#' @export
#'
#' @include cascade_defaults.R
#'
#' @description
#' Builds the single parameter object that drives the whole pipeline. Every
#' exported entry point calls this as its first statement, so a bad argument is
#' reported at the door with the offending name and value, instead of surfacing
#' several calls deeper as a `switch()` returning `NULL`.
#'
#' It also resolves the two parameters that exist in both grid points and
#' minutes: supply `sd_max_minutes` or `smoothing_width_minutes` and the point
#' value is computed from the resolved grid step, so the filter keeps the same
#' physical meaning when you change `frequency` or `resample`.
#'
#' @param detector Detector key. One of `"cad"`, `"bpi"`, `"pda"`.
#' @param chromatogram Which processed trace to detect peaks on. One of
#'   `"baselined"`, `"improved"`, `"original"`.
#' @param headers Named character vector mapping detector keys to the
#'   `chromatogramId` strings in the mzML.
#' @param time_min,time_max Analysis window, minutes.
#' @param shift Detector-vs-MS time offset, minutes.
#' @param frequency Detector acquisition rate, Hz.
#' @param resample Resampling factor (>= 1).
#' @param improve_signal Apply Fourier filtering and derivative sharpening?
#' @param fourier_components Fraction of Fourier coefficients kept, in (0, 1].
#' @param intensity_floor Small positive offset applied when the trace has
#'   non-positive values.
#' @param smoothing_width Running-mean window, in grid points.
#' @param smoothing_width_minutes Same, in minutes. Overrides `smoothing_width`
#'   when not `NULL`.
#' @param k2,k4 Divisors of the 2nd- and 4th-derivative sharpening terms. Higher
#'   is weaker.
#' @param sigma Sharpening gain. Higher is stronger.
#' @param baseline_method Method passed to \code{\link[baseline]{baseline}}.
#' @param fit Peak model. One of `"egh"`, `"gaussian"`, `"raw"`.
#' @param sd_max Maximum fitted peak width, in grid points.
#' @param sd_max_minutes Same, in minutes. Overrides `sd_max` when not `NULL`.
#' @param max_iter Maximum NLS iterations.
#' @param min_peak_height Detection sensitivity, as a fraction of the trace
#'   maximum, in \[0, 1).
#' @param min_area Relative area cutoff, in \[0, 1].
#' @param min_area_absolute Absolute area cutoff. Applied in addition to
#'   `min_area`.
#' @param min_intensity Absolute cutoff on the MZmine feature intensity.
#' @param intensity_threshold Fraction of peak height below which points are
#'   dropped before shape comparison, in \[0, 1].
#' @param noise_threshold Deprecated and inert. See details.
#'
#' @details
#' `noise_threshold` never had any effect. It was forwarded only to an internal
#' purity estimator that is force-disabled for single-channel input, which is
#' always the case in this pipeline. Supplying it now warns.
#'
#' @return An object of class `cascade_params`: a named list with every
#'   parameter resolved to a final value.
#'
#' @examples
#' \dontrun{
#' cascade_params(detector = "cad", sigma = 0.08, sd_max_minutes = 0.5)
#' }
cascade_params <- function(
  detector = cascade_defaults$detector,
  chromatogram = cascade_defaults$chromatogram,
  headers = cascade_defaults$headers,
  time_min = cascade_defaults$time_min,
  time_max = cascade_defaults$time_max,
  shift = cascade_defaults$shift,
  frequency = cascade_defaults$frequency,
  resample = cascade_defaults$resample,
  improve_signal = cascade_defaults$improve_signal,
  fourier_components = cascade_defaults$fourier_components,
  intensity_floor = cascade_defaults$intensity_floor,
  smoothing_width = cascade_defaults$smoothing_width,
  smoothing_width_minutes = cascade_defaults$smoothing_width_minutes,
  k2 = cascade_defaults$k2,
  k4 = cascade_defaults$k4,
  sigma = cascade_defaults$sigma,
  baseline_method = cascade_defaults$baseline_method,
  fit = cascade_defaults$fit,
  sd_max = cascade_defaults$sd_max,
  sd_max_minutes = cascade_defaults$sd_max_minutes,
  max_iter = cascade_defaults$max_iter,
  min_peak_height = cascade_defaults$min_peak_height,
  min_area = cascade_defaults$min_area,
  min_area_absolute = cascade_defaults$min_area_absolute,
  min_intensity = cascade_defaults$min_intensity,
  intensity_threshold = cascade_defaults$intensity_threshold,
  noise_threshold = NULL
) {
  ## ---- enumerated parameters --------------------------------------------
  ## match.arg() with an explicit choices vector, so a typo names itself
  ## instead of silently producing NULL from a downstream switch().
  detector <- match_choice(detector, "detector")
  chromatogram <- match_choice(chromatogram, "chromatogram")
  fit <- match_choice(fit, "fit")
  baseline_method <- match_choice(baseline_method, "baseline_method")

  ## ---- headers and channel availability ---------------------------------
  if (!is.character(headers) || is.null(names(headers))) {
    stop(
      "`headers` must be a *named* character vector, e.g. c(cad = \"UV#1_CAD_1_0\").",
      call. = FALSE
    )
  }
  if (!detector %in% names(headers)) {
    stop(
      sprintf(
        "`detector = \"%s\"` has no entry in `headers`. Available: %s.",
        detector,
        paste(sprintf("\"%s\"", names(headers)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## ---- numeric ranges ----------------------------------------------------
  check_number(time_min, "time_min", min = 0)
  check_number(time_max, "time_max", min = 0)
  if (time_min >= time_max) {
    stop(
      sprintf(
        "`time_min` (%s) must be strictly less than `time_max` (%s).",
        format(time_min),
        format(time_max)
      ),
      call. = FALSE
    )
  }
  check_number(shift, "shift")
  check_number(frequency, "frequency", min = 0, inclusive_min = FALSE)
  check_number(resample, "resample", min = 1)
  check_number(
    fourier_components,
    "fourier_components",
    min = 0,
    max = 1,
    inclusive_min = FALSE
  )
  check_number(intensity_floor, "intensity_floor", min = 0)
  check_number(k2, "k2", min = 0, inclusive_min = FALSE)
  check_number(k4, "k4", min = 0, inclusive_min = FALSE)
  check_number(sigma, "sigma", min = 0)
  check_number(max_iter, "max_iter", min = 1)
  check_number(min_peak_height, "min_peak_height", min = 0, max = 1)
  check_number(min_area, "min_area", min = 0, max = 1)
  check_number(min_area_absolute, "min_area_absolute", min = 0)
  check_number(min_intensity, "min_intensity", min = 0)
  check_number(intensity_threshold, "intensity_threshold", min = 0, max = 1)

  if (!is.logical(improve_signal) || length(improve_signal) != 1L) {
    stop("`improve_signal` must be a single TRUE or FALSE.", call. = FALSE)
  }

  ## ---- unit resolution ---------------------------------------------------
  ## sd_max and smoothing_width are consumed in grid points. Expressing them in
  ## minutes and converting here keeps their physical meaning stable when
  ## frequency or resample change.
  step <- grid_step(frequency = frequency, resample = resample)

  if (!is.null(smoothing_width_minutes)) {
    check_number(
      smoothing_width_minutes,
      "smoothing_width_minutes",
      min = 0,
      inclusive_min = FALSE
    )
    smoothing_width <- max(1L, as.integer(round(smoothing_width_minutes / step)))
  }
  check_number(smoothing_width, "smoothing_width", min = 1)

  if (!is.null(sd_max_minutes)) {
    check_number(
      sd_max_minutes,
      "sd_max_minutes",
      min = 0,
      inclusive_min = FALSE
    )
    sd_max <- sd_max_minutes / step
  }
  check_number(sd_max, "sd_max", min = 0, inclusive_min = FALSE)

  ## ---- deprecated --------------------------------------------------------
  if (!is.null(noise_threshold)) {
    warning(
      "`noise_threshold` is deprecated and has no effect: it was only forwarded ",
      "to a purity estimator that is disabled for single-channel input. ",
      "Use `min_peak_height` to control detection sensitivity.",
      call. = FALSE
    )
  }

  structure(
    list(
      detector = detector,
      chromatogram = chromatogram,
      headers = headers,
      time_min = time_min,
      time_max = time_max,
      shift = shift,
      frequency = frequency,
      resample = resample,
      grid_step_minutes = step,
      improve_signal = improve_signal,
      fourier_components = fourier_components,
      intensity_floor = intensity_floor,
      smoothing_width = smoothing_width,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      baseline_method = baseline_method,
      fit = fit,
      sd_max = sd_max,
      max_iter = max_iter,
      min_peak_height = min_peak_height,
      min_area = min_area,
      min_area_absolute = min_area_absolute,
      min_intensity = min_intensity,
      intensity_threshold = intensity_threshold
    ),
    class = "cascade_params"
  )
}

#' Match a value against its allowed set
#'
#' @param value The supplied value.
#' @param name Parameter name, used in the error message and to look up the
#'   allowed set in [cascade_choices].
#'
#' @return The matched value.
#'
#' @keywords internal
#'
#' @examples NULL
match_choice <- function(value, name) {
  choices <- cascade_choices[[name]]
  if (length(value) != 1L || !is.character(value) || !value %in% choices) {
    stop(
      sprintf(
        "`%s` must be one of %s, not %s.",
        name,
        paste(sprintf("\"%s\"", choices), collapse = ", "),
        if (is.character(value) && length(value) == 1L) {
          sprintf("\"%s\"", value)
        } else {
          paste0("`", paste(format(value), collapse = ", "), "`")
        }
      ),
      call. = FALSE
    )
  }
  value
}

#' Check a scalar numeric parameter
#'
#' @param value The supplied value.
#' @param name Parameter name, used in the error message.
#' @param min,max Optional bounds.
#' @param inclusive_min,inclusive_max Are the bounds inclusive?
#'
#' @return `NULL`, invisibly. Called for the error side effect.
#'
#' @keywords internal
#'
#' @examples NULL
check_number <- function(
  value,
  name,
  min = -Inf,
  max = Inf,
  inclusive_min = TRUE,
  inclusive_max = TRUE
) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop(
      sprintf(
        "`%s` must be a single non-missing number, not %s.",
        name,
        paste(format(value), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  too_low <- if (inclusive_min) value < min else value <= min
  too_high <- if (inclusive_max) value > max else value >= max
  if (too_low || too_high) {
    stop(
      sprintf(
        "`%s` must be in %s%s, %s%s, but is %s.",
        name,
        if (inclusive_min) "[" else "(",
        format(min),
        format(max),
        if (inclusive_max) "]" else ")",
        format(value)
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Print a cascade_params object
#'
#' @param x A `cascade_params` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
#'
#' @examples NULL
print.cascade_params <- function(x, ...) {
  cat("<cascade_params>\n")
  for (nm in names(x)) {
    v <- x[[nm]]
    cat(sprintf(
      "  %-24s %s\n",
      nm,
      paste(format(v), collapse = ", ")
    ))
  }
  invisible(x)
}

#' Coerce parameters to a flat data frame
#'
#' @description
#' Used to write the provenance sidecar next to the exported tables.
#'
#' @param params A `cascade_params` object.
#'
#' @return A two-column data frame: `parameter`, `value`.
#'
#' @keywords internal
#'
#' @examples NULL
params_as_table <- function(params) {
  flat <- purrr::imap(params, function(v, nm) {
    if (length(v) > 1L) {
      ## e.g. `headers`, which is a named vector: keep the names visible.
      value <- paste(
        if (!is.null(names(v))) paste0(names(v), "=", v) else format(v),
        collapse = "; "
      )
    } else {
      value <- paste(format(v), collapse = "")
    }
    data.frame(parameter = nm, value = value, stringsAsFactors = FALSE)
  })
  do.call(rbind, unname(flat))
}
