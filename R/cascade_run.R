#' Prepare a CASCADE run
#'
#' @export
#'
#' @include cascade_params.R
#' @include load_chromatograms.R
#' @include load_features.R
#' @include load_ms_data.R
#' @include load_name.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#'
#' @description
#' Runs the front half of the pipeline once: load the inputs, prepare the feature
#' table, condition the detector trace, detect peaks, and match features to them.
#'
#' This used to be copy-pasted into each of the four entry points, which is how
#' they drifted apart — most visibly, `process_compare_peaks()` forwarded only six
#' of the signal-processing parameters, so a peak set validated in
#' `check_peaks_integration()` could not be reproduced on export. Both now call
#' this function, so tuning transfers by construction.
#'
#' Pass the result to [check_peaks_integration()] and [process_compare_peaks()] to
#' guarantee they operate on exactly the same peaks, and to skip recomputing this
#' stage:
#'
#' ```r
#' run <- cascade_run(file, features, detector = "cad", sigma = 0.08)
#' check_peaks_integration(run)   # look at it
#' process_compare_peaks(run)     # export exactly what you looked at
#' ```
#'
#' @param file Path to the mzML.
#' @param features Path to the MZmine feature table.
#' @param params A `cascade_params` object. Build one with [cascade_params()];
#'   the default is the package defaults.
#' @param show_example Use the bundled example data instead of `file`/`features`?
#' @param example_polarity Polarity of the bundled example chromatograms.
#' @param subsample Optional integer: keep this many randomly-chosen features
#'   (seeded, so it is reproducible). Used to keep the examples fast.
#' @param shapes Compute the per-peak normalised chromatogram shapes needed for
#'   MS comparison? `FALSE` skips work that a purely visual QC pass never reads.
#'   [process_compare_peaks()] requires `TRUE`.
#' @param load_ms Load the raw MS1 scans? Only [process_compare_peaks()] needs
#'   them, and they are the most expensive input to open.
#'
#' @return An object of class `cascade_run`: a named list with
#'   `params`, `name`, `chromatograms`, `peaks`, `ms_data` and `shapes`.
#'
#' @examples
#' \dontrun{
#' run <- cascade_run(show_example = TRUE)
#' }
cascade_run <- function(
  file = NULL,
  features = NULL,
  params = cascade_params(),
  show_example = FALSE,
  example_polarity = "pos",
  subsample = NULL,
  shapes = TRUE,
  load_ms = FALSE
) {
  if (!inherits(params, "cascade_params")) {
    stop(
      "`params` must be a `cascade_params` object, as built by `cascade_params()`.",
      call. = FALSE
    )
  }

  ## ---- 1. Load raw inputs ------------------------------------------------
  ## The MS1 scans are opened only when asked for: `onDisk` keeps an index in
  ## RAM and reads spectra on demand, but even opening the file is slow.
  ms_data <- NULL
  if (load_ms) {
    message("loading MS data")
    ms_data <- file |>
      load_ms_data(show_example = show_example)
  }

  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(
      show_example = show_example,
      headers = params$headers,
      example_polarity = example_polarity
    )

  ## Sample name = basename(file). Keys this sample's columns in the
  ## multi-sample MZmine export ("datafile:<name>:rt_range:min", ...).
  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  if (!is.null(subsample)) {
    message("selecting ", subsample, " random features")
    set.seed(42)
    feature_table <- feature_table |>
      tidytable::slice_sample(n = subsample)
  }

  message("preparing features")
  df_features <- feature_table |>
    prepare_features(min_intensity = params$min_intensity, name = name)

  ## ---- 2. Condition the chosen detector trace ----------------------------
  ## `params$detector` was validated against `names(params$headers)` in
  ## cascade_params(), so this subset cannot come back empty the way the old
  ## unchecked `switch()` could.
  channel <- params$headers[params$detector]
  trace <- chromatograms_all[names(channel)]

  if (length(trace) == 0L || is.null(trace[[1]])) {
    stop(
      sprintf(
        "Detector channel \"%s\" (id \"%s\") was not found in the file. Check `headers`.",
        params$detector,
        channel
      ),
      call. = FALSE
    )
  }

  message("preprocessing chromatograms")
  chromatograms_list <- preprocess_chromatograms(
    detector = params$detector,
    name = name,
    list = trace,
    shift = params$shift,
    fourier_components = params$fourier_components,
    time_min = params$time_min,
    time_max = params$time_max,
    frequency = params$frequency,
    resample = params$resample,
    intensity_floor = params$intensity_floor,
    k2 = params$k2,
    k4 = params$k4,
    sigma = params$sigma,
    smoothing_width = params$smoothing_width,
    baseline_method = params$baseline_method,
    improve_signal = params$improve_signal
  )

  ## ---- 3. Detect peaks and match them to features ------------------------
  ## Two views of the same trace: `_long` (row-bound, with interval columns)
  ## drives the joins, the wide form feeds the NLS fitter.
  message("preprocessing peaks")
  peaks <- preprocess_peaks(
    detector = params$detector,
    df_features = df_features,
    df_long = pick_chromatogram(
      chromatograms_list,
      params$chromatogram,
      long = TRUE
    ) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    df_xy = pick_chromatogram(
      chromatograms_list,
      params$chromatogram,
      long = FALSE
    ),
    min_area = params$min_area,
    min_area_absolute = params$min_area_absolute,
    shift = params$shift,
    name = name,
    sd_max = params$sd_max,
    max_iter = params$max_iter,
    fit = params$fit,
    min_peak_height = params$min_peak_height,
    intensity_threshold = params$intensity_threshold,
    shapes = shapes
  )

  structure(
    list(
      params = params,
      name = name,
      chromatograms = chromatograms_list,
      peaks = peaks,
      ms_data = ms_data,
      shapes = shapes
    ),
    class = "cascade_run"
  )
}

#' Select one of the three processed trace versions
#'
#' @description
#' Replaces the `switch(chromatogram, "original" = ..., ...)` blocks that were
#' repeated twice per entry point. The value has already been validated by
#' [cascade_params()], so an unmatched key here is an internal error rather than
#' a user mistake.
#'
#' @param chromatograms_list The list returned by [preprocess_chromatograms()].
#' @param which One of `"original"`, `"improved"`, `"baselined"`.
#' @param long Return the long form (`TRUE`) or the wide single-sample data
#'   frame (`FALSE`)?
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @examples NULL
pick_chromatogram <- function(chromatograms_list, which, long = TRUE) {
  key <- paste0(
    "chromatograms_",
    which,
    if (long) "_long" else ""
  )
  out <- chromatograms_list[[key]]
  if (is.null(out)) {
    stop(sprintf("Internal error: no chromatogram named `%s`.", key), call. = FALSE)
  }
  if (long) out else out[[1]]
}

#' Print a cascade_run object
#'
#' @param x A `cascade_run` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
#'
#' @examples NULL
print.cascade_run <- function(x, ...) {
  n_peaks <- length(x$peaks$list_df_features_with_peaks_long)
  n_out <- nrow(x$peaks$df_features_without_peaks)
  cat("<cascade_run>\n")
  cat("  sample        ", x$name, "\n")
  cat("  detector      ", x$params$detector, "\n")
  cat("  trace         ", x$params$chromatogram, "\n")
  cat("  matched units ", n_peaks, "(feature x peak)\n")
  cat("  unmatched     ", if (is.null(n_out)) 0L else n_out, "features\n")
  cat("  shapes        ", x$shapes, "\n")
  cat("  MS1 loaded    ", !is.null(x$ms_data), "\n")
  invisible(x)
}
