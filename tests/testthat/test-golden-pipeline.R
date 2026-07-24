library(testthat)

## Golden-output regression tests.
##
## Purpose: make every future change to the pipeline provably behaviour
## preserving. These snapshot the numbers the pipeline produces on the bundled
## example data. If a refactor is meant to be invisible, these stay green; if a
## default is deliberately reconciled, the snapshot moves and that change must be
## recorded in NEWS.md.
##
## Determinism: the example path subsamples features with a fixed seed and the
## pipeline uses no other randomness, so repeated runs agree.

## Helpers must be defined before the test_that() calls that use them:
## test_that() runs as the file is sourced, so anything declared below it does
## not exist yet.
temp_export_dir <- function() {
  d <- file.path(tempdir(), paste0("cascade-test-", as.integer(Sys.time())))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

peak_fingerprint <- function(run) {
  ## A compact, order-stable summary of what peak detection found.
  ##
  ## Deliberately PEAK-level, not (feature x peak)-level. The long table repeats
  ## each peak once per feature that overlaps it, so snapshotting it stores the
  ## same ~50 peaks thousands of times -- on the example data that is a 400 kB
  ## base64 blob that churns on every run. The peaks are what peak detection
  ## actually produces; the expansion is a join artifact, and n_matches below
  ## still catches any change to the join.
  ##
  ## Rounded because NLS fits are not bit-reproducible across BLAS builds.
  df <- run$peaks$list_df_features_with_peaks_long |>
    tidytable::bind_rows()
  if (nrow(df) == 0L) {
    return(list(n_matches = 0L, n_peaks = 0L))
  }
  peaks <- df |>
    tidytable::distinct(peak_id, rt_min, rt_apex, rt_max, integral)
  peaks <- peaks[order(peaks$peak_id), ]
  list(
    n_matches = nrow(df),
    n_peaks = nrow(peaks),
    peak_id = as.integer(peaks$peak_id),
    rt_min = round(as.numeric(peaks$rt_min), 4),
    rt_apex = round(as.numeric(peaks$rt_apex), 4),
    rt_max = round(as.numeric(peaks$rt_max), 4),
    integral = signif(as.numeric(peaks$integral), 6)
  )
}

test_that(desc = "example run peak detection is stable", code = {
  run <- cascade_run(show_example = TRUE, shapes = FALSE)
  expect_s3_class(run, "cascade_run")
  expect_snapshot_value(peak_fingerprint(run), style = "serialize")
})

test_that(desc = "example run peak detection is stable off-default", code = {
  ## A non-default DSP setting, so the snapshot is sensitive to the default
  ## reconciliations rather than only to the code path.
  run <- cascade_run(
    show_example = TRUE,
    shapes = FALSE,
    params = cascade_params(sigma = 0.08, fourier_components = 0.02)
  )
  expect_snapshot_value(peak_fingerprint(run), style = "serialize")
})

test_that(desc = "check_peaks_integration returns an inspectable figure", code = {
  p <- check_peaks_integration(show_example = TRUE)
  ## Bug fixed by the refactor: the QC figure used to discard its peak table.
  expect_false(is.null(attr(p, "peaks")))
  expect_s3_class(attr(p, "params"), "cascade_params")
})

test_that(desc = "process_compare_peaks exports and returns its tables", code = {
  dir <- temp_export_dir()
  res <- process_compare_peaks(show_example = TRUE, export_dir = dir)

  ## Bug fixed by the refactor: this function used to end on fwrite() and so
  ## returned NULL, making `x <- process_compare_peaks(...)` an empty variable.
  expect_false(is.null(res))
  expect_true(all(
    c("informed", "not_informed", "files", "params") %in% names(res)
  ))
  expect_true(all(file.exists(res$files)))

  ## Both tables must share one schema so they can be row-bound downstream.
  expect_identical(names(res$informed), names(res$not_informed))

  ## The provenance sidecar must round-trip every resolved parameter.
  sidecar <- tidytable::fread(res$files[["params"]])
  expect_equal(nrow(sidecar), length(res$params))

  expect_snapshot_value(
    list(
      cols = names(res$informed),
      n_informed = nrow(res$informed),
      n_not_informed = nrow(res$not_informed)
    ),
    style = "serialize"
  )
})

test_that(desc = "QC and export agree on the same peaks", code = {
  ## The point of the refactor. Before it, process_compare_peaks() forwarded
  ## only six of the signal-processing parameters and ran the rest at helper
  ## defaults, so a peak set validated with a non-default sigma could not be
  ## reproduced on export. This test fails on the pre-refactor code.
  params <- cascade_params(sigma = 0.08, fourier_components = 0.02)

  qc_run <- cascade_run(
    show_example = TRUE,
    params = params,
    shapes = FALSE
  )
  export_run <- cascade_run(
    show_example = TRUE,
    params = params,
    shapes = TRUE,
    load_ms = TRUE
  )

  expect_equal(
    peak_fingerprint(qc_run),
    peak_fingerprint(export_run)
  )
})

test_that(desc = "a QC-only run is refused by the scorer, with a useful message", code = {
  run <- cascade_run(show_example = TRUE, shapes = FALSE)
  expect_error(
    process_compare_peaks(run),
    "shapes = FALSE"
  )
})
