library(testthat)

## These tests cover the parameter layer only. They need no MS data, no example
## files and no heavy dependencies, so they run fast and catch the class of bug
## the refactor was written to prevent: a default drifting between layers, or a
## bad argument slipping through to fail somewhere unrelated.

test_that(desc = "defaults are defined exactly once and reused", code = {
  ## Every function that takes one of these must inherit it from
  ## cascade_defaults rather than carrying its own literal. If someone
  ## reintroduces a literal, the value can drift again -- which is exactly what
  ## this refactor removed.
  shared <- c(
    "frequency",
    "resample",
    "time_min",
    "time_max",
    "shift",
    "fourier_components",
    "intensity_floor",
    "smoothing_width",
    "k2",
    "k4",
    "sigma",
    "baseline_method",
    "fit",
    "sd_max",
    "max_iter",
    "min_area",
    "min_intensity",
    "intensity_threshold"
  )
  expect_true(all(shared %in% names(cascade_defaults)))

  targets <- list(
    check_peaks_integration = check_peaks_integration,
    process_compare_peaks = process_compare_peaks,
    preprocess_chromatograms = preprocess_chromatograms,
    preprocess_peaks = preprocess_peaks,
    improve_signal = improve_signal,
    improve_signals_progress = improve_signals_progress,
    peaks_progress = peaks_progress,
    join_peaks = join_peaks
  )

  for (fname in names(targets)) {
    fmls <- formals(targets[[fname]])
    for (p in intersect(names(fmls), shared)) {
      expr <- fmls[[p]]
      ## Must literally be `cascade_defaults$<p>`, not a bare value.
      expect_true(
        is.call(expr) && identical(deparse(expr), paste0("cascade_defaults$", p)),
        info = sprintf(
          "%s() declares a literal default for `%s` instead of cascade_defaults$%s",
          fname,
          p,
          p
        )
      )
    }
  }
})

test_that(desc = "the six historically drifted defaults now agree", code = {
  ## Before the refactor: frequency was 1 up top and 2 underneath, time_max was
  ## 32.5 up top and Inf underneath, max_iter was 1000 up top and 100 in
  ## get_peaks(). The effective default depended on which layer you entered at.
  drifted <- c("frequency", "time_min", "time_max", "shift", "min_area", "max_iter")
  layers <- list(
    check_peaks_integration,
    preprocess_chromatograms,
    preprocess_peaks,
    improve_signal,
    improve_signals_progress,
    get_peaks
  )
  for (p in drifted) {
    values <- lapply(layers, function(f) {
      fmls <- formals(f)
      if (!p %in% names(fmls)) {
        return(NULL)
      }
      eval(fmls[[p]], envir = asNamespace("cascade"))
    })
    values <- Filter(Negate(is.null), values)
    expect_true(
      length(unique(values)) == 1L,
      info = sprintf("`%s` still differs between layers: %s", p, toString(values))
    )
  }
})

test_that(desc = "enumerated parameters are validated by name", code = {
  expect_error(cascade_params(detector = "kad"), "`detector` must be one of")
  expect_error(
    cascade_params(chromatogram = "baslined"),
    "`chromatogram` must be one of"
  )
  expect_error(cascade_params(fit = "egc"), "`fit` must be one of")
  expect_error(
    cascade_params(baseline_method = "nope"),
    "`baseline_method` must be one of"
  )
  ## The mistyped value must appear in the message, so the user can see it.
  expect_error(cascade_params(detector = "kad"), "kad")
})

test_that(desc = "numeric parameters are range checked", code = {
  expect_error(cascade_params(time_min = 10, time_max = 5), "strictly less than")
  expect_error(cascade_params(fourier_components = 0), "fourier_components")
  expect_error(cascade_params(fourier_components = 1.5), "fourier_components")
  expect_error(cascade_params(frequency = 0), "frequency")
  expect_error(cascade_params(resample = 0.5), "resample")
  expect_error(cascade_params(min_area = 2), "min_area")
  expect_error(cascade_params(intensity_threshold = -1), "intensity_threshold")
  expect_error(cascade_params(sigma = "loud"), "single non-missing number")
  expect_error(cascade_params(improve_signal = "yes"), "single TRUE or FALSE")
})

test_that(desc = "detector must exist in headers", code = {
  expect_error(
    cascade_params(detector = "pda", headers = c(cad = "UV#1_CAD_1_0")),
    "has no entry in `headers`"
  )
  expect_error(
    cascade_params(headers = c("UV#1_CAD_1_0")),
    "must be a \\*named\\* character vector"
  )
})

test_that(desc = "minute-denominated widths convert to grid points", code = {
  ## 1 Hz, resample 1 -> one point per second -> 0.5 min = 30 points.
  p <- cascade_params(frequency = 1, resample = 1, sd_max_minutes = 0.5)
  expect_equal(p$sd_max, 30)

  ## Doubling the grid density must double the point count for the same
  ## physical width. This is the trap the _minutes variants remove: a bare
  ## sd_max = 30 would silently mean 0.25 min here.
  p2 <- cascade_params(frequency = 1, resample = 2, sd_max_minutes = 0.5)
  expect_equal(p2$sd_max, 60)

  p3 <- cascade_params(
    frequency = 1,
    resample = 1,
    smoothing_width_minutes = 10 / 60
  )
  expect_equal(p3$smoothing_width, 10)
})

test_that(desc = "deprecated arguments warn but still work", code = {
  expect_warning(
    p <- cascade_params(noise_threshold = 0.001),
    "`noise_threshold` is deprecated"
  )
  expect_false("noise_threshold" %in% names(p))
})

test_that(desc = "params round-trip to a provenance table", code = {
  p <- cascade_params()
  tbl <- params_as_table(p)
  expect_true(all(c("parameter", "value") %in% names(tbl)))
  expect_equal(nrow(tbl), length(p))
  ## `headers` is a named vector: it must survive as readable text.
  expect_true(any(grepl("cad=", tbl$value, fixed = TRUE)))
})
