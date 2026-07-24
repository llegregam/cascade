# cascade

## Parameter uniformity refactor

### Breaking / behaviour-changing

Six parameters had a different default depending on which layer you called. They
now resolve to a single value, taken from the entry-point defaults. Results can
move **only if you were calling the helper functions directly** — calls through
the exported entry points are unaffected, since those already used these values.

| Parameter | Old (helper layer) | New | Affected functions |
|---|---|---|---|
| `frequency` | `2` | `1` | `preprocess_chromatograms()`, `improve_signal()`, `improve_signals_progress()` |
| `time_min` | `0` | `0.5` | same |
| `time_max` | `Inf` | `32.5` | same |
| `shift` | `0` | `0.05` | `preprocess_chromatograms()`, `preprocess_peaks()` |
| `min_area` | `0` | `0.005` | `preprocess_peaks()` |
| `max_iter` | `100` | `1000` | `get_peaks()` (was `max.iter`) |

- `check_chromatograms_alignment()`: `pda_shift` now defaults to the same value
  as `cad_shift` (`0.05`). It previously defaulted to `0.1`, which was arbitrary
  rather than an instrument property.
- `plot_peak_detection()`: the "detected minimum (end)" marker took its height
  from the peak's **start** retention time. It now uses `rt_max`. Marker
  placement only — no effect on integration.
- `get_peaks()`: a call to `remove_bad_peaks()` was missing its required `n`
  argument, which errored at peak detection. Restored.
- `process_compare_peaks()` now returns its results instead of `NULL`.

### Dependencies

- `tima` moved from `Imports` to `Suggests`. It is used only by the annotation
  path (`generate_ids()`, `generate_tables()`), to download LOTUS from Zenodo.
  The peak-integration and alignment workflow does not need it, and the package
  now loads without it. Both call sites go through the new guarded
  `get_lotus()`, which explains how to install `tima` if you do need it.
- `minpack.lm` added to `Imports`, and `nlsLM()` is now called as
  `minpack.lm::nlsLM()`. It was previously called unqualified and undeclared.
  Because every fit is wrapped in `try()`, a missing `minpack.lm` did not error
  — it silently fell back to the starting estimates, so peaks were reported but
  never actually fitted.

### Deprecated

- `process_compare_peaks(type =)` → `chromatogram =`, matching
  `check_peaks_integration()`. The old name still works and warns.
- `get_peaks(sd.max =, max.iter =)` → `sd_max =`, `max_iter =`. A deliberate
  divergence from upstream `chromatographR`, noted in the file header. Old names
  warn.
- `signal_sharpening(Smoothing_width =, Baseline_adjust =)` →
  `smoothing_width =`, `baseline_adjust =`.
- `noise_threshold` is deprecated everywhere and warns if set. It never had any
  effect: it reached only a purity estimator that is disabled for
  single-channel input, which is always the case in this pipeline. Use the new
  `min_peak_height` to control detection sensitivity.
- `extract_ms_progress(nrows =)` is accepted and ignored; it was never read.

### New

- `cascade_defaults`: one object holding every shared default. Function
  signatures now read from it, so a value is defined exactly once and cannot
  drift again.
- `cascade_params()`: validates every argument at the entry point and reports
  the offending name and value. A mistyped `chromatogram` or `detector` used to
  produce a `NULL` from an unchecked `switch()` and fail somewhere unrelated.
- `cascade_run()` (exported): runs the load → prepare → condition → detect front
  half once. Pass the result to both `check_peaks_integration()` and
  `process_compare_peaks()` so the peaks you inspect are the peaks you export.
- `process_compare_peaks()` now exposes the full signal-processing and
  peak-fitting surface. It previously forwarded only six parameters and ran the
  rest at helper defaults, so a peak set tuned in `check_peaks_integration()`
  could not be reproduced on export.
- `min_peak_height`: detection sensitivity as a fraction of the trace maximum.
  This is `amp_thresh` inside `get_peaks()`, which no exported function reached.
- `sd_max_minutes` and `smoothing_width_minutes`: the same quantities in
  minutes. `sd_max` and `smoothing_width` are in **grid points**, so their
  physical meaning silently changed whenever `frequency` or `resample` did.
- `min_area_absolute`: an absolute area cutoff, applied alongside the existing
  relative `min_area`.
- `process_compare_peaks()` writes a `<name>_params_<detector>.tsv` sidecar
  recording every resolved parameter.
- `check_peaks_integration()` attaches its peak table to the returned figure as
  the `peaks` attribute.

# cascade 0.0.9001

- Adapt to `tima` updates
- Added helper functions to visualize
  [TIMA](https://taxonomicallyinformedannotation.github.io/tima/) results
- Switched documentation from `pkgdown` to `altdoc`
- Updated minimal R version to `4.4.0` (and related Bioconductor dependencies)

# cascade 0.0.9000

- Initial version.
