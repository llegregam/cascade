# CASCADE parameter-uniformity refactor — actionable plan

**Status:** IMPLEMENTED (2026-07-24). See §7 for what landed and what still needs
running in a full environment.
**Written:** 2026-07-24
**Scope:** the parameter surface of the four exported entry points and the helper
chain beneath them. No change to the scientific method.

---

## 0. How to use this file

The phases below were executed. This file is now the record of *why* each change
was made — the evidence in §2 is the justification for the design in §3, and §7
records the outcome. Keep it: if a future change reintroduces a literal default
or a second parameter surface, §2 is the argument against it.

Companion documents:
- `notes/peak_pipeline.md` — line-by-line walkthrough (describes the **pre**-refactor
  code; still accurate for the science, stale for the call structure).
- `notes/pipeline_concepts.md` — the science behind the pipeline.

---

## 1. Codebase orientation

Four exported entry points share one helper chain:

| Entry point | File | Params | Output |
|---|---|---|---|
| `check_chromatograms_alignment` | `R/check_chromatograms_alignment.R` | 22 | plotly figure (detector alignment QC) |
| `check_peaks_integration` | `R/check_peaks_integration.R` | 26 | plotly figure (peak-detection QC) |
| `process_compare_peaks` | `R/process_compare_peaks.R` | **15** | 2 TSVs on disk |
| `generate_pseudochromatograms` | `R/generate_pseudochromatograms.R` | 24 | plots + tables |

Shared helper chain (roughly in call order):

```
load_chromatograms / load_features / load_name / load_ms_data
        └─ prepare_features
        └─ preprocess_chromatograms
                 └─ improve_signals_progress → improve_signal
                          └─ filter_fft
                          └─ signal_sharpening → second_der → deriv / middle_pts
                 └─ baseline_chromatogram
        └─ preprocess_peaks
                 └─ peaks_progress → get_peaks        (vendored from chromatographR)
                 └─ join_peaks
                 └─ normalize_chromato / prepare_peaks / prepare_rt / prepare_mz
        └─ extract_ms_progress → transform_ms         (process_compare_peaks only)
        └─ extract_ms_peak → compare_peaks            (process_compare_peaks only)
```

---

## 2. Evidence — what is actually wrong

### 2.1 The same parameter has a different default at every layer

| Parameter | Entry points | Mid layer | Leaf |
|---|---|---|---|
| `frequency` | `1` — `check_peaks_integration.R:71` | `2` — `preprocess_chromatograms.R:36`, `improve_signals_progress.R:25` | `2` — `improve_signal.R:27` |
| `time_min` | `0.5` — `check_peaks_integration.R:69` | `0` — `preprocess_chromatograms.R:41` | `0` — `improve_signal.R:29` |
| `time_max` | `32.5` — `check_peaks_integration.R:70` | `Inf` — `preprocess_chromatograms.R:42` | `Inf` — `improve_signal.R:30` |
| `shift` | `0.05` — `check_peaks_integration.R:66` | `0` — `preprocess_chromatograms.R:40`, `preprocess_peaks.R:34` | — |
| `min_area` | `0.005` — `check_peaks_integration.R:64` | `0` — `preprocess_peaks.R:35` | — |
| `max_iter` | `1000` — `check_peaks_integration.R:80` | `1000` — `preprocess_peaks.R:37` | **`100`** — `get_peaks.R:30` (`max.iter`) |

Nothing enforces agreement. The effective default depends on which function you
call directly — a trap for anyone using the helpers outside the entry points.

### 2.2 Four entry points, four parameter surfaces

`process_compare_peaks` — the only one that writes results — exposes 15
parameters. It forwards to `preprocess_chromatograms` only
`shift, fourier_components, time_min, time_max, frequency, resample`
(`process_compare_peaks.R`), and to `preprocess_peaks` only
`min_area, shift, name, detector`. Everything else silently runs at helper
defaults:

```
improve_signal = TRUE, intensity_floor = 0.001, smoothing_width = 8,
k2 = 250, k4 = 1250000, sigma = 0.05, baseline_method = "peakDetection",
fit = "egh", sd_max = 50, max_iter = 1000, intensity_threshold = 0.1
```

**Consequence:** a peak set validated in `check_peaks_integration` with non-default
DSP settings cannot be reproduced by `process_compare_peaks`. The QC tool and the
production tool are not the same pipeline.

### 2.3 Naming divergence for identical concepts

| Concept | Names in use |
|---|---|
| which trace version | `chromatogram` (check) vs `type` (compare) |
| detector time offset | `shift` vs `cad_shift` / `pda_shift` |
| fit controls | `sd_max` / `max_iter` vs `sd.max` / `max.iter` (`get_peaks`) |

### 2.4 Implicit and mixed units

| Parameter | Unit | Problem |
|---|---|---|
| `time_min`, `time_max`, `shift` | minutes | fine |
| `prepare_rt` output | **seconds** | silent ×60 at the MSnbase boundary |
| `frequency` | Hz | fine |
| `sd_max`, `smoothing_width` | **grid points** | physical meaning changes when `frequency`/`resample` change |
| `min_area` | **relative fraction** | name implies absolute |
| `fourier_components`, `intensity_threshold` | fraction | fine |
| `min_intensity` | absolute counts | fine |

`sd_max` is the worst offender: with `fit = "raw"`, `get_peaks` reports
`sd = end - start` in grid points (`get_peaks.R`, `fitpk_raw`), so `sd_max = 50`
silently discards every peak wider than ~50 s at 1 Hz.

### 2.5 Dead and misleading parameters

- **`noise_threshold` is inert.** Only forwarded to `get_purity()`, which is
  force-disabled when `ncol(x) == 1` — always true, since `peaks_progress` builds
  a single-column matrix.
- **`intensity_threshold` is inert in `check_peaks_integration`.** It drives
  `list_chromato_peaks`, which the QC plot never reads. Computed and discarded.
- **`nrows` is unused** in `extract_ms_progress`'s inner map.
- **Detection sensitivity is unreachable.** `find_peaks`'s `amp_thresh`,
  `slope_thresh`, `smooth_type`, `smooth_window` are not plumbed to any public API.
  Users are forced to tune sensitivity indirectly via `fourier_components`.

### 2.6 No argument validation

No `match.arg` or `stopifnot` in any entry point. `chromatogram = "baslined"`
makes `switch()` return `NULL`, which fails several calls deeper with an
unrelated message. Same for a mistyped `detector`.

### 2.7 Concrete bugs found while reading

| # | Location | Bug |
|---|---|---|
| 1 | `check_peaks_integration.R:151-175` | `detector` is never passed to `preprocess_peaks`, so it always logs `"preprocessing cad peaks"` regardless of the real detector |
| 2 | `plot_peak_detection.R:45` | end-marker height uses `fun(rt_min)` instead of `fun(rt_max)` — copy-paste slip from the start-marker block |
| 3 | `check_peaks_integration.R:1` | roxygen title reads `"Check chromatograms alignment"` — leftover from `check_chromatograms_alignment` |
| 4 | `process_compare_peaks.R` (end) | function ends on `fwrite`, so it returns `NULL`; `x <- process_compare_peaks(...)` yields an empty variable |
| 5 | all four entry points | local variables named `switch` and `list` shadow the base-R functions |
| 6 | `process_compare_peaks.R` roxygen | missing `@include load_ms_data.R` although it calls `load_ms_data` (cosmetic: affects Collate order only) |

---

## 3. Target design

**One definition per parameter, one name per concept, one unit per quantity, and
one code path shared by QC and production.**

Concretely:

1. Defaults live in exactly one object (`cascade_defaults`), referenced by every
   signature. Drift becomes impossible rather than merely discouraged.
2. Every entry point validates its arguments at the door.
3. The load → prepare → preprocess → detect front half exists once, as an internal
   function, not copy-pasted four times.
4. `process_compare_peaks` accepts the same surface as `check_peaks_integration`.
5. Parameters are expressed in physical units (minutes), converted to grid points
   internally.

---

## 4. Phases

### Phase 0 — Golden-output regression test (do this first)

**Goal:** make every later phase provably behaviour-preserving.

- [ ] `tests/testthat/test-golden-pipeline.R`: run `check_peaks_integration(show_example = TRUE)`
      and `process_compare_peaks(show_example = TRUE, export_dir = <tempdir>)`.
- [ ] Snapshot with `testthat::expect_snapshot_value()`:
      - the peak table from `preprocess_peaks` (`peak_id, rt_min, rt_apex, rt_max, integral`)
      - both exported TSVs (all 11 columns)
      - the number of features matched vs unmatched
- [ ] Run at defaults **and** at one non-default DSP setting, so the snapshot is
      sensitive to the reconciliations in Phase 1.

**Acceptance:** `devtools::test()` green; snapshots committed.

**Note:** `process_compare_peaks(show_example = TRUE)` subsamples 10 features with
`set.seed(42)`, so it is already deterministic. `check_peaks_integration` has no
sampling. Neither uses randomness elsewhere.

---

### Phase 1 — Single source of truth for defaults (non-breaking)

**1.1 Central defaults object**
- [ ] New `R/cascade_defaults.R` exporting an internal named list covering every
      shared parameter: `frequency, resample, time_min, time_max, shift,
      fourier_components, intensity_floor, k2, k4, sigma, smoothing_width,
      baseline_method, improve_signal, fit, sd_max, max_iter, min_area,
      min_intensity, intensity_threshold, headers`.
- [ ] Replace every literal default in every signature with
      `sigma = cascade_defaults$sigma`, etc. Applies to all four entry points plus
      `preprocess_chromatograms`, `preprocess_peaks`, `improve_signal`,
      `improve_signals_progress`, `peaks_progress`, `baseline_chromatogram`,
      `normalize_chromato`, `transform_ms`.
- [ ] `headers` currently repeats its default literal in 4+ files — collapse to
      `cascade_defaults$headers`.

**1.2 Reconcile the six drifted defaults** (table §2.1) to the entry-point values.
Notably `get_peaks`'s `max.iter: 100 → 1000` and `improve_signal`'s
`frequency: 2 → 1`. Record each in `NEWS.md`.

**1.3 Argument validation**
- [ ] New internal `validate_params()`, called as the first statement of each entry point:
      - `match.arg` on `detector`, `chromatogram`/`type`, `fit`, `baseline_method`
      - `stopifnot` on `0 < fourier_components <= 1`, `time_min < time_max`,
        `frequency > 0`, `resample >= 1`, `0 <= min_area <= 1`,
        `0 <= intensity_threshold <= 1`
      - check `detector` is a name of `headers`, and that the requested channel was
        actually found in the mzML (currently an empty subset fails much later)
- [ ] Error messages must name the offending parameter and its received value.

**1.4 Fix bugs 1–3 and 6** from §2.7. Leave 4 and 5 for Phase 2/3.

**Acceptance:** Phase 0 snapshots unchanged except where a reconciled default is
expected to move them — each such change explicitly listed in `NEWS.md`.
No signature changes; existing user scripts keep working.

---

### Phase 2 — Unify naming and returns (one deprecation cycle)

**2.1 Names**
- [ ] `process_compare_peaks`: `type` → `chromatogram`; accept `type` with
      `lifecycle::deprecate_warn()`.
- [ ] `get_peaks`: `sd.max` → `sd_max`, `max.iter` → `max_iter`, same treatment.
      (Vendored from chromatographR — note the divergence in the file header.)
- [ ] Consider `min_area` → `min_area_relative` (see 4.3), aliased.

**2.2 Useful return values**
- [ ] `process_compare_peaks` returns
      `invisible(list(informed = <df>, not_informed = <df>, files = <chr>, params = <list>))`.
      Fixes bug 4 and lets users inspect results without re-reading the TSVs.
- [ ] `check_peaks_integration` attaches the peak table as an attribute on the
      returned plotly object, so QC results are inspectable.

**2.3 Provenance sidecar**
- [ ] Write the resolved parameter list to `<name>_params_<detector>.json`
      alongside the TSVs. Turns "we integrated with cascade" into a reproducible
      method section.

**Acceptance:** deprecation warnings fire and are tested; snapshots unchanged.

---

### Phase 3 — Close the QC → production gap (highest value)

**3.1 Extract the shared front half**
- [ ] New internal `prepare_run()` containing load → `prepare_features` →
      `preprocess_chromatograms` → `preprocess_peaks`, returning a `cascade_run`
      object (chromatograms, peaks, features matched/unmatched, resolved params).
- [ ] Rewrite all four entry points to call it. This removes the fourfold
      copy-paste, including the `switch`/`list` shadowing (bug 5).

**3.2 Give `process_compare_peaks` the full surface**
- [ ] Forward every DSP and fit parameter through `prepare_run()`. Nearly free once
      3.1 lands. **This is the change that makes tuning transferable.**

**3.3 Accept a prepared run directly**
```r
run <- cascade_run(file, features, detector = "cad", sigma = 0.08, ...)
check_peaks_integration(run)   # look at it
process_compare_peaks(run)     # export exactly what you looked at
```
- [ ] Add `run` as the first argument of both, with the existing flat form kept.
      Guarantees by construction that validated peaks are exported peaks, and
      skips recomputing the front half.

**Acceptance:** snapshots unchanged. New test: a non-default `sigma` produces
identical peak tables via `check_peaks_integration` and `process_compare_peaks`
(this test **fails** on today's code — that is the point).

---

### Phase 4 — Make parameters mean what they say

- [ ] **4.1 Expose detection sensitivity.** Plumb `find_peaks`'s `amp_thresh` up as
      `min_peak_height` (fraction of max intensity). The missing knob.
- [ ] **4.2 Unit-aware parameters.** Accept `sd_max` and `smoothing_width` in
      **minutes**; convert to grid points internally from the resolved step
      `1 / (frequency * 60 * resample)`. Removes the silent-rescaling trap.
- [ ] **4.3 Split `min_area`** into `min_area_relative` (current) and
      `min_area_absolute`, both applied; `min_area` deprecated alias for the former.
- [ ] **4.4 Remove dead parameters.** Deprecate then remove `noise_threshold`.
      Either use or stop computing `list_chromato_peaks` in
      `check_peaks_integration`. Drop `nrows` from `extract_ms_progress`.

**Acceptance:** new parameters covered by tests; deprecations in `NEWS.md`.

---

## 5. Breaking-change policy

- Phases 0–1: no breaking changes.
- Phase 2: renames warn for one minor release, then error with a pointer.
- Phase 3: purely additive to the public API.
- Phase 4: `noise_threshold` removal is the only true break — it is inert, so no
  result changes.

Every default reconciliation gets its own `NEWS.md` bullet stating old value, new
value, and which results can move.

---

## 7. Outcome (2026-07-24)

### What landed

| Phase | Status | Where |
|---|---|---|
| 0 — golden-output test | Done | `tests/testthat/test-golden-pipeline.R` |
| 1.1 — central defaults | Done | `R/cascade_defaults.R` |
| 1.2 — reconcile drift | Done | all six, recorded in `NEWS.md` |
| 1.3 — validation | Done | `R/cascade_params.R`, `tests/testthat/test-params.R` |
| 1.4 — bugs 1,2,3,6 | Done | + a 7th found during the work (below) |
| 2.1 — naming | Done | `type`→`chromatogram`, `sd.max`→`sd_max`, `Smoothing_width`→`smoothing_width` |
| 2.2 — return values | Done | `process_compare_peaks()` returns a list; QC figure carries `attr(,"peaks")` |
| 2.3 — provenance sidecar | Done | `<name>_params_<detector>.tsv` (TSV not JSON, to avoid a new dependency) |
| 3.1 — shared front half | Done | `R/cascade_run.R` |
| 3.2 — full surface | Done | `process_compare_peaks()` went 15 → 31 parameters |
| 3.3 — run passthrough | Done | accepted via `run =` **or** positionally in `file` |
| 4.1 — `min_peak_height` | Done | plumbed to `find_peaks()`'s `amp_thresh` |
| 4.2 — unit-aware params | Done | added `sd_max_minutes`, `smoothing_width_minutes` |
| 4.3 — `min_area` split | Done | added `min_area_absolute` |
| 4.4 — dead parameters | Done | `noise_threshold` deprecated, `nrows` ignored, `shapes=FALSE` skips unused work |

### Deviations from the plan, and why

1. **`sd_max` was not redefined to mean minutes.** Reinterpreting `sd_max = 50`
   from 50 points to 50 minutes would silently disable the filter for every
   existing script. Added `sd_max_minutes` / `smoothing_width_minutes` as
   separate overriding parameters instead — non-breaking and explicit about unit.
2. **Provenance sidecar is TSV, not JSON.** `jsonlite` is not in `Imports`;
   `tidytable::fwrite` already is. Same information, no new dependency.
3. **Deprecation uses `warning()`, not `lifecycle`.** `lifecycle` is in
   `Suggests`, so every call would need a `requireNamespace()` guard. `R/cascade_defaults.R`
   has a four-line `deprecate_arg()` instead.
4. **`run` is not the first positional argument.** Making it first would break
   every existing positional `f(file, features)` call. It is a named argument,
   and `file` additionally accepts a `cascade_run` (a run object is never a path,
   so this is unambiguous). Both `f(run = run)` and `f(run)` work.

### Bug 7, found during the work

`R/get_peaks.R` called `remove_bad_peaks(pks)` while the function requires
`remove_bad_peaks(pks, n)`, where `n` bounds the valid apex-index range. This
errors at peak detection. The `n = nrow(chrom_list[[sample]])` argument was
restored. This was introduced between the plan being written and the refactor
starting — worth a look at how, since it would have broken every run.

### Dependency changes that made the pipeline runnable

Two problems blocked the package from even loading, both resolved:

1. **`tima` moved from `Imports` to `Suggests`.** It was a hard dependency but is
   used in exactly two places, both `get_last_version_from_zenodo()` calls that
   download LOTUS for the *annotation* path (`generate_ids()`,
   `generate_tables()`). The peak-integration and alignment workflow never
   touches it. Both call sites now go through `R/get_lotus.R`, which checks for
   the package and explains how to install it if genuinely needed. The package
   now loads without `tima`.

2. **`minpack.lm` added to `Imports`, and `nlsLM()` qualified.** It was called
   unqualified in `get_peaks()` and was not declared anywhere. Because every fit
   is wrapped in `try()`, a missing `nlsLM` did **not** error — it silently fell
   back to the starting estimates, so peaks were still reported, just never
   actually fitted. This is now `minpack.lm::nlsLM()` with a declared dependency.

### Verification status — all green

Run with R 4.5.3 on Windows, `pkgload::load_all()` + `testthat::test_dir()`:

```
failed : 0
error  : 0
passed : 122
skipped: 0
```

**Static checks:**
- All 83 R files parse; `Collate` lists every file exactly once.
- Every call site in `R/` and `scripts/cascade_wrapper_LLE.r` uses argument names
  that exist in the callee (catches renames that missed a caller).
- A static audit confirms no shared parameter carries a literal default, except
  genuine name collisions (`query_wikidata`'s HTTP `headers`,
  `prepare_hierarchy`'s `detector = "ms"` label, neutral `shift = 0` in plotting
  helpers).

**Executed against the bundled example data:**
- `cascade_run(show_example = TRUE)` detects 49 peaks over 0.883–31.317 min and
  matches 3800 (feature × peak) units.
- `check_peaks_integration()` returns a plotly widget carrying a 3800-row `peaks`
  attribute and a `cascade_params` object.
- `process_compare_peaks(show_example = TRUE)` writes all three files, returns
  13 informed / 2 not-informed features with identical schemas, and scores
  spanning −0.895 to 0.975.
- Validation fires as designed:
  `chromatogram = "baslined"` →
  ``Error : `chromatogram` must be one of "baselined", "improved", "original", not "baslined".``
- A `shapes = FALSE` run is refused by the scorer with a message telling you how
  to rebuild it.

**The key claim, tested directly:** at a non-default `sigma = 0.08` /
`fourier_components = 0.02`, the QC run and the export run produce **identical
peak apices**. This is the test that fails on pre-refactor code, because
`process_compare_peaks()` used to ignore both parameters.

Golden snapshots are committed at `tests/testthat/_snaps/golden-pipeline.md`.
They were written by this first run, so they encode *post*-refactor behaviour —
if you want proof that nothing moved, diff against a pre-refactor checkout
rather than trusting the snapshot alone.

`NAMESPACE` was hand-edited for the new exports; run `devtools::document()` to
regenerate it and the `man/*.Rd` files properly.

### Observation worth following up

Comparison scores are correlations and can be **negative** (the example run
reaches −0.895). `prepare_comparison()` filters with `comparison_score >=
min_similarity_prefilter` and carries a `## TODO check negative values` note.
A strongly negative score means the EIC is *anti*-correlated with the detector
peak, which is arguably more informative than a score near zero. Decide whether
those should be filtered on absolute value.

---

## 6. Open questions for the maintainer

1. **`get_peaks` is vendored** from chromatographR (parallelism stripped for
   Windows). Should Phase 2 renames apply to it, or should it stay verbatim to ease
   future upstream syncs? *Suggestion: rename, and document the divergence in the
   file header.*
2. **`pda_shift` defaults to `0.1` while `cad_shift` defaults to `0.05`** in
   `check_chromatograms_alignment`. Is that a real instrument difference worth
   keeping, or arbitrary?
3. **Should `cascade_run` be exported?** Phase 3.3 works either way; exporting makes
   the QC → export handoff explicit but adds public API surface.
4. **Is `docs/` intentionally gitignored** (it is Quarto's `output-dir`)? These notes
   live in `notes/` for that reason.
