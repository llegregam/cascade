# Cascade peak pipeline — annotated walkthrough

Reference documentation for the two exported entry points
[`check_peaks_integration()`](../R/check_peaks_integration.R) and
[`process_compare_peaks()`](../R/process_compare_peaks.R), plus every helper they
call. Written for a computational-MS audience who thinks in Python; R idioms are
translated to `pandas`/`numpy`/`scipy` as we go.

---

## 0. Orientation

Both functions answer the same practical question — *"do the peaks I detect in an
auxiliary detector trace (CAD / PDA / BPI) actually correspond to the LC–MS
features I care about?"* — but stop at different points:

| | `check_peaks_integration` | `process_compare_peaks` |
|---|---|---|
| Purpose | **Visual QC** of peak detection | **Quantitative scoring** + export |
| Output | one interactive `plotly` figure | two `.tsv` files (features *with* / *without* a matched peak) |
| Reads MS1 data? | No — only the detector chromatogram | Yes — pulls EICs from the raw MS1 to score each match |
| Cost | seconds | minutes (EIC extraction dominates) |

They share the **entire front half** of the pipeline (load → prepare features →
preprocess chromatogram → detect peaks → join peaks to features). `process_compare_peaks`
then adds an MS1-extraction + shape-comparison + export tail.

### Dataflow

```
                    ┌─────────────── shared front half ───────────────┐
 mzML ──load_chromatograms──► detector trace (rtime, intensity)
                                        │
                                 preprocess_chromatograms
                                 (FFT low-pass → derivative sharpen → baseline)
                                        │
 features.csv ──load_features──► prepare_features (filter DETECTED, min_intensity)
                                        │                    │
                                        ▼                    ▼
                                 preprocess_peaks ──► get_peaks (EGH/Gaussian NLS fit)
                                        │                    │
                                        │            join_peaks (integrate, min_area)
                                        │                    │
                                        │            foverlaps → features matched to peaks
                    └────────────────────┬───────────────────────────┘
                                         │
             ┌───────────────────────────┴────────────────────────────┐
             ▼                                                          ▼
   check_peaks_integration                                   process_compare_peaks
   normalize + plot_peak_detection                 extract_ms_progress (EIC per feature)
       → plotly figure                             extract_ms_peak → MSnbase Chromatogram
                                                   compare_peaks (compareChromatograms)
                                                   → *_featuresInformed_*.tsv
                                                   → *_featuresNotInformed_*.tsv
```

---

## 1. R-for-Python primer (idioms used everywhere)

| R | Python analogue | Notes |
|---|---|---|
| `x \|> f(y)` | `f(x, y)` | Native pipe: LHS becomes the **first** positional arg. `x \|> f(y=2)` is `f(x, y=2)`. |
| `f(x)` where `x` on its own line then `\|>` | method chaining | Cascade uses `tidytable`, a `data.table`-backed clone of `dplyr`/`pandas`. |
| `tidytable::mutate(df, a = b/2)` | `df.assign(a=df.b/2)` | Adds/overwrites a column. |
| `tidytable::filter(df, x > 0)` | `df[df.x > 0]` | Row filter. |
| `tidytable::select(df, new = old)` | `df.rename(columns={'old':'new'})[[...]]` | Select **and** rename in one call. |
| `switch(key, "a"=1, "b"=2)` | `{"a":1,"b":2}[key]` | Value dispatch. Also used as control flow. |
| `purrr::map(xs, f)` | `[f(x) for x in xs]` / `list(map(f, xs))` | Always returns a list. `.progress=TRUE` adds a progress bar. |
| `~rtime` (a "formula") | `lambda df: df.rtime` | In `plotly`/`plot_ly`, `~col` means "this column of the bound data". |
| `obj@slot` | `obj.attr` | `@` accesses an **S4** object slot (MSnbase objects are S4). |
| `.id = "id"` in `bind_rows` | `pd.concat(..., keys=...)` then reset_index | Adds a column recording which list element each row came from. |
| `[[1]]` | `x[0]` | Extract **one** element (unwraps the list). `[1]` keeps it wrapped. |
| `1:n`, `seq_along(x)`, `seq_len(n)` | `range`-like | **1-indexed**, inclusive. `seq_along(x)` = `range(1, len(x)+1)`. |
| `diff(y)` | `np.diff(y)` | Consecutive differences (length `n-1`). |
| `stats::approxfun(x,y)` | `scipy.interpolate.interp1d(x,y)` | Returns a callable interpolator. |

**Index base matters here.** Peak `start`/`end`/`rt` are integer positions into the
chromatogram vector and are **1-based**; keep that in mind whenever you see them
compared to `length(...)`.

---

## 2. `check_peaks_integration()` — line by line

Signature: [check_peaks_integration.R:54-85](../R/check_peaks_integration.R#L54-L85). ~30
parameters, all with defaults, grouped as: I/O (`file`, `features`), detector selection
(`detector`, `chromatogram`, `headers`), feature filtering (`min_intensity`), DSP knobs
(`fourier_components`, `k2/k4/sigma/smoothing_width`, `baseline_method`, `improve_signal`,
`intensity_floor`), time window (`time_min/max`, `frequency`, `resample`, `shift`), and
peak-fit knobs (`fit`, `sd_max`, `max_iter`, `noise_threshold`, `min_area`,
`intensity_threshold`).

> ⚠️ The roxygen `@title` reads *"Check chromatograms alignment"* — a copy-paste leftover
> from `check_chromatograms_alignment()`. The function has nothing to do with alignment;
> it's peak-integration QC.

```r
chromatograms_all <- file |> load_chromatograms(show_example, headers)   # L87-88
```
Open the mzML, pull the three detector chromatograms named in `headers`. → §4.1.

```r
name <- file |> load_name(show_example)                                  # L91-92
```
`basename(file)` — the sample name string, later used to address the per-sample columns
in the feature table. → §4.2.

```r
feature_table <- features |> load_features(show_example)                 # L95-96
df_features   <- feature_table |> prepare_features(min_intensity, name)  # L99-100
```
Read the MZmine "comprehensive export" CSV and reshape it to the columns cascade needs,
keeping only `DETECTED` features above `min_intensity`. → §4.3.

```r
switch <- switch(detector, "bpi"=headers["bpi"], "cad"=headers["cad"], "pda"=headers["pda"])
list   <- chromatograms_all[switch |> names()]                           # L103-109
```
Pick the one detector requested and subset the chromatogram list to it. Note two
variables named `switch` and `list` — they **shadow** the base-R functions of the same
name. Harmless inside this scope, but jarring; treat as legacy style.

```r
chromatograms_list <- preprocess_chromatograms(...)                      # L110-127
```
The DSP core: optional FFT low-pass + derivative sharpening (`improve_signal`), then
baseline subtraction. Returns a named list with **three versions** of the trace
(`original`, `improved`, `baselined`), each in both "wide" (list of data.frames) and
"long" (`_long`, row-bound with an `id`) form. → §4.4.

```r
peaks <- preprocess_peaks(                                               # L129-153
  df_features = df_features,
  df_long = switch(chromatogram, ...) |> mutate(intensity = intensity / max(intensity)),
  df_xy   = switch(chromatogram, ...)[[1]],
  ...
)
```
The heart of it. `switch(chromatogram, ...)` selects which of the three trace versions to
work on (default `"baselined"`). `df_long` is normalized to a 0–1 peak height; `df_xy` is
the same trace as an (rtime, intensity) matrix used for fitting. `preprocess_peaks` detects
peaks, integrates them, and overlaps them with the features. → §4.5.

```r
chromatogram_normalized <- switch(chromatogram, ...) |>                  # L155-167
  bind_rows() |> mutate(intensity = intensity / max(intensity))
```
Rebuilds the chosen trace as one long, height-normalized table for plotting.
Quirk: the `"original"` branch also does `filter(row_number() %% 10 == 1)` — it keeps
**every 10th sample** to thin an un-resampled raw trace. The `improved`/`baselined`
branches don't (they're already resampled).

```r
peaks_normalized <- peaks$list_df_features_with_peaks_long |>            # L169-174
  bind_rows() |>
  mutate(intensity = intensity_max / max(intensity_max),
         peak_max  = peak_max  / max(peak_max))
```
Height-normalize the detected-peak apex markers to the same 0–1 axis.

```r
approx_f <- stats::approxfun(x = ...$rtime, y = ...$intensity)           # L176-181
```
A linear interpolator over the normalized trace — `scipy.interpolate.interp1d`. Used by
the plot to look up the trace's y-value at a peak's start/end retention time.

```r
chromatogram_normalized |> plot_peak_detection(df2 = peaks_normalized, fun = approx_f)  # L183-184
```
Draw the figure: signal line + apex/start/end markers. This return value **is** the
function's output (an interactive plotly widget). → §4.9.

---

## 3. `process_compare_peaks()` — line by line

Signature: [process_compare_peaks.R:36-56](../R/process_compare_peaks.R#L36-L56). Same
front-half knobs, minus the QC-only ones, plus `export_dir`.

**Lines 57-124 are the shared front half** — identical in spirit to §2:
`load_ms_data` (new — §4.2), `load_chromatograms`, `load_name`, `load_features`,
`prepare_features`, `preprocess_chromatograms`, `preprocess_peaks`. Two differences worth
noting:

- **L57-59** additionally loads the raw MS1 data (`ms_data`), because this function needs
  to pull EICs later. → §4.2.
- **L74-79**: in example mode it `slice_sample(n=10)` (seeded) to keep the demo fast —
  `df.sample(n=10, random_state=42)`.

Then the MS-scoring tail:

```r
list_ms_chromatograms <- seq_along(peaks_prelist$list_df_features_with_peaks_long) |>   # L132-141
  extract_ms_progress(ms_data, rts = list_rtr, mzs = list_mzr,
                      nrows = ... |> map(nrow))
```
For **each feature that got matched to a peak**, extract the MS1 EIC over that feature's
(rt-window × m/z-window). This is the "longest step" (the code warns: ≈1 min per 1000
features per worker). → §4.6.

```r
list_ms_peaks <- list_ms_chromatograms |> purrr::map(extract_ms_peak)                   # L144-147
```
Convert each extracted EIC into an `MSnbase::Chromatogram` S4 object so it can be compared
with the detector peak. → §4.7.

```r
list_comparison_score <- seq_along(list_ms_peaks) |>                                     # L150-155
  purrr::map(compare_peaks, list_ms_peaks = list_ms_peaks, peaks_prelist = peaks_prelist)
comparison_scores <- list_comparison_score |> purrr::flatten()                           # L158-159
```
Score the **shape similarity** between the detector-trace peak and the MS1 EIC for each
match, using `MSnbase::compareChromatograms(method="closest")`. → §4.8. `flatten()`
collapses the nested list to one score per match.

```r
list_df_features_with_scores <- Map(function(df, score){ df$comparison_score <- score; df },  # L164-171
                                    peaks_prelist$list_df_features_with_peaks_long,
                                    comparison_scores)
```
`Map` = element-wise `zip`: staple each score back onto its feature-with-peak table.

```r
if (length(list_df_features_with_scores) == 0) {                                         # L173-180
  df_features_with_peaks <- peaks_prelist$df_features_without_peaks[0, ]   # empty, right schema
  df_features_with_peaks$comparison_score <- numeric(0)
} else {
  df_features_with_peaks <- bind_rows(list_df_features_with_scores)
}
```
Guard for the "no feature matched any peak" case, preserving column names so the later
`select`/export doesn't blow up. `[0, ]` is an empty slice that keeps the schema — like
`df.iloc[0:0]`.

```r
df_features_with_peaks_scored <- df_features_with_peaks |>                                # L182-197
  select(sample = id, peak_id, peak_rt_min = rt_min, ..., comparison_score) |> distinct()
df_features_without_peaks_scored <- peaks_prelist$df_features_without_peaks |>            # L199-215
  mutate(comparison_score = NA) |> select(same columns) |> distinct()
```
Two tidy output tables with human-readable column names — matched features (with a real
score) and unmatched features (score `NA`).

```r
check_export_dir(export_dir)                                                             # L217-218
```
Create the export directory tree if missing. → §4.10.

```r
df_features_with_peaks_scored    |> fwrite(".../<name>_featuresInformed_<detector>.tsv")    # L222-232
df_features_without_peaks_scored |> fwrite(".../<name>_featuresNotInformed_<detector>.tsv") # L234-244
```
Write both TSVs. The filename is built from the sample name with
`gsub("\\.[^.]+$", "", name)` — a regex that strips the file extension
(`re.sub(r"\.[^.]+$", "", name)`), then `paste`d with the detector tag. The function
returns invisibly; its product is the two files on disk.

---

## 3bis. Complete parameter reference

Every parameter of both entry points, what it physically means, and where in the
pipeline it actually bites.

### 3bis.1 I/O and mode

| Parameter | Default | In | Meaning |
|---|---|---|---|
| `file` | `NULL` | both | Path to the `.mzML`. Read twice: once by `mzR` for the detector chromatograms, once (compare only) by `MSnbase` for the MS1 scans. Ignored when `show_example = TRUE`. |
| `features` | `NULL` | both | Path to the MZmine "comprehensive export" CSV. |
| `show_example` | `FALSE` | both | Use the bundled demo `.rds` files in `inst/extdata` instead of `file`/`features`. `process_compare_peaks` additionally subsamples 10 features (seeded 42). |
| `export_dir` | `"data/interim/peaks"` | compare | Where the two `.tsv` files land. Created if missing. |
| `detector` | `"cad"` | both | Which detector channel to analyse: `"cad"`, `"pda"` or `"bpi"`. Selects a key of `headers`, and is also pasted into the output filenames. |
| `headers` | `c(bpi=..., pda=..., cad=...)` | both | Named vector mapping the short detector key → the vendor's `chromatogramId` string in the mzML. **Change this if your instrument exports different channel IDs** — everything else keys off the short names. |
| `chromatogram` / `type` | `"baselined"` | check / compare | Which of the three processed trace versions to run peak detection on: `"original"` (raw), `"improved"` (FFT + sharpened), `"baselined"` (improved + baseline-subtracted). Same parameter, different name in the two functions. |

### 3bis.2 Feature filtering

| Parameter | Default | Meaning |
|---|---|---|
| `min_intensity` | `1E4` | Absolute cutoff on the feature's `intensity_range:max` in this sample. Applied in `prepare_features` (§4.3), *before* any peak matching. Features below it never enter the pipeline. |
| `min_area` | `0.005` | **Relative** cutoff on integrated peak area, applied in `join_peaks` (§4.5.2): a peak is kept if `integral / sum(all integrals in the sample) >= min_area`. So `0.005` = "at least 0.5 % of the total integrated detector signal". Not an absolute area — see gotcha #3. |
| `intensity_threshold` | `0.1` | check only. After each peak is cut out and min-max scaled to [0, 1] (`normalize_chromato`, §4.5.3), points below this fraction of the peak height are dropped before shape comparison. `0` keeps everything. Trims the tails so the comparison focuses on the peak body. |

### 3bis.3 Time axis

| Parameter | Default | Meaning |
|---|---|---|
| `time_min` | `0.5` | Start of the analysis window, **minutes**. Applied when resampling in `improve_signal` (§4.4.1). Usual purpose: cut the void volume / injection front. |
| `time_max` | `32.5` | End of the analysis window, minutes. The effective end is `min(max(rtime), time_max)`. |
| `frequency` | `1` | Detector acquisition frequency in **Hz**. Together with `resample` it defines the uniform resampling grid step: `1 / (frequency * 60 * resample)` minutes. Set it to your detector's true rate — if it's wrong, the resampled grid over- or under-samples the trace. |
| `resample` | `1` | Multiplier on that grid density. `2` = twice as many points as the native rate (interpolated, adds no information but gives the NLS fitter more points per peak). |
| `shift` | `0.05` | Time offset in minutes added to the detector `rtime`, to compensate the dead volume between the MS source and the auxiliary detector. Applied twice, deliberately: once to the trace (`preprocess_chromatograms` L86) so peaks land at MS-comparable times, and once when building the EIC extraction window (`prepare_rt`, §4.5). **Instrument-specific — measure it once on a standard and set it.** |

### 3bis.4 Signal improvement (DSP) — `check_peaks_integration` only

`process_compare_peaks` does not forward these; `preprocess_chromatograms` runs
with its own defaults there (note `frequency` defaults to `2` in that function,
vs `1` in the entry points).

| Parameter | Default | Meaning |
|---|---|---|
| `improve_signal` | `TRUE` | Master switch. `FALSE` skips FFT + sharpening entirely; the "improved" trace then just = the original cropped to `[time_min, time_max]`. Set `FALSE` if your detector trace is already clean, or to check that sharpening isn't inventing peaks. |
| `fourier_components` | `0.01` | Fraction of Fourier coefficients kept by the brick-wall low-pass in `filter_fft` (§4.4.1a). `0.01` = keep the lowest 1 % of frequencies at each end of the spectrum. **Lower = more smoothing**, and too low will merge closely-eluting peaks. |
| `intensity_floor` | `0.001` | If the trace has any intensity ≤ 0, everything is shifted up by `abs(min) + intensity_floor` so values are strictly positive (`improve_signal` L39-44). Purely numerical hygiene for the log/derivative steps. |
| `smoothing_width` | `8` | Window `k` of the centered running mean (`caTools::runmean`) applied at four points in `signal_sharpening`. Higher = smoother but lower resolution. In grid points, not minutes — its physical width depends on `frequency`/`resample`. |
| `sigma` | `0.05` | Global gain on the sharpening correction. **Higher = stronger sharpening**, i.e. narrower peaks and more risk of ringing/negative lobes. |
| `k2` | `250` | Divisor on the 2nd-derivative term: the correction is `- (sigma/k2) * f''`. **Higher `k2` = weaker** sharpening from the second derivative. |
| `k4` | `1250000` | Divisor on the 4th-derivative term: `+ (sigma/k4) * f''''`. Higher = weaker. The 2nd-derivative term does the narrowing; the 4th-derivative term corrects the side lobes the 2nd introduces. |
| `baseline_method` | `"peakDetection"` | Passed to `baseline::baseline()`. Alternatives include `"als"`, `"rollingBall"`, `"modpolyfit"`, `"irls"`, `"fillPeaks"`, `"medianWindow"`, `"rfbaseline"`, `"shirley"`, `"TAP"`, `"lowpass"`. Only affects the `"baselined"` trace. |

> The sharpening formula is `f_sharp = f − (σ/k₂)·f'' + (σ/k₄)·f''''`, the classic
> even-derivative resolution enhancement. Practical tuning order: get
> `fourier_components` right first (smooth but peaks still resolved), then raise
> `sigma` until peaks separate, then back off if you see negative dips flanking
> large peaks.

### 3bis.5 Peak fitting — `check_peaks_integration` only

| Parameter | Default | Meaning |
|---|---|---|
| `fit` | `"egh"` | NLS peak model. `"egh"` = Exponentially-modified Gaussian Hybrid, which has a `tau` parameter for chromatographic **tailing** — the right default for real LC peaks. `"gaussian"` = symmetric Gaussian. `"raw"` = no fit at all; bounds come straight from the derivative crossings and area is trapezoidal on the raw signal (fastest, most robust, least precise apex). |
| `sd_max` | `50` | Drops fitted peaks whose fitted `sd` exceeds this. Units are **grid points at fit time** (the filter runs before index→time conversion), so its physical meaning again depends on `frequency`/`resample`. Guards against a diverging fit swallowing the whole chromatogram. |
| `max_iter` | `1000` | `nls.control(maxiter=)` for `minpack.lm::nlsLM`. With `warnOnly = TRUE`, a peak that doesn't converge returns its last iterate rather than erroring — so raising this rarely fixes anything; a bad fit usually means bad bounds. |
| `noise_threshold` | `0.001` | **Effectively inert.** It is only forwarded to `get_purity()`, which is force-disabled for single-column input — and `peaks_progress` always builds a single-column matrix. See gotcha #1. |

---

## 4. The shared toolbox

Each helper documented once; both entry points reference these.

### 4.1 `load_chromatograms()` — [R](../R/load_chromatograms.R)

Opens the mzML and returns a named list (`bpi`/`pda`/`cad`) of chromatogram data.frames.

- **L21-52** example branch: `readRDS` a bundled `.rds`, then a "dirty fix" that renames
  the first two columns to `rtime`/`intensity`.
- **L53-66** real branch:
  - `mzR::openMSfile(file)` → a file handle (lazy; doesn't read everything).
  - `mzR::chromatogramHeader()` → metadata table of all stored chromatograms.
  - `file_headers$chromatogramIndex[file_headers$chromatogramId %in% headers]` — boolean
    mask + index, i.e. "the indices whose `chromatogramId` is one of my requested header
    strings". Pure `numpy`-style fancy indexing.
  - `mzR::chromatograms(chrom = indices)` pulls just those. Result is `names()`-tagged with
    `bpi/pda/cad`.

The `headers` vector maps a friendly detector key to the vendor's channel string, e.g.
`"cad" = "UV#1_CAD_1_0"`. If your instrument exports different IDs, this is the knob to change.

### 4.2 `load_name()` / `load_features()` / `load_ms_data()`

- **`load_name`** [R](../R/load_name.R): just `basename(file)` (or a hardcoded default in
  example mode). This string is the join key into the feature table's per-sample columns.
- **`load_features`** [R](../R/load_features.R): `tidytable::fread(file)` — a fast CSV
  reader (`pd.read_csv`), or a bundled `.rds` in example mode.
- **`load_ms_data`** [R](../R/load_ms_data.R) *(only `process_compare_peaks`)*:
  `MSnbase::readMSData(file, mode="onDisk", msLevel.=1)`. `onDisk` means MSnbase keeps only
  an index in memory and reads spectra on demand — essential for full LC–MS files. Only
  MS1 is loaded.

### 4.3 `prepare_features()` — [R](../R/prepare_features.R)

Reshapes the wide MZmine export into cascade's canonical feature schema.

- **L12-27** `select(all_of(c(feature_id="id", rt="rt", mz="mz", area="area", status=..., rt_1=..., rt_2=..., mz_min=..., mz_max=..., intensity_min=..., intensity_max=...)))`.
  The per-sample columns are addressed by pasting the sample `name` into the MZmine
  column convention, e.g. `paste0("datafile:", name, ":rt_range:min")`. This is why
  `load_name` matters — it selects *this sample's* columns out of a multi-sample export.
- **L28** `filter(status == "DETECTED")` — drop gap-filled/predicted features.
- **L30** cast everything to numeric.
- **L33-34** `filter(intensity_max >= min_intensity)` — the intensity cutoff.
- **L38** `setkey(df_features, rt_1, rt_2)` — sets the `data.table` sort key on the RT
  range so the later `foverlaps` interval join is O(log n). `rt_1`/`rt_2` are the
  feature's RT-range min/max (equal endpoints are fine; foverlaps treats them as a closed
  interval).

### 4.4 `preprocess_chromatograms()` — [R](../R/preprocess_chromatograms.R)

Produces three versions of the trace. Given `list` (one detector's chromatogram) and `name`:

- **L52** `chromatograms_original <- list`.
- **L54-69** if `improve_signal`: `improve_signals_progress` runs the FFT + sharpening
  chain (§4.4.1). Else **L70-79** just crops to `[time_min, time_max]`.
- **L84-89** `chromatograms_original_long`: `bind_rows(..., .id="id")` stacks the (single)
  sample into long form with an `id` column; `rtime + shift` applies the detector time
  offset; `intensity - min(intensity)` floors the baseline to 0; `rt_1 = rt_2 = rtime`
  creates degenerate intervals so the trace can be interval-joined against features later.
- **L91-95** `chromatograms_improved_long`: same, on the improved trace (no min-subtraction).
- **L97-99** `baseline_chromatogram` (§4.4.2) applied to the improved trace → baselined.
- **L101-105** `chromatograms_baselined_long`: long form, min-subtracted.
- **L107-123** returns all six objects (three versions × {wide, long}) in a named list.

Note the naming: `_long` tables are what feed the interval joins; the plain (wide) list
elements `[[1]]` are the `(rtime, intensity)` matrices fed to the NLS peak fitter.

#### 4.4.1 The signal-improvement chain

`improve_signals_progress` [R](../R/improve_signals_progress.R) is just a progress-barred
`map` over `improve_signal` [R](../R/improve_signal.R):

1. **L39-44** shift intensities strictly positive (`min + intensity_floor` if any ≤ 0).
2. **L46-53** `filter_fft` (§4.4.1a) → a Fourier low-pass smoothed trace.
3. **L55-70** build an interpolator over the smoothed trace and **resample onto a uniform
   time grid** `seq(time_min, ..., by = 1/(frequency*60*resample))`. This is the step that
   makes downstream index↔time conversion valid (uniform Δt).
4. **L72-79** `signal_sharpening` (§4.4.1b) — even-derivative peak sharpening.
5. **L81-88** drop the first 4 points (lost to the derivative/smoothing stencils) and
   return `(rtime, intensity)`.

**4.4.1a `filter_fft`** [R](../R/filter_fft.R): classic FFT brick-wall low-pass.
```r
temp <- fft(x); keep <- round(length(temp)*components)
temp[keep:(length(temp)-keep)] <- 0            # zero the high-freq middle band
Re(fft(temp, inverse=TRUE)) / length(temp)     # inverse, normalize
```
Because a real signal's FFT is conjugate-symmetric, the **low** frequencies sit at both
ends of the coefficient vector; zeroing the middle keeps a fraction `components` (default
1%) of the lowest frequencies at each end. `numpy` equivalent: `np.fft.fft` / zero a slice
/ `np.fft.ifft(...).real`.

**4.4.1b `signal_sharpening`** [R](../R/signal_sharpening.R): even-derivative sharpening,
the standard `f - a·f'' + b·f''''` resolution-enhancement:
- `smooth_1`, `smooth_2` — successive centered running means (`caTools::runmean`).
- `deriv_2 = second_der(...)` — 2nd derivative; `smooth_3` smooths it.
- `deriv_4 = second_der(smooth_3)` — 2nd derivative of the 2nd derivative = 4th derivative;
  `smooth_4` smooths it.
- **L66-68** `sharpened = smooth_1 − (sigma/k2)·smooth_3 + (sigma/k4)·smooth_4`, with the
  `[5:]`/`[3:]` offsets realigning vectors shortened by each `diff`. Larger `k2`/`k4` =
  weaker sharpening; larger `sigma` = stronger.
- `second_der(x,y) = deriv(middle_pts(x), deriv(x,y))` where `deriv = diff(y)/diff(x)` and
  `middle_pts(x) = x[-1] - diff(x)/2` (midpoints), i.e. a finite-difference 2nd derivative
  valid on a non-uniform grid. See [second_der.R](../R/second_der.R),
  [deriv.R](../R/deriv.R), [middle_pts.R](../R/middle_pts.R).

#### 4.4.2 `baseline_chromatogram()` — [R](../R/baseline_chromatogram.R)

Wraps `baseline::baseline(spectra = t(intensity), method = method)` (default
`"peakDetection"`), replaces `NA`→0 first, and writes the `@corrected` slot back as the new
intensity. `t()` transposes the vector to the 1×N row-matrix that `baseline` expects.

### 4.5 `preprocess_peaks()` — [R](../R/preprocess_peaks.R)

Detects peaks and reconciles them with features. Given `df_xy` (matrix for fitting),
`df_long` (normalized long trace), and `df_features`:

- **L44-50** `peaks_progress` (§4.5.1) → fitted peak table (rt/start/end/height/area…).
- **L53-54** `bind_rows(peaks, .id="id")` → one long peak table.
- **L57-62** `join_peaks` (§4.5.2): overlap peaks with the trace to **integrate** each peak
  (sum of intensities inside its RT bounds) and drop peaks below `min_area`.
- **L64** `setkey(df_peaks, rt_min, rt_max)` then **L67-68**
  `foverlaps(df_features, df_peaks)` — the interval join: match each feature's RT range to
  any peak whose RT range overlaps it. Think `pandas` `merge` but on **overlapping
  intervals** rather than equality.
- **L70-73** `df_features_with_peaks`: rows where `peak_id` is not `NA` (a match), drop the
  raw `rt_1/rt_2`, dedupe.
- **L76-78** `df_features_without_peaks`: the `NA`-`peak_id` rows (features with no peak).
- **L81-99** split the matched table by sample (`id`) and then by `peak_id`, and
  `flatten()` to one list-element per (sample, peak) — this is
  `list_df_features_with_peaks_long`, the unit of work for MS extraction/scoring.
- **L102-107** `normalize_chromato` (§4.5.3) — for each peak, cut the trace to its RT window
  and min-max normalize, dropping points below `intensity_threshold`.
- **L110-113** `prepare_peaks` [R](../R/prepare_peaks.R): wrap each cut trace as an
  `MSnbase::Chromatogram` with rtime **rescaled to [0,1]** (so shape comparison is
  RT-scale-invariant; see the linked xcms issue).
- **L115-120** `prepare_rt` [R](../R/prepare_rt.R): the peak's RT window as a 1×2 matrix in
  **seconds** (`(rt + shift)*60`) — MSnbase works in seconds, cascade in minutes.
- **L122-126** `prepare_mz` [R](../R/prepare_mz.R): the feature's m/z window as an
  `(mz_min, mz_max)` matrix, one row per feature. `rts`+`mzs` together define each EIC box.
- **L128-142** returns the five artifacts (`list_df_features_with_peaks_long`,
  `list_chromato_peaks`, `list_rtr`, `list_mzr`, `df_features_without_peaks`).

#### 4.5.1 `peaks_progress()` → `get_peaks()`

`peaks_progress` [R](../R/peaks_progress.R) reshapes `df_xy` to a single-column,
row-names=rtime matrix (the format `get_peaks` wants), fits, then renames the fit columns
to cascade's convention (`peak_id`, `peak_max`, `rt_apex`, `rt_min`, `rt_max`).

**`get_peaks()`** [R](../R/get_peaks.R) is a large routine adapted from the
[`chromatographR`](https://github.com/ethanbass/chromatographR) package (parallelism
removed for Windows). The essentials:

- `find_peaks(y)` — locate candidate apices by first-derivative **zero-crossings**
  (`sign(d[i]) > sign(d[i+1])`), optionally with slope/amplitude thresholds, and derive
  each peak's lower/upper bounds from the neighbouring opposite crossings.
- `fitpk_egh` / `fitpk_gaussian` / `fitpk_raw` — fit each candidate with non-linear least
  squares (`minpack.lm::nlsLM`). EGH = Exponentially-Modified Gaussian Hybrid, the default,
  which models chromatographic **tailing** (its `tau` parameter). `nls.control(warnOnly=TRUE)`
  means a non-converging fit returns its last iterate instead of erroring.
- `sd.max` drops implausibly wide fits; `remove_bad_peaks` drops rows with a broken
  position or an out-of-range apex index; `convert_indices_to_times` maps integer
  scan indices → minutes using the (now uniform) time grid.

> **Gotchas in `get_peaks`, worth knowing before you trust a peak count:**
> 1. **`noise_threshold` is effectively a no-op.** It's only forwarded to `get_purity()`,
>    and purity estimation is force-disabled whenever the input has one column
>    (`ncol(x)==1`) — which is always true here, since `peaks_progress` builds a
>    single-column matrix. The real detection sensitivity lever is `amp_thresh` inside
>    `find_peaks`, which isn't plumbed through to the public API.
> 2. **`remove_bad_peaks` was recently fixed.** An earlier version compared the *peak count*
>    to scan indices and silently kept only the earliest peak; the current code checks
>    `complete.cases(rt,start,end)` and `1 ≤ rt ≤ n`. If you pull an older cascade build and
>    see "only one peak detected", that's the bug.

#### 4.5.2 `join_peaks()` — [R](../R/join_peaks.R)

`foverlaps(peaks, chromatograms)` attaches every trace sample inside a peak's RT window to
that peak, then per `(peak_id, id)` computes `integral = sum(intensity)` (a discrete
integration), dedupes to one row per peak, and — **per sample** — keeps peaks with
`integral / sum(integral) >= min_area`.

> `min_area` is a **relative** fraction of the sample's total integrated signal, not an
> absolute area. With one dominant peak, genuine small peaks can fall below a seemingly
> tiny cutoff like `0.001`.

#### 4.5.3 `normalize_chromato()` — [R](../R/normalize_chromato.R)

For one peak: `filter(rtime in [rt_min, rt_max])`, min-max scale intensity to [0,1], then
`filter(intensity >= intensity_threshold)`. `intensity_threshold=0` keeps all points.

### 4.6 `extract_ms_progress()` — [R](../R/extract_ms_progress.R) *(compare only)*

For each matched peak, pull the MS1 EIC:
- `safe_chromatogram` wraps `MSnbase::chromatogram(ms_data, rt=..., mz=...)` in a
  `max_attempts=10` retry loop (transient mzML read failures are common); on final failure
  it warns and returns `NULL`.
- `rt = rts[[x]]` is the peak's RT box (seconds), `mz = Reduce(rbind, mzs[[x]])` stacks the
  feature m/z windows into one matrix — `functools.reduce(np.vstack, ...)`.
- The extracted chromatogram is passed through `transform_ms` (§4.6a).

> The inner map takes `nrows` but never uses it — vestigial.

**4.6a `transform_ms`** [R](../R/transform_ms.R): for each sub-chromatogram, build a
`(intensity, rtime)` frame from the S4 `@intensity`/`@rtime` slots, drop `NA`, **min-max
normalize both intensity and rtime to [0,1]**, keep `intensity >= min_intensity` (default
0.1), and sort by rtime. The `custom_min/custom_max` helpers return `Inf` on empty input to
avoid `min()`-of-empty warnings. RT is normalized so peaks are compared on **shape**, not
absolute time (see the linked xcms issue #593).

### 4.7 `extract_ms_peak()` — [R](../R/extract_ms_peak.R) *(compare only)*

Maps each normalized EIC frame back into an `MSnbase::Chromatogram(intensity=, rtime=)` S4
object — the type `compareChromatograms` requires.

### 4.8 `compare_peaks()` — [R](../R/compare_peaks.R) *(compare only)*

For peak `x`: if it has extracted MS peaks, `MSnbase::compareChromatograms(detector_peak,
ms_peak, method="closest")` returns a similarity score per candidate; degenerate cases
(length ≤ 1, or no MS peaks) score `0`. `method="closest"` aligns the two traces by nearest
retention time before correlating — appropriate since both were RT-normalized to [0,1].

### 4.9 `plot_peak_detection()` — [R](../R/plot_peak_detection.R) *(check only)*

Builds the plotly figure: a signal **line** (`df1`), an apex marker at
`(rt_apex, peak_max)`, and start/end markers on a secondary y-axis using `fun` (the
interpolator) to read the trace height at the peak bounds. Dual y-axes let the 0–1 peak
markers overlay the 0–1 signal cleanly; legend/ticks are hidden for a clean QC look.

> **Likely bug:** the "detected minimum (end)" trace ([L42-51](../R/plot_peak_detection.R#L42-L51))
> plots `x = ~rt_max` but `y = ~ fun(rt_min)` — the **start** RT is used for the end
> marker's height. Almost certainly a copy-paste slip from the start-marker block above;
> the end marker's y should be `fun(rt_max)`. Cosmetic (marker height only), but it makes
> the end triangle sit at the wrong height on tailing peaks.

### 4.10 `check_export_dir()` — [R](../R/check_export_dir.R) *(compare only)*

Creates up to three ancestor levels of `dir` if missing (nested `dirname` calls), using
`ifelse` per level. Effectively `os.makedirs(dir, exist_ok=True)` — just spelled out by
hand.

---

## 5. Quick gotcha index

| # | Where | Issue |
|---|---|---|
| 1 | `get_peaks` | `noise_threshold` is a no-op (purity disabled for single-column input); use would-be `amp_thresh`, not exposed. §4.5.1 |
| 2 | `get_peaks` (old builds) | `remove_bad_peaks` once kept only the earliest peak — fixed in current tree. §4.5.1 |
| 3 | `join_peaks` | `min_area` is a **relative** fraction of total signal, not absolute area. §4.5.2 |
| 4 | `plot_peak_detection` | end-marker height uses `fun(rt_min)` instead of `fun(rt_max)`. §4.9 |
| 5 | `check_peaks_integration` | roxygen title is wrong ("Check chromatograms alignment"). §2 |
| 6 | both entry points | `switch`/`list` local variables shadow base-R functions. §2 |
| 7 | `check_peaks_integration` | `"original"` plot branch thins to every 10th point; `improved`/`baselined` don't. §2 |
| 8 | `extract_ms_progress` | inner map receives `nrows` but never uses it. §4.6 |
