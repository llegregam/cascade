# The ideas behind CASCADE

*Written for someone who knows metabolomics and chromatography, but not this
codebase. This document is about **why** the pipeline does what it does. For the
line-by-line code walkthrough see `notes/peak_pipeline.md`.*

---

## 1. The problem CASCADE exists to solve

You run an extract on LC–MS/MS, process it in MZmine, annotate it with TIMA, and
you get a feature table: a few thousand `(m/z, RT, area)` triplets, many with a
putative structure attached.

Now someone asks the question that actually matters for a natural-products paper:

> **Which of these compounds are actually abundant in the extract?**

And the feature table cannot answer it. MS peak area is not proportional to
amount. It is proportional to *amount × ionisation efficiency*, and ionisation
efficiency varies by **three to four orders of magnitude** across compound
classes for reasons that have nothing to do with concentration:

- a quaternary alkaloid ionises beautifully in ESI+ and will dominate the TIC at
  trace level;
- a saponin or a fatty acid may be the bulk of the dry mass and barely show up;
- ion suppression means a compound's response depends on **what co-elutes with it**,
  so the same compound gives different areas in different samples;
- adduct and in-source fragment distributions split one compound's signal across
  several features.

So the MS gives you **identity with untrustworthy quantity**.

The Charged Aerosol Detector gives you the opposite. CAD nebulises the eluent,
dries it to particles, charges them in a corona discharge, and measures the
current. Response depends on the **mass of non-volatile material** in the
droplet, not on chemistry. It is quasi-universal: any non-volatile analyte
responds, and two different compounds at the same mass concentration give
responses within roughly a factor of two of each other — versus a factor of
10 000 in ESI. So CAD gives you **quantity without identity**.

**CASCADE's whole purpose is to bolt these two together**: take the abundance
information from the CAD trace and the identity information from the MS features,
and decide *which annotated compound is responsible for each CAD peak*. That is
what turns "we detected 3000 features" into "these six compounds are 80 % of the
extract, and here is what they are."

The same machinery works for PDA (UV absorbance, universal only for chromophores)
and for the MS BPI trace itself; CAD is just the one where the quantitative claim
is strongest.

---

## 2. Why this is harder than "match the retention times"

The naive approach — for each CAD peak, take the MS feature at the same RT — fails
for four reasons, and the whole pipeline is built around them.

**2.1 The detectors do not share a time base.**
The CAD sits downstream of the MS split (or vice versa), connected by real tubing
with real volume. At 0.4 mL/min, 20 µL of extra tubing is a 3-second offset.
That offset is a fixed property of your plumbing, and it is exactly the width of a
narrow UHPLC peak — so ignoring it means matching every peak to its neighbour.
Hence the `shift` parameter, and hence the alignment QC step that exists purely to
let you measure it.

**2.2 One CAD peak contains many MS features.**
In a plant extract, a single chromatographic peak routinely contains the
protonated molecule, the sodium adduct, an ammonium adduct, two in-source
fragments, an isotopologue MZmine failed to group, and one or two genuinely
co-eluting isomers. RT overlap alone cannot tell you which of these *is* the
compound responsible for the CAD response.

**2.3 The CAD trace is a bulk measurement.**
It has no mass axis. You cannot ask it "how much of *this* compound is here" —
only "how much material elutes here in total."

**2.4 Both traces are noisy, drifting, and differently shaped.**
CAD is noisier than UV, drifts with gradient composition, and has broader peaks
than the MS EIC because of the extra dead volume and the nebulisation process.

---

## 3. The core idea: match by shape, not just by time

The insight that makes CASCADE work is this:

> If compound X is what produces a given CAD peak, then the **shape** of X's
> extracted ion chromatogram must match the **shape** of that CAD peak.

Retention-time overlap gets you a shortlist. Peak-shape correlation tells you
which member of the shortlist is real.

This works because co-eluting-but-distinct compounds are almost never *perfectly*
co-eluting. Two compounds 3 seconds apart inside a 20-second CAD peak produce EICs
whose apices sit at different positions within that peak. Adducts and fragments of
the *same* compound, on the other hand, share the apex exactly — which is why they
all score well, and why downstream steps use annotation and taxonomy to pick one.

So the pipeline is:

```
CAD trace ──► clean it up ──► find peaks ──► integrate them
                                               │
                                               │  RT-overlap join
                                               ▼
MS features ──────────────────────────► shortlist per peak
                                               │
                                               │  for each candidate: pull its MS1 EIC
                                               ▼
                                    compare EIC shape to CAD peak shape
                                               │
                                               ▼
                                    similarity score per (peak, feature)
```

Two practical details matter here.

**Both traces are rescaled before comparison.** Intensity is min–max normalised to
[0, 1], and *so is retention time within the peak window*. That makes the
comparison scale-invariant: you are asking "is this the same shape?", not "is this
the same height at the same absolute second?" A CAD peak and an EIC of the same
compound have systematically different widths (different dead volumes, different
detector time constants), so comparing raw shapes would penalise correct matches.

**The score is a chromatogram correlation**, computed by
`MSnbase::compareChromatograms(method = "closest")`, which aligns the two traces
by nearest retention time and correlates. Scores near 1 mean the EIC tracks the
CAD peak; scores near 0 mean it does not.

---

## 4. Why the signal processing is necessary (and what each step costs you)

Peak detection in this pipeline is **derivative-based**: an apex is a
zero-crossing of the first derivative, and a peak's start and end are the
neighbouring crossings in the opposite direction. That is a clean, assumption-free
definition — and it is catastrophically sensitive to noise, because every noise
wiggle is also a zero-crossing. Raw CAD data would give you thousands of "peaks".

So the trace is conditioned first. The order is deliberate.

### 4.1 Fourier low-pass — separate signal from noise by frequency

A chromatographic peak 20 seconds wide is a **low-frequency** feature. Detector
noise is **high-frequency**. In the frequency domain they barely overlap, so an
FFT, zero the high-frequency coefficients, inverse FFT gives you a dramatically
cleaner trace with the peaks essentially intact.

`fourier_components` is the fraction of coefficients kept (0.01 = the lowest 1 %).

- **Too high** → noise survives → spurious peaks.
- **Too low** → you also remove the frequencies that *distinguish* two closely
  eluting peaks, and they merge into one.

This is a brick-wall filter, which is theoretically the sharpest possible cutoff
but can produce Gibbs ringing — small oscillations flanking sharp transitions.
In practice chromatographic peaks are smooth enough that this is rarely visible,
but it is the reason not to set the cutoff aggressively low.

### 4.2 Resampling onto a uniform grid

After filtering, the trace is interpolated onto a strictly uniform time grid of
step `1 / (frequency × 60 × resample)` minutes.

This is not cosmetic. Everything downstream — the finite-difference derivatives,
the conversion of peak indices back to retention times, the trapezoidal
integration — assumes constant Δt. Detector files often have small timing jitter,
and this step removes it.

It also means **`frequency` must match your detector's real acquisition rate**. If
you declare 1 Hz for a detector running at 5 Hz, you discard 80 % of your points
and narrow peaks lose the sampling density needed to be detected at all. As a rule
of thumb you want ≥ 15–20 points across the narrowest peak you care about.

### 4.3 Derivative sharpening — buying resolution with signal-to-noise

This is the step most people have not met before. The trace is transformed as

$$f_{\text{sharp}} = f - \frac{\sigma}{k_2} f'' + \frac{\sigma}{k_4} f''''$$

**Why this narrows peaks:** for a Gaussian, the second derivative is negative at
the apex and positive in the wings. Subtracting a scaled `f''` therefore *adds* to
the apex and *subtracts* from the wings — the peak gets taller and narrower while
its position and area are approximately preserved. The fourth-derivative term
corrects the negative side lobes that the second-derivative term introduces.

This is the classical even-derivative resolution enhancement (the same idea as
Savitzky–Golay derivative sharpening), and it can resolve shoulders that are not
separable by eye.

**What it costs:** differentiation amplifies high-frequency content — each
derivative multiplies noise power by roughly ω². A fourth derivative is a
ferocious noise amplifier. That is *why* the FFT low-pass runs first and why
running means are interleaved between every derivative step. Smooth, then sharpen.
Never the reverse.

**Tuning:** `sigma` is the overall gain (higher = stronger); `k2` and `k4` are
divisors (higher = weaker). Practical order: get `fourier_components` right first,
then raise `sigma` until shoulders separate, then back off as soon as you see
negative dips flanking large peaks or single peaks splitting in two. Those are the
signature of over-sharpening, and they will corrupt your integration.

> **Judgement call:** sharpening is a *hypothesis* about your data — that an
> apparent single peak is really two. If quantification is the goal and the
> shoulders do not matter, `improve_signal = FALSE` is a defensible, more
> conservative choice. Always compare the two before committing.

### 4.4 Baseline correction — because you cannot integrate on a slope

In reversed-phase gradient elution, the mobile phase composition changes
continuously, and both CAD and UV respond to it. The CAD baseline typically rises
through the gradient because organic solvent nebulises more efficiently. Left
uncorrected, this drift is added to every peak's integral, inflating late-eluting
peaks relative to early ones.

The default `"peakDetection"` method identifies peak-free regions and interpolates
the baseline through them. Alternatives (`"als"` — asymmetric least squares,
`"rollingBall"`, `"modpolyfit"`) differ in how they distinguish "baseline" from
"broad peak", which is genuinely ambiguous for a badly resolved region. If your
corrected trace has peaks sitting below zero, or broad humps that were clearly
real got flattened, that is the knob to change.

The pipeline keeps all three versions of the trace — `original`, `improved`
(filtered + sharpened), `baselined` — precisely so you can compare them.

---

## 5. From peaks to numbers

### 5.1 Peak fitting: why EGH and not Gaussian

Once candidate peaks are located, each is fitted by non-linear least squares.
The default model is the **Exponentially-modified Gaussian Hybrid (EGH)**.

Real chromatographic peaks are not Gaussian; they tail. Causes include secondary
interactions with residual silanols, extra-column volume, slow mass transfer, and
column overload. A symmetric Gaussian fitted to a tailing peak systematically
places the apex too early and **under-integrates the tail** — precisely the region
where a minor co-eluting compound would hide.

EGH adds a `tau` parameter describing the exponential decay of the tail, so it
follows a real peak much more faithfully. `"gaussian"` is available for genuinely
symmetric peaks; `"raw"` skips fitting entirely and integrates the signal between
the derivative-derived boundaries — the most robust and least presumptuous option,
at the cost of a less precise apex.

### 5.2 Integration and the relative-area threshold

A peak's `integral` is the sum of trace intensities between its start and end —
discrete trapezoidal integration. For CAD this quantity is, to first order,
proportional to the **mass** of material eluting in that window. That is the
semi-quantitative claim at the heart of the method.

Peaks are then filtered by `min_area`, and this is the parameter most likely to
surprise you: it is a **fraction of the sample's total integrated signal**, not an
absolute area. `min_area = 0.005` means "keep peaks worth at least 0.5 % of
everything integrated in this run."

The consequence is that **one dominant compound suppresses everything else**. In an
extract where a single peak is 60 % of the CAD signal, genuine minor constituents
can fall below even a seemingly tiny relative cutoff. If peaks you can clearly see
in the plot are missing from the output, this is the first thing to lower.

### 5.3 Matching features to peaks: an interval join

Each MZmine feature carries an RT range (`rt_range:min`, `rt_range:max`), and each
detected peak has `[start, end]`. The match is an **interval overlap** — implemented
as a `data.table::foverlaps` join, which is why the tables are sorted by their
range columns beforehand.

Every feature therefore lands in exactly one of two buckets:

- **informed** — its RT range overlaps at least one detector peak; it gets a
  similarity score;
- **not informed** — no overlap; score is `NA`.

A feature can be "not informed" for two very different reasons: either it is a
genuine trace compound below the CAD detection limit (interesting and expected),
or your `shift` is wrong (a bug). If *almost everything* is uninformed, suspect the
shift before you suspect the biology.

---

## 6. From scores to "major versus minor metabolites"

The similarity scores are the raw material; the final call combines three
independent lines of evidence.

**Chemical evidence — the shape score.** Two thresholds are applied
(`min_similarity_prefilter`, default 0.6, then `min_similarity_filter`, default
0.8). Features passing both are strong candidates for *being* the compound behind
the peak.

**Annotation confidence — the TIMA score.** Each feature's structural annotation
carries a confidence (`min_confidence`, default 0.4). Below it, the candidate is
labelled `notConfident` rather than silently trusted. Annotations are also reduced
to NPClassifier levels — pathway, superclass, class — so that even where the exact
structure is uncertain, the *chemical class* is usually right. For a
natural-products story, "this major peak is a sesquiterpene lactone" is often the
defensible claim.

**Taxonomic plausibility.** Where several candidates survive for one peak, the
pipeline prefers the one whose structure has actually been reported in the
organism under study. This is the same biological-prior logic TIMA uses. It is a
prior, not proof — it will preferentially reject genuinely novel compounds, which
is exactly the wrong behaviour if novelty is what you are hunting. Know which mode
you are in.

The output is a peak-level view of the extract: for each major CAD peak, the most
plausible identity and its share of the total detected mass. That is the
"pseudochromatogram" — a chromatogram redrawn with chemical classes rather than
raw signal.

---

## 7. Caveats worth stating in a methods section

**CAD response is not linear.** It follows roughly a power law, `response ∝ mass^b`
with `b` typically between 0.5 and 1. Comparing raw areas between compounds
overestimates the small ones. The package carries a `predict_response()` helper
implementing an inverse power-law correction, though it is not wired into the main
path.

**CAD response depends on mobile-phase composition.** Nebulisation and transport
efficiency improve with organic content, so the *same* mass gives a larger response
late in a gradient than early. Rigorous work either uses an inverse gradient or
applies a composition-dependent correction — hence the `acn`-aware terms in
`predict_response()`. Without one of those, comparing an early-eluting glycoside to
a late-eluting aglycone is approximate.

**Only non-volatile analytes respond.** Anything lost during nebulisation and
drying is invisible to CAD. "Not detected by CAD" ≠ "not present."

**One CAD peak may be several compounds.** The score identifies which annotated
feature best explains the peak; it does not prove that feature accounts for *all*
of it. A high score plus a clean peak shape is strong evidence; a high score on a
broad, tailing, obviously composite peak is not.

**The shape comparison is scale-free.** Because both axes are normalised, a
feature whose EIC has the right shape but is 1000× weaker still scores well. That
is intended — the score is about *co-elution identity*, not abundance — but it
means the score must always be read alongside the peak area, never instead of it.

---

## 8. Where to look next

| To understand... | Read |
|---|---|
| what each function does, line by line | `notes/peak_pipeline.md` |
| what each parameter does and how to tune it | the cheat-sheet in `scripts/cascade_wrapper_LLE.r` |
| known parameter inconsistencies and the plan to fix them | `notes/refactor_plan.md` |
| a worked end-to-end run | `vignettes/articles/II-processing.qmd` |
