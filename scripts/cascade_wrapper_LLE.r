
#===================================== INSTALLATION =========================================================================================================================================================
#
# Where should `cascade` come from? Set INSTALL_FROM and run this block once.
#
#   "loadall"  Load the source tree directly, WITHOUT installing. Picks up every
#              edit including uncommitted ones, instantly. This is the right
#              choice while you are actively changing the package.
#
#   "local"    Build and install from your local working copy. Also picks up
#              uncommitted edits, but takes a minute and needs a session restart.
#              Use when you want a properly installed package to test against.
#
#   "github"   Install from your fork on GitHub.
#              *** ONLY SEES CHANGES YOU HAVE COMMITTED AND PUSHED. ***
#              If you have local edits that are not pushed, this will silently
#              give you the OLD code. Use for other machines / collaborators.
#
#   "release"  Upstream r-universe build (Adafede's). Does NOT contain your fork's
#              changes at all. This was the original line in this script.
#
INSTALL_FROM <- "loadall"

# Path to your local clone (the folder containing DESCRIPTION), used by
# "loadall" and "local".
cascade_source_dir <- "C:/Users/legregam/OneDrive - Université de Genève/Documents/git_projects/cascade"

# Your fork, used by "github". Add a branch with "llegregam/cascade@my-branch".
cascade_github_repo <- "llegregam/cascade"

# cascade depends on Bioconductor packages (MSnbase, mzR, baseline, BiocParallel).
# These repos let a plain install resolve them without needing BiocManager.
#
# NOTE: `tima` is no longer required. It used to be a hard dependency, which
# meant cascade would not even load without it, but it is only used by the
# annotation path (generate_ids / generate_tables, to download LOTUS). Since you
# are only doing peak integration and alignment, it is now in Suggests and you
# can ignore it entirely.
options(repos = c(
  "https://adafede.r-universe.dev",
  "https://bioc.r-universe.dev",
  "https://cloud.r-project.org"
))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

switch(
  INSTALL_FROM,
  "loadall" = {
    if (!requireNamespace("pkgload", quietly = TRUE)) install.packages("pkgload")
    message("Loading cascade from source: ", cascade_source_dir)
    pkgload::load_all(cascade_source_dir)
  },
  "local" = {
    message("Installing cascade from local source: ", cascade_source_dir)
    remotes::install_local(
      cascade_source_dir,
      dependencies = TRUE,
      upgrade = "never",
      force = TRUE
    )
    # An already-attached package is NOT refreshed by reinstalling.
    message("Installed. RESTART R before continuing (Session > Restart R).")
  },
  "github" = {
    message("Installing cascade from GitHub: ", cascade_github_repo)
    message("NOTE: this only sees COMMITTED AND PUSHED changes.")
    remotes::install_github(
      cascade_github_repo,
      dependencies = TRUE,
      upgrade = "never",
      force = TRUE
    )
    message("Installed. RESTART R before continuing (Session > Restart R).")
  },
  "release" = {
    install.packages("cascade")
    message("Installed. RESTART R before continuing (Session > Restart R).")
  },
  stop("INSTALL_FROM must be one of: loadall, local, github, release")
)

# ===========================================================================================================================================================================================================

#===================================== CASCADE STARTS HERE===================================================================================================================================================
# `loadall` already put cascade on the search path, so library() would be a
# no-op there; it is needed after any of the installing modes.
if (INSTALL_FROM != "loadall") {
  library(cascade)
}
library(dplyr)

# Confirm you are running the refactored version: these are new, so if either
# line errors you are still on an old build (most likely you chose "github"
# without pushing, or skipped the session restart).
stopifnot(
  "cascade is out of date: cascade_run() not found. See INSTALL_FROM above." =
    exists("cascade_run"),
  "cascade is out of date: min_peak_height missing from defaults." =
    "min_peak_height" %in% names(cascade_defaults)
)

# Check documentation for the package
#  ?cascade
#  ?cascade_run
#  ?cascade_params
#  help(package = "cascade")

# =================================== PATHS ==================================================================================================================================================================

data_path <- "J:/COMMON FASIE-FATHO/Loic/Cascade/"              #"path/to/your/data"
filename_pos <- "JA_LL_arnica_1000mg_10g_pos.mzML"          #your indexed mzml in pos mode
filename_neg <- "JA_LL_arnica_1000mg_10g_neg.mzML"          #your indexed mzml in neg mode
feature_table_pos <- "arnica_pos_full_feature_table.csv"                #your MZmine features table in pos mode
feature_table_neg <- "arnica_neg_full_feature_table.csv"                #your MZmine features table in neg mode
annotations_results_pos <- "arnica_pos_results.tsv"      #your TimaR annotations results table in pos mode
annotations_results_neg <- "arnica_neg_results.tsv"      #your TimaR annotations results table in neg mode

# ===========================================================================================================================================================================================================
# Do not modify this section
data_in_path       <- file.path(data_path, "in/")
file_negative      <- file.path(data_in_path, filename_neg)
file_positive      <- file.path(data_in_path, filename_pos)
features_pos       <- file.path(data_in_path, feature_table_pos)
features_neg       <- file.path(data_in_path, feature_table_neg)
annotations_pos    <- file.path(data_in_path, annotations_results_pos)
annotations_neg    <- file.path(data_in_path, annotations_results_neg)

export_path <- file.path(data_path, "results")
if (!dir.exists(export_path)) {
  dir.create(export_path)
}
# ================================= Validate the file paths =================================================================================================================================================
stopifnot(
  "Negative mzML file does not exist in the specified data_path." = file.exists(file_negative),
  "Positive mzML file does not exist in the specified data_path." = file.exists(file_positive),
  "Features table in pos does not exist in the specified data_path." = file.exists(features_pos),
  "Features table in neg does not exist in the specified data_path." = file.exists(features_neg),
  "Annotation results in pos do not exist in the specified data_path." = file.exists(annotations_pos),
  "Annotation results in neg do not exist in the specified data_path." = file.exists(annotations_neg)
)

# ===========================================================================================================================================================================================================
# ================================= HEADERS (CRITICAL) ======================================================================================================================================================

# Adjust these strings if your mzML uses different channel IDs.
headers <- c(
  bpi = "BasePeak_0",
  pda = "PDA#1_TotalAbsorbance_0",
  cad = "UV#1_CAD_1_0"
)

# ===========================================================================================================================================================================================================
# ================================= WHAT CHANGED IN THIS SCRIPT =============================================================================================================================================
#
# The pipeline was refactored so that QC and export share one code path. Three
# things follow from that, and this script is written to exploit all three.
#
# 1. cascade_run() does the work ONCE.
#    Loading, conditioning the trace, detecting peaks and matching features to
#    them used to be repeated inside every entry point. It is now one function.
#    You build a run, LOOK at it, then EXPORT it -- and the exported peaks are
#    literally the objects you inspected, not peaks recomputed from the same
#    settings and hoping they agree.
#
#    Previously process_compare_peaks() accepted only 15 parameters and silently
#    ran sigma / smoothing_width / fit / sd_max / baseline_method at package
#    defaults. So if you tuned the QC view, the exported TSVs did NOT match it.
#    That gap is closed; the run object closes it by construction.
#
# 2. Parameters are validated at the door.
#    A typo like chromatogram = "baslined" now errors immediately, naming the
#    parameter and the bad value, instead of failing several calls deeper with
#    an unrelated message.
#
# 3. Defaults are defined exactly once, in `cascade_defaults`.
#    Print it to see every knob and its current value:  print(cascade_defaults)
#
# NEW KNOBS worth knowing about:
#    min_peak_height          detection sensitivity, as a fraction of the trace
#                             maximum. THIS is the sensitivity control; you no
#                             longer have to abuse fourier_components for it.
#    sd_max_minutes           sd_max in minutes instead of grid points, so it
#    smoothing_width_minutes  keeps its physical meaning if you change
#                             frequency or resample.
#    min_area_absolute        absolute area cutoff alongside the relative one.
#
# DEPRECATED: `type =` is now `chromatogram =`; `noise_threshold` is gone (it
# never did anything -- use min_peak_height).
#
# For the full parameter reference see notes/pipeline_concepts.md and
# notes/peak_pipeline.md.
# ===========================================================================================================================================================================================================

# ================================= VARIABLES ===============================================================================================================================================================

# 0.1. Detector-vs-MS time offsets, minutes. Read these off the alignment plot
#      in section 1. This is the single most important parameter for matching:
#      if it is off by more than a peak width, nothing matches and every feature
#      lands in the "NotInformed" file.
cad_shift <- -0.045
pda_shift <- -0.04
ms_shift  <- 0.0

# 0.2. Retention-time window to consider, minutes.
time_min <- 0.7
time_max <- 27

# 0.3. Acquisition grid. `frequency` MUST match your detector's real rate in Hz;
#      it sets the resampling step 1 / (frequency * 60 * resample) minutes.
frequency <- 1
resample  <- 1

# ===========================================================================================================================================================================================================
# ================================= ALIGNMENT ===============================================================================================================================================================

# 1. Verify the detector traces line up, and read off cad_shift / pda_shift.
result_plot <- check_chromatograms_alignment(
  file_negative = file_negative,
  file_positive = file_positive,
  headers = headers,
  show_example = FALSE,
  time_min = time_min,
  time_max = time_max,
  cad_shift = cad_shift,             # what you are tuning here
  pda_shift = pda_shift,             # what you are tuning here
  # --- signal conditioning ---------------------------------------------------
  improve_signal = TRUE,
  fourier_components = 0.1,          # keep lowest 10 % of frequencies (lower = smoother)
  frequency = frequency,
  resample = 10,                     # denser grid, purely to make the overlay smooth
  intensity_floor = 0.001,
  smoothing_width = 20,              # grid points
  k2 = 250,
  k4 = 1250000,
  sigma = 0.001,
  baseline_method = "peakDetection",
  # --- display ---------------------------------------------------------------
  chromatograms = c("bpi_pos", "cad_pos", "bpi_neg", "pda_neg"),
  type = "baselined",
  normalize_intensity = TRUE,
  normalize_time = FALSE
)

print("Alignment check completed. Resulting plot:")
print(result_plot)

# ===========================================================================================================================================================================================================
# ================================= TUNING PROFILE ==========================================================================================================================================================

# One place to define how the trace is processed, so every detector below is
# treated identically unless you deliberately override it. Previously these
# values were copy-pasted into six separate calls and had already diverged
# (min_area was 0.002 in the QC calls but 0.001 in the export calls, so the
# peaks being validated were not the peaks being exported).
#
# `cascade_params()` validates immediately -- if something here is wrong you
# find out now, not thirty minutes into an EIC extraction.

make_params <- function(detector, shift, min_area, ...) {
  cascade_params(
    detector = detector,
    chromatogram = "baselined",
    headers = headers,
    # --- (A) grid ------------------------------------------------------------
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    shift = shift,
    # --- (B) signal conditioning ---------------------------------------------
    improve_signal = TRUE,
    fourier_components = 0.01,       # tune this FIRST: lower = smoother
    intensity_floor = 0.001,
    smoothing_width_minutes = 8 / 60, # ~8 s, stated in minutes so it survives
    #                                   a change of frequency/resample
    k2 = 250,
    k4 = 1250000,
    sigma = 0.05,                    # raise to sharpen; back off on ringing
    baseline_method = "peakDetection",
    # --- (C) peak detection and fitting --------------------------------------
    fit = "egh",                     # EGH models tailing; correct for real LC peaks
    sd_max_minutes = 50 / 60,        # ~50 s max fitted width, in minutes
    max_iter = 1000,
    min_peak_height = 0,             # raise (e.g. 0.01) to ignore small wiggles
    # --- (D)/(E) filtering and matching --------------------------------------
    min_area = min_area,             # RELATIVE: fraction of total integrated signal
    min_area_absolute = 0,
    min_intensity = 1e4,             # ABSOLUTE cutoff on the MZmine feature
    intensity_threshold = 0.1,       # trims peak tails before shape comparison
    ...
  )
}

params_cad_pos <- make_params("cad", cad_shift, min_area = 0.002)
params_cad_neg <- make_params("cad", cad_shift, min_area = 0.002)
params_pda_neg <- make_params("pda", pda_shift, min_area = 0.002)
params_bpi_pos <- make_params("bpi", ms_shift, min_area = 0.001)
params_bpi_neg <- make_params("bpi", ms_shift, min_area = 0.001)

# ===========================================================================================================================================================================================================
# ================================= BUILD THE RUNS ==========================================================================================================================================================

# 2. Each run loads its inputs and detects peaks ONCE.
#
#    shapes / load_ms control how much work is done:
#      shapes  = TRUE  computes the per-peak normalised shapes used for MS
#                      comparison. Required by process_compare_peaks().
#      load_ms = TRUE  opens the raw MS1 scans. Also required for scoring, and
#                      the slowest input to open.
#    For a trace you only want to LOOK at (the two BPI ones below), both stay
#    FALSE and the run is much faster.

# 2.1. CAD on the POS file -- inspected AND exported, so it needs the full run.
run_cad_pos <- cascade_run(
  file = file_positive,
  features = features_pos,
  params = params_cad_pos,
  shapes = TRUE,
  load_ms = TRUE
)
print(run_cad_pos)

# 2.2. CAD on the NEG file -- also inspected and exported.
run_cad_neg <- cascade_run(
  file = file_negative,
  features = features_neg,
  params = params_cad_neg,
  shapes = TRUE,
  load_ms = TRUE
)
print(run_cad_neg)

# 2.3. PDA on the NEG file -- QC only for now, so skip the expensive parts.
run_pda_neg <- cascade_run(
  file = file_negative,
  features = features_neg,
  params = params_pda_neg,
  shapes = FALSE,
  load_ms = FALSE
)

# 2.4./2.5. BPI traces -- QC only.
run_bpi_pos <- cascade_run(
  file = file_positive,
  features = features_pos,
  params = params_bpi_pos,
  shapes = FALSE,
  load_ms = FALSE
)
run_bpi_neg <- cascade_run(
  file = file_negative,
  features = features_neg,
  params = params_bpi_neg,
  shapes = FALSE,
  load_ms = FALSE
)

# ===========================================================================================================================================================================================================
# ================================= INTEGRATION QC ==========================================================================================================================================================

# 3. Look at what was detected. Passing the run means these figures show exactly
#    the peaks the export will use -- no risk of the two drifting apart.

result_peaks_cad_pos <- check_peaks_integration(run_cad_pos)
result_peaks_cad_pos

result_peaks_cad_neg <- check_peaks_integration(run_cad_neg)
result_peaks_cad_neg

result_peaks_pda_neg <- check_peaks_integration(run_pda_neg)
result_peaks_pda_neg

result_peaks_bpi_pos <- check_peaks_integration(run_bpi_pos)
result_peaks_bpi_pos

result_peaks_bpi_neg <- check_peaks_integration(run_bpi_neg)
result_peaks_bpi_neg

# The peak table is now attached to the figure, so you can check the numbers
# behind the picture without re-running anything:
peaks_cad_pos <- attr(result_peaks_cad_pos, "peaks")
head(peaks_cad_pos)
nrow(peaks_cad_pos)   # how many (feature x peak) matches were made

# If peaks you can clearly see are MISSING, in order of likelihood:
#   1. min_area too high  -- it is RELATIVE, so one dominant peak suppresses the rest
#   2. sd_max too low     -- broad late peaks get dropped by the width filter
#   3. fourier_components too low -- neighbouring peaks merged into one
# If you get SPURIOUS peaks: lower fourier_components, or raise min_peak_height.

# ===========================================================================================================================================================================================================
# ================================= MS & CAD LINKS CONSTRUCTION =============================================================================================================================================

# 4. Score each matched feature by comparing its MS1 EIC shape to the detector
#    peak shape, and export.
#
#    Passing the SAME run object used for QC above is the whole point: what you
#    validated is what gets written. This step is the slow one (roughly 1 minute
#    per 1000 features) because it extracts one EIC per matched feature.
#
#    It now RETURNS its results (it used to return NULL), so you can inspect
#    them without re-reading the TSVs.

# 4.1. CAD and POS
compare_cad_pos <- process_compare_peaks(
  run_cad_pos,
  export_dir = export_path
)

# 4.2. CAD and NEG
compare_cad_neg <- process_compare_peaks(
  run_cad_neg,
  export_dir = export_path
)

# ===========================================================================================================================================================================================================
# ================================= INSPECT THE RESULTS =====================================================================================================================================================

# Files written: featuresInformed / featuresNotInformed, plus a params sidecar
# recording every resolved setting that produced them (reproducibility).
compare_cad_pos$files
list.files(export_path, full.names = TRUE)

# Features a peak informs, with a shape-similarity score in [0, 1].
head(compare_cad_pos$informed)

# How many features got a peak, and how many did not.
nrow(compare_cad_pos$informed)
nrow(compare_cad_pos$not_informed)

# Sanity check. If almost everything is "not informed", suspect `shift` before
# you suspect the biology -- a wrong offset means no interval overlaps at all.
with(
  compare_cad_pos,
  cat(sprintf(
    "informed: %d / %d (%.1f%%)\n",
    nrow(informed),
    nrow(informed) + nrow(not_informed),
    100 * nrow(informed) / max(1, nrow(informed) + nrow(not_informed))
  ))
)

# The strongest matches: high score AND a large peak area. Remember the score is
# scale-free (it compares SHAPE, not abundance), so it must always be read
# alongside peak_area, never instead of it.
compare_cad_pos$informed |>
  filter(!is.na(comparison_score)) |>
  arrange(desc(comparison_score)) |>
  select(feature_id, feature_mz, feature_rt, peak_area, comparison_score) |>
  head(20)

# The exact settings that produced this, for your methods section:
print(compare_cad_pos$params)
