#!/usr/bin/env Rscript
# litmus_process.R — SUSHI wrapper for the LITMUS Processing step.
# Invoked by LitmusProcessApp.rb (SUSHI Ruby frontend).
#
# Contract:
#   Input — a SUSHI dataset.tsv with one row per acquisition file.
#           Required columns (ezRun convention; SUSHI strips type suffix):
#             Name, Raw [File], Polarity [Factor], Sample Type [Factor]
#           Recommended metadata columns (read verbatim, never parsed from filenames):
#             Grouping [Factor], Sample Name [Characteristic],
#             Sample Id [B-Fabric], Resource [B-Fabric], Container [B-Fabric],
#             Injection Order [Characteristic]
#
#   Output — a single directory containing:
#             feature_table.csv   (samples × features intensity matrix)
#             feature_meta.csv    (feature_id, mz, rt, [lipid_class, annotation...])
#             sample_meta.csv     (filename, group, sample_type, polarity, batch, ...)
#             parameters.json     (provenance)
#             processing.log      (engine stderr/stdout)
#             00index.html        (rendered Rmd report)
#
# Design:
#   - Reads grouping from dataset.tsv metadata columns. NEVER parses filenames.
#   - Filters samples by polarity + sample_type (both from metadata).
#   - Dispatches to processing_engine.R::run_engine(engine, input_files, out_dir, params).
#   - Falls back to synthetic "demo mode" when input files are absent, so
#     SUSHI wiring (epilogue, report render, gstore copy) can be validated
#     end-to-end without a live mass-spec instrument.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1]) && nzchar(a[1])) a else b

# ------------------------------------------------------------------
# 1. CLI
# ------------------------------------------------------------------
option_list <- list(
  make_option("--dataset",         type = "character", help = "SUSHI dataset.tsv (required)"),
  make_option("--output",          type = "character", help = "Output directory (required)"),
  make_option("--engine",          type = "character", default = "xcms",
              help = "Engine: xcms | mzmine4 | msdial5 | openms | masster | masscube [default: %default]"),
  make_option("--polarity",        type = "character", default = "pos",
              help = "Polarity filter: pos | neg | both [default: %default]"),
  make_option("--sample_types",    type = "character", default = "sample,qc",
              help = "Comma-separated Sample Type values to include [default: %default]"),
  make_option("--mz_tol_ppm",      type = "double",    default = 5,      help = "MS1 mass tolerance ppm"),
  make_option("--rt_peakwidth_sec",type = "character", default = "5,40", help = "CentWave peakwidth range (sec)"),
  make_option("--snthresh",        type = "double",    default = 3,      help = "Signal-to-noise threshold"),
  make_option("--noise",           type = "double",    default = 1000,   help = "Intensity noise floor"),
  make_option("--min_frac",        type = "double",    default = 0.5,    help = "Min fraction of samples per feature"),
  make_option("--annotation",      type = "character", default = "none", help = "Lipid annotation: none | goslin | lipidoracle"),
  make_option("--demo_mode",       type = "character", default = "auto", help = "auto | always | never"),
  make_option("--litmus_path",     type = "character", help = "LITMUS staging directory (required)"),
  make_option("--gstore_base",     type = "character", default = "/srv/gstore/projects",
              help = "gstore root to prepend to relative [File] paths")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot("--dataset is required"     = !is.null(opt$dataset))
stopifnot("--output is required"      = !is.null(opt$output))
stopifnot("--litmus_path is required" = !is.null(opt$litmus_path))
stopifnot("dataset.tsv not found"     = file.exists(opt$dataset))
stopifnot("litmus_path not found"     = dir.exists(opt$litmus_path))

dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(opt$output, "processing.log")
log_con  <- file(log_file, open = "wt")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  writeLines(msg, log_con)
  flush(log_con)
}

log_msg("=== LITMUS Process (SUSHI) ===")
log_msg("  engine:        ", opt$engine)
log_msg("  polarity:      ", opt$polarity)
log_msg("  sample_types:  ", opt$sample_types)
log_msg("  mz_tol_ppm:    ", opt$mz_tol_ppm)
log_msg("  demo_mode:     ", opt$demo_mode)
log_msg("  litmus_path:   ", opt$litmus_path)

# ------------------------------------------------------------------
# 2. Load LITMUS processing engine (source verbatim, no modification)
# ------------------------------------------------------------------
engine_lib <- file.path(opt$litmus_path, "scripts", "R", "processing_engine.R")
if (file.exists(engine_lib)) {
  # processing_engine.R sources output_contracts.R if present; shield with try()
  # so a missing contract file (not required for XCMS path) doesn't abort.
  try(source(engine_lib), silent = FALSE)
  log_msg("[OK] sourced processing_engine.R")
} else {
  log_msg("[WARN] processing_engine.R not found at ", engine_lib,
          " — demo mode only.")
}

# ------------------------------------------------------------------
# 3. Read dataset.tsv, apply metadata filters (NEVER parse filenames)
# ------------------------------------------------------------------
ds <- read.table(opt$dataset, sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE,
                 quote = "", comment.char = "")

# SUSHI strips [File]/[Factor]/[Characteristic]/[B-Fabric] suffixes in some
# contexts but the staged dataset.tsv on the compute node still carries them.
# Build a lookup that tolerates both forms.
col <- function(df, canonical) {
  # returns the first matching column name in df for the canonical key,
  # or NA_character_ if none found.
  candidates <- c(
    canonical,
    paste0(canonical, " [File]"),
    paste0(canonical, " [Factor]"),
    paste0(canonical, " [Characteristic]"),
    paste0(canonical, " [B-Fabric]")
  )
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) hit[1] else NA_character_
}
require_col <- function(df, canonical) {
  c <- col(df, canonical)
  if (is.na(c)) stop("Required column not found in dataset.tsv: '", canonical,
                     "' (looked for variants with [File]/[Factor]/...).")
  c
}

c_name       <- require_col(ds, "Name")
c_raw        <- require_col(ds, "Raw")
c_polarity   <- require_col(ds, "Polarity")
c_smptype    <- require_col(ds, "Sample Type")
c_grouping   <- col(ds,         "Grouping")
c_injord     <- col(ds,         "Injection Order")
c_sampname   <- col(ds,         "Sample Name")
c_sampid     <- col(ds,         "Sample Id")
c_resource   <- col(ds,         "Resource")
c_container  <- col(ds,         "Container")

log_msg("[OK] dataset.tsv: ", nrow(ds), " rows, ", ncol(ds), " cols")
log_msg("     Raw col:       '", c_raw, "'")
log_msg("     Polarity col:  '", c_polarity, "'")
log_msg("     Sample Type:   '", c_smptype, "'")
log_msg("     Grouping col:  '", ifelse(is.na(c_grouping), "(not present — grouping is 'all')", c_grouping), "'")

# Apply polarity filter (metadata-driven)
if (opt$polarity != "both") {
  keep_pol <- ds[[c_polarity]] == opt$polarity
  log_msg("[filter] polarity == '", opt$polarity, "': ",
          sum(keep_pol), " / ", nrow(ds), " rows")
  ds <- ds[keep_pol, , drop = FALSE]
}

# Apply sample-type filter (metadata-driven)
allowed_types <- trimws(strsplit(opt$sample_types, ",", fixed = TRUE)[[1]])
keep_type <- ds[[c_smptype]] %in% allowed_types
log_msg("[filter] Sample Type in {", paste(allowed_types, collapse = ","), "}: ",
        sum(keep_type), " / ", nrow(ds), " rows")
ds <- ds[keep_type, , drop = FALSE]
if (nrow(ds) == 0) stop("No samples remain after filtering. Check --polarity and --sample_types.")

# Resolve Raw file paths: absolute → verbatim, relative → gstore-prefixed
raw_paths <- sapply(ds[[c_raw]], function(p) {
  p <- trimws(as.character(p))
  if (nchar(p) == 0) return(NA_character_)
  if (startsWith(p, "/")) p else file.path(opt$gstore_base, p)
}, USE.NAMES = FALSE)
ds[[".raw_abs"]] <- raw_paths

# Check input-file existence (drives demo-mode fallback)
exists_vec <- file.exists(raw_paths)
n_present  <- sum(exists_vec)
n_total    <- length(raw_paths)
log_msg("[inputs] ", n_present, " / ", n_total, " raw files exist on disk")
if (n_present > 0 && n_present < n_total) {
  log_msg("[WARN] ", n_total - n_present, " files missing — those rows will be dropped or demo-filled")
}

# Minimum number of real input files required for any feature-detection
# engine (XCMS / MZmine / MS-DIAL / OpenMS / MASSter / MassCube). Below this
# threshold, alignment + grouping + minfrac filtering are not meaningful.
# Codifies the lesson from the 2026-04-28 SUSHI submission failure where the
# UI-selected subset shrank to 1 file after polarity + sample-type filtering
# and XCMS aborted with "Need >=2 input files for feature detection".
MIN_FILES_FOR_REAL_MODE <- 2L

use_demo <- switch(opt$demo_mode,
  "always" = TRUE,
  "never"  = FALSE,
  "auto"   = (n_present < MIN_FILES_FOR_REAL_MODE)
)

if (opt$demo_mode == "never" && n_present < MIN_FILES_FOR_REAL_MODE) {
  stop("demo_mode=never but only ", n_present, " input file(s) survived",
       " polarity + sample-type filtering. Feature detection needs at least ",
       MIN_FILES_FOR_REAL_MODE, ". Either select more samples in the SUSHI UI,",
       " loosen --sample_types, or set --demo_mode auto/always.")
}
if (opt$demo_mode == "never" && n_present < n_total) {
  stop("demo_mode=never and ", n_total - n_present, " files missing — aborting.")
}

if (opt$demo_mode == "auto" && n_present > 0 && n_present < MIN_FILES_FOR_REAL_MODE) {
  log_msg("[auto-demo] only ", n_present, " input file(s) survived filtering",
          " (need >= ", MIN_FILES_FOR_REAL_MODE, " for real engine) —",
          " falling back to SYNTHETIC DEMO so the SUSHI job still",
          " produces a complete result bundle.")
}

log_msg("[mode] ",
        if (use_demo) "SYNTHETIC DEMO (no real MS processing)"
        else          paste0("REAL processing on ", n_present, " files"))

# ------------------------------------------------------------------
# 4. Build sample_meta
# ------------------------------------------------------------------
sample_meta <- data.frame(
  sample_id    = ds[[c_name]],
  filename     = basename(raw_paths),
  sample_type  = ds[[c_smptype]],
  polarity     = ds[[c_polarity]],
  stringsAsFactors = FALSE
)
sample_meta$group         <- if (!is.na(c_grouping))   ds[[c_grouping]]   else "all"
sample_meta$sample_name   <- if (!is.na(c_sampname))   ds[[c_sampname]]   else sample_meta$sample_id
sample_meta$sample_bfabric<- if (!is.na(c_sampid))     ds[[c_sampid]]     else ""
sample_meta$resource      <- if (!is.na(c_resource))   ds[[c_resource]]   else ""
sample_meta$container     <- if (!is.na(c_container))  ds[[c_container]]  else ""
sample_meta$injection_order <- if (!is.na(c_injord))   suppressWarnings(as.integer(ds[[c_injord]])) else NA_integer_
sample_meta$batch           <- 1L  # single-batch assumption for PoC

# ------------------------------------------------------------------
# 5. Run the engine (or synthesize)
# ------------------------------------------------------------------
engine_params <- list(
  ppm       = opt$mz_tol_ppm,
  peakwidth = as.numeric(strsplit(opt$rt_peakwidth_sec, ",", fixed = TRUE)[[1]]),
  snthresh  = opt$snthresh,
  noise     = opt$noise,
  minfrac   = opt$min_frac
)

if (use_demo) {
  # Synthetic matrix for wiring validation.
  # 40 features × actual sample count; lipid-like m/z in 700-900 range,
  # RT in 1-15 min, log-normal intensities with a spike of "differential"
  # features between the first two groups present in sample_meta.
  set.seed(42)
  n_feat    <- 40
  n_samp    <- nrow(sample_meta)
  mz_vals   <- round(runif(n_feat, 700, 900), 4)
  rt_vals   <- round(runif(n_feat, 1, 15), 3)
  feat_ids  <- sprintf("F%04d_%s", seq_len(n_feat), formatC(mz_vals, format = "f", digits = 2))

  base_intensity <- matrix(
    rlnorm(n_samp * n_feat, meanlog = 12, sdlog = 1.2),
    nrow = n_samp, ncol = n_feat,
    dimnames = list(sample_meta$sample_id, feat_ids)
  )
  # Inject a mild group effect on the first 5 features (for the QC plot in the report)
  distinct_groups <- unique(sample_meta$group)
  if (length(distinct_groups) >= 2) {
    g1_rows <- sample_meta$group == distinct_groups[1]
    base_intensity[g1_rows, 1:5] <- base_intensity[g1_rows, 1:5] * 2.5
  }
  feature_table <- base_intensity

  feature_meta <- data.frame(
    feature_id = feat_ids,
    mz         = mz_vals,
    rt         = rt_vals,
    stringsAsFactors = FALSE
  )

  engine_used <- paste0(opt$engine, " (DEMO_SYNTHETIC)")
  log_msg("[demo] generated synthetic matrix: ", nrow(feature_table),
          " samples × ", ncol(feature_table), " features")

} else {
  # Drop rows with missing files
  keep <- exists_vec
  raw_paths   <- raw_paths[keep]
  sample_meta <- sample_meta[keep, , drop = FALSE]
  if (length(raw_paths) < 2) {
    stop("Need >=2 input files for feature detection; found ", length(raw_paths))
  }

  # MS-DIAL/XCMS/MassCube expect mzML/mzXML (or vendor formats on Windows).
  # For XCMS (Bioconductor) we assume .mzML; .raw would need msconvert first.
  ext <- tolower(tools::file_ext(raw_paths[1]))
  if (opt$engine == "xcms" && ext == "raw") {
    log_msg("[WARN] XCMS needs .mzML input; you passed .raw files.",
            " Consider ProteoWizard msconvert or a pre-step SUSHIApp.")
  }

  log_msg("[engine] run_engine(engine='", opt$engine,
          "', n_files=", length(raw_paths), ")")
  t0 <- Sys.time()
  result <- run_engine(opt$engine, raw_paths, opt$output, params = engine_params)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_msg("[engine] completed in ", round(elapsed, 1), " sec")

  feature_table <- result$feature_table
  feature_meta  <- result$feature_meta
  # Merge engine-side sample_meta (filename, batch) with user-side metadata
  engine_used <- result$engine
}

# ------------------------------------------------------------------
# 6. Lipid annotation (optional post-processing)
# ------------------------------------------------------------------
if (opt$annotation == "goslin" && exists("goslin_standardize", mode = "function")) {
  try({
    lipid_ids <- feature_meta$feature_id
    # demo path: no real lipid names yet; skip unless feature_id is already Liebisch
    log_msg("[annotation] goslin requested but feature IDs are feature_id placeholders; skipped.")
  }, silent = FALSE)
} else if (opt$annotation == "lipidoracle" && exists("lipidoracle_annotate", mode = "function")) {
  log_msg("[annotation] lipidoracle requested — requires .mgf input from MASSter/MZmine; skipped in this PoC.")
}

# ------------------------------------------------------------------
# 7. Write outputs
# ------------------------------------------------------------------
# Matrix with sample ID as the first column (keeps lipid/feature names verbatim
# via check.names=FALSE semantics on the read side).
out_mat <- data.frame(
  Sample = rownames(feature_table),
  feature_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(out_mat,      file.path(opt$output, "feature_table.csv"), row.names = FALSE)
write.csv(feature_meta, file.path(opt$output, "feature_meta.csv"),  row.names = FALSE)
write.csv(sample_meta,  file.path(opt$output, "sample_meta.csv"),   row.names = FALSE)

params_json <- list(
  tool         = "LitmusProcessApp",
  version      = "0.1-poc",
  timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  dataset_tsv  = normalizePath(opt$dataset),
  output_dir   = normalizePath(opt$output),
  engine       = engine_used,
  demo_mode    = use_demo,
  filters      = list(polarity = opt$polarity, sample_types = allowed_types),
  engine_params = engine_params,
  annotation   = opt$annotation,
  n_samples    = nrow(feature_table),
  n_features   = ncol(feature_table),
  groups       = sort(unique(sample_meta$group)),
  litmus_path  = normalizePath(opt$litmus_path),
  sessionInfo  = capture.output(sessionInfo())
)
write(toJSON(params_json, pretty = TRUE, auto_unbox = TRUE),
      file.path(opt$output, "parameters.json"))
log_msg("[OK] wrote feature_table.csv, feature_meta.csv, sample_meta.csv, parameters.json")

# ------------------------------------------------------------------
# 8. Render Rmd report
# ------------------------------------------------------------------
args_raw   <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args_raw, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
rmd_src <- file.path(script_dir, "litmus_process_report.Rmd")

if (file.exists(rmd_src)) {
  suppressPackageStartupMessages(library(rmarkdown))
  rmd_target <- file.path(opt$output, "litmus_process_report.Rmd")
  file.copy(rmd_src, rmd_target, overwrite = TRUE)
  render(rmd_target,
         output_file = "00index.html",
         output_dir  = opt$output,
         params = list(
           params_json   = file.path(opt$output, "parameters.json"),
           feature_table = file.path(opt$output, "feature_table.csv"),
           feature_meta  = file.path(opt$output, "feature_meta.csv"),
           sample_meta   = file.path(opt$output, "sample_meta.csv")
         ),
         quiet = TRUE)
  log_msg("[OK] rendered ", file.path(opt$output, "00index.html"))
} else {
  log_msg("[WARN] Rmd template not found at ", rmd_src, " — skipping HTML report")
}

log_msg("=== DONE ===")
close(log_con)
