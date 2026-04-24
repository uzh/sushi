#!/usr/bin/env Rscript
# litmus_normalize.R — PoC standalone wrapper for LITMUS normalize_lipidomics.R
# Invoked by LitmusNormalizeApp.rb (SUSHI Ruby frontend).
#
# Bypasses ezRun for the PoC. Once validated, the same logic can be lifted into
#   ~/git/ezRun/R/app-LitmusNormalize.R  (ezMethodLitmusNormalize)
#
# Arguments (all required except --groups):
#   --input        <path>   Intensity CSV (samples in rows, features in columns)
#   --output       <dir>    Output directory (created if missing)
#   --method       <str>    Normalization method: PQN | TIC | Median | Quantile | IS | none
#   --transform    <str>    Transform: log2 | log10 | none
#   --scaling      <str>    Scaling: pareto | auto | none
#   --litmus_path  <path>   LITMUS repo root (provides normalize_lipidomics.R)
#   --groups       <path>   (Optional) groups CSV (SampleID,Group)
#
# Outputs (in --output):
#   normalized_matrix.csv   Final matrix (samples × features) after the pipeline
#   qc_metrics.csv          Per-sample/feature QC metrics
#   parameters.json         Parameters used (for provenance / reproducibility)
#   00index.html            Rendered Rmd report

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# base-R fallback for rlang's null-default operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# -------------------------------------------------------------------
# 1. CLI argument parsing
# -------------------------------------------------------------------
option_list <- list(
  make_option("--input",       type = "character", help = "Intensity CSV (required)"),
  make_option("--output",      type = "character", help = "Output directory (required)"),
  make_option("--method",      type = "character", default = "PQN",     help = "Normalization method [default: %default]"),
  make_option("--transform",   type = "character", default = "log2",    help = "Transform [default: %default]"),
  make_option("--scaling",     type = "character", default = "pareto",  help = "Scaling [default: %default]"),
  make_option("--litmus_path", type = "character", default = "/srv/sushi/LITMUS_mcp_server",
              help = "LITMUS repo root [default: %default]"),
  make_option("--groups",      type = "character", default = NA_character_, help = "Groups CSV (optional)")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot("--input is required"       = !is.null(opt$input)  && nzchar(opt$input))
stopifnot("--output is required"      = !is.null(opt$output) && nzchar(opt$output))
stopifnot("input file does not exist" = file.exists(opt$input))
stopifnot("litmus_path does not exist" = dir.exists(opt$litmus_path))

dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

cat("=== LITMUS Normalize (SUSHI PoC) ===\n")
cat("  input:       ", opt$input, "\n")
cat("  output:      ", opt$output, "\n")
cat("  method:      ", opt$method, "\n")
cat("  transform:   ", opt$transform, "\n")
cat("  scaling:     ", opt$scaling, "\n")
cat("  groups:      ", opt$groups, "\n")
cat("  litmus_path: ", opt$litmus_path, "\n")

# -------------------------------------------------------------------
# 2. Source LITMUS core library
# -------------------------------------------------------------------
norm_lib <- file.path(opt$litmus_path, "scripts", "R", "normalize_lipidomics.R")
stopifnot("normalize_lipidomics.R not found" = file.exists(norm_lib))
source(norm_lib)
cat("[OK] sourced normalize_lipidomics.R\n")

# -------------------------------------------------------------------
# 3. Read input CSV (samples in rows) — preserve lipid names verbatim
# -------------------------------------------------------------------
raw_df <- read.csv(opt$input, check.names = FALSE, stringsAsFactors = FALSE)

# First column is the sample ID by convention
sample_col <- names(raw_df)[1]
sample_ids <- as.character(raw_df[[sample_col]])

# Detect a Group column if present
group_col_idx <- which(tolower(names(raw_df)) %in% c("group", "groups", "condition", "class"))
groups <- if (length(group_col_idx) > 0) {
  as.character(raw_df[[group_col_idx[1]]])
} else if (!is.na(opt$groups) && file.exists(opt$groups)) {
  gdf <- read.csv(opt$groups, check.names = FALSE, stringsAsFactors = FALSE)
  # Assume columns: SampleID (or similar), Group
  id_col  <- which(tolower(names(gdf)) %in% c("sampleid", "sample", "sample_id", "name"))[1]
  grp_col <- which(tolower(names(gdf)) %in% c("group", "groups", "condition", "class"))[1]
  stopifnot("Groups CSV must have id + group columns" = !is.na(id_col) && !is.na(grp_col))
  gdf[[grp_col]][match(sample_ids, gdf[[id_col]])]
} else {
  rep("all", length(sample_ids))
}

# Drop sample id and group columns to build the numeric matrix
non_numeric_cols <- c(sample_col, names(raw_df)[group_col_idx])
mat_df <- raw_df[, !names(raw_df) %in% non_numeric_cols, drop = FALSE]

# Coerce to numeric matrix (samples × features)
mat <- as.matrix(mat_df)
storage.mode(mat) <- "numeric"
rownames(mat) <- sample_ids

cat("[OK] loaded matrix:", nrow(mat), "samples ×", ncol(mat), "features\n")
cat("     groups found: ", paste(unique(groups), collapse = ", "), "\n")

# -------------------------------------------------------------------
# 4. Run the normalize → transform → scale pipeline
# -------------------------------------------------------------------
# All LITMUS functions are matrix-in/matrix-out, dimensions preserved.
# Pipeline order: raw → normalize → transform → scale (per CLAUDE.md rules)

current_mat   <- mat
current_scale <- "raw"
steps         <- list()

# -- step 1: normalization --
norm_fn <- switch(opt$method,
  "PQN"      = normalize_pqn,
  "TIC"      = normalize_tic,
  "Median"   = normalize_median,
  "Quantile" = normalize_quantile,
  "IS"       = normalize_is,
  "none"     = identity,
  stop("Unknown normalization method: ", opt$method)
)
if (opt$method != "none") {
  current_mat <- norm_fn(current_mat)
  steps[[length(steps) + 1]] <- list(step = "normalize", method = opt$method)
  cat("[OK] normalized via", opt$method, "\n")
}

# -- step 2: transform --
if (opt$transform == "log2") {
  current_mat   <- transform_log2(current_mat)
  current_scale <- "log"
  steps[[length(steps) + 1]] <- list(step = "transform", method = "log2")
  cat("[OK] log2-transformed\n")
} else if (opt$transform == "log10") {
  # log10 equivalent: manual since transform_log2 is the canonical wrapper
  offset <- min(current_mat[current_mat > 0], na.rm = TRUE) / 2
  current_mat <- log10(current_mat + offset)
  current_scale <- "log"
  steps[[length(steps) + 1]] <- list(step = "transform", method = "log10")
  cat("[OK] log10-transformed (offset =", signif(offset, 3), ")\n")
}

# -- step 3: scaling --
if (opt$scaling == "pareto") {
  current_mat   <- scale_pareto(current_mat)
  current_scale <- "centered"
  steps[[length(steps) + 1]] <- list(step = "scale", method = "pareto")
  cat("[OK] Pareto-scaled\n")
} else if (opt$scaling == "auto") {
  current_mat   <- scale_auto(current_mat)
  current_scale <- "centered"
  steps[[length(steps) + 1]] <- list(step = "scale", method = "auto")
  cat("[OK] auto-scaled (z-score)\n")
}

# -------------------------------------------------------------------
# 5. Compute QC metrics (scale-aware)
# -------------------------------------------------------------------
# Per-feature: mean, SD, CV% (only on raw scale), %missing
feature_stats <- data.frame(
  feature  = colnames(current_mat),
  mean     = apply(current_mat, 2, mean, na.rm = TRUE),
  sd       = apply(current_mat, 2, sd,   na.rm = TRUE),
  pct_na   = colMeans(is.na(current_mat)) * 100,
  stringsAsFactors = FALSE
)
if (current_scale == "raw") {
  feature_stats$cv_pct <- 100 * feature_stats$sd / feature_stats$mean
}

# Per-sample: mean, SD
sample_stats <- data.frame(
  sample = rownames(current_mat),
  group  = groups,
  mean   = rowMeans(current_mat, na.rm = TRUE),
  sd     = apply(current_mat, 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------
# 6. Write outputs
# -------------------------------------------------------------------
# 6a. Normalized matrix (restore sample names in first column)
out_norm <- data.frame(Sample = rownames(current_mat), current_mat,
                       check.names = FALSE, stringsAsFactors = FALSE)
write.csv(out_norm, file.path(opt$output, "normalized_matrix.csv"),
          row.names = FALSE)

# 6b. QC metrics
write.csv(feature_stats,
          file.path(opt$output, "qc_metrics_features.csv"),
          row.names = FALSE)
write.csv(sample_stats,
          file.path(opt$output, "qc_metrics_samples.csv"),
          row.names = FALSE)
# Compatibility alias for the Ruby next_dataset output schema
write.csv(feature_stats,
          file.path(opt$output, "qc_metrics.csv"),
          row.names = FALSE)

# 6c. Parameters JSON (for provenance)
params_json <- list(
  tool              = "LitmusNormalizeApp",
  version           = "0.1-poc",
  timestamp         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  input             = normalizePath(opt$input),
  output_dir        = normalizePath(opt$output),
  method            = opt$method,
  transform         = opt$transform,
  scaling           = opt$scaling,
  final_scale       = current_scale,
  pipeline_steps    = steps,
  n_samples         = nrow(current_mat),
  n_features        = ncol(current_mat),
  groups            = unique(groups),
  litmus_path       = normalizePath(opt$litmus_path),
  sessionInfo       = capture.output(sessionInfo())
)
write(toJSON(params_json, pretty = TRUE, auto_unbox = TRUE),
      file.path(opt$output, "parameters.json"))

# -------------------------------------------------------------------
# 7. Render Rmd report (00index.html)
# -------------------------------------------------------------------
# Locate this script's directory so the sibling Rmd can be found
args_raw   <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args_raw, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
rmd_src <- file.path(script_dir, "litmus_normalize_report.Rmd")

if (file.exists(rmd_src)) {
  suppressPackageStartupMessages(library(rmarkdown))
  # Copy Rmd to output so rendering paths resolve
  rmd_target <- file.path(opt$output, "litmus_normalize_report.Rmd")
  file.copy(rmd_src, rmd_target, overwrite = TRUE)
  render(rmd_target,
         output_file = "00index.html",
         output_dir  = opt$output,
         params      = list(
           params_json = file.path(opt$output, "parameters.json"),
           norm_csv    = file.path(opt$output, "normalized_matrix.csv"),
           qc_feats    = file.path(opt$output, "qc_metrics_features.csv"),
           qc_samples  = file.path(opt$output, "qc_metrics_samples.csv"),
           final_scale = current_scale
         ),
         quiet = TRUE)
  cat("[OK] rendered report:", file.path(opt$output, "00index.html"), "\n")
} else {
  cat("[WARN] Rmd template not found at", rmd_src, "— skipping HTML report\n")
}

cat("=== DONE ===\n")
