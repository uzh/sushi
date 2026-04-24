# LITMUS Normalize — SUSHI PoC

**Branch**: `litmus-sushi-app`
**Instance**: `fgcz-h-083:/srv/sushi/masa_test_sushi_20260416/`
**Created**: 2026-04-23

First proof-of-concept converting a LITMUS Shiny workbench module into a
SUSHI batch App. Validates that LITMUS's matrix-in/matrix-out R library
(`normalize_lipidomics.R`) can be wrapped as a `SushiFabric::SushiApp`
without modifying the LITMUS core.

## Files

### In-repo (branch `litmus-sushi-app`)
| Path | Role |
|---|---|
| `master/lib/LitmusNormalizeApp.rb` | SUSHI Ruby frontend (≈100 lines) |
| `master/lib/litmus_normalize.R` | Standalone Rscript wrapping normalize_lipidomics.R |
| `master/lib/litmus_normalize_report.Rmd` | HTML report template (rendered via rmarkdown) |
| `test_data/` | Local reference copies for off-cluster testing only |

### SLURM-accessible staging (for real job submissions)
Data and parameter files for SLURM jobs live on the shared analysis filesystem:

| Path | Role |
|---|---|
| `/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423/lipidomics_demo.csv` | 30 samples × 37 lipids (LITMUS demo data) |
| `/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423/groups_demo.csv` | Sample–group metadata |
| `/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423/litmus_dataset.tsv` | SUSHI Dataset TSV (absolute paths to the CSVs above) |
| `/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423/litmus_params.tsv` | SUSHI parameterset TSV (PQN + log2 + pareto) |

`/srv/sushi/` (for the R scripts and LITMUS core) is already on the shared
filesystem and visible from compute nodes — no staging needed for code.

## Architecture

```
SUSHI Web UI / CLI
  ↓
LitmusNormalizeApp (Ruby, DATASET mode, Dev/R module)
  ↓ commands → sbatch script
Rscript --vanilla litmus_normalize.R \
  --input <intensity CSV> \
  --output <result dir> \
  --method PQN --transform log2 --scaling pareto \
  --litmus_path /srv/sushi/LITMUS_mcp_server
  ↓
source(normalize_lipidomics.R)  # LITMUS core, unchanged
  ↓
normalize_pqn → transform_log2 → scale_pareto
  ↓
normalized_matrix.csv + qc_metrics.csv + 00index.html
```

**Zero modification to LITMUS** — the R script just sources the existing
`scripts/R/normalize_lipidomics.R` and calls its public functions.

## Validate (CLI test mode)

```bash
cd /srv/sushi/masa_test_sushi_20260416/master

TARGET=/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423

bundle exec sushi_fabric \
  --class LitmusNormalizeApp \
  --dataset "$TARGET/litmus_dataset.tsv" \
  --parameterset "$TARGET/litmus_params.tsv" \
  --project 35611
```

Expected: all SUSHI checks PASS, generated command printed.

## Submit (real run via SLURM)

Add `--run` to the same command:

```bash
bundle exec sushi_fabric \
  --class LitmusNormalizeApp \
  --dataset "$TARGET/litmus_dataset.tsv" \
  --parameterset "$TARGET/litmus_params.tsv" \
  --project 35611 \
  --run
```

This submits an sbatch job to SLURM. Output appears under
`p35611/LitmusNormalize_<id>_<timestamp>/LitmusNormalize_Result/`.

## End-to-end smoke test (direct R script)

```bash
mkdir -p /tmp/litmus_smoke
Rscript --vanilla /srv/sushi/masa_test_sushi_20260416/master/lib/litmus_normalize.R \
  --input /srv/sushi/masa_test_sushi_20260416/test_data/lipidomics_demo.csv \
  --output /tmp/litmus_smoke \
  --method PQN \
  --transform log2 \
  --scaling pareto \
  --litmus_path /srv/sushi/LITMUS_mcp_server
```

Expected outputs in `/tmp/litmus_smoke/`:
- `normalized_matrix.csv` (30 samples × 37 features, Pareto-scaled)
- `qc_metrics.csv` + `qc_metrics_features.csv` + `qc_metrics_samples.csv`
- `parameters.json` (provenance)
- `00index.html` (~1 MB self-contained HTML report)

## SUSHI web UI submission (real run)

1. Log into test SUSHI at `http://fgcz-h-083:5000` (via VPN)
2. Navigate to Project `p35611`
3. Register a new DataSet from `/srv/GT/analysis/masaomi/2026/FGCZ/litmus_sushiapp_test_20260423/litmus_dataset.tsv`
4. Run application → `LitmusNormalize` (category: Lipidomics)
5. Select normalization method, transform, scaling; submit

Job goes via SLURM → output appears in project gstore.

## Gotchas discovered

1. **`required_columns` strips type suffix**: `@required_columns = ['Intensity']`
   matches `Intensity [File]` column in the TSV. Do NOT include `[File]`.
2. **Parameterset TSV has NO header row** (SUSHI's `CSV.readlines` reads every
   row as data; header row `Parameter\tValue` would trigger `eval("Value")` →
   NameError).
3. **String values in parameterset TSV should NOT be quoted** — `data_type(k) == String`
   skips eval and preserves the verbatim value including any quote chars.
4. **Array-typed params** (dropdowns): supply just the selected value
   (e.g. `normalization_method\tPQN`), not the quoted literal.
5. **Standalone Rscript bypasses ezRun** for PoC; once validated, move the
   body into `~/git/ezRun/R/app-LitmusNormalize.R` as `ezMethodLitmusNormalize`
   and replace `commands` with `run_RApp("EzAppLitmusNormalize")`.

## Next steps

- [ ] Deploy to real SUSHI test instance and run --run submission
- [ ] Register the `.rb` + `.R` + `.Rmd` triple as a reusable pattern in
      LITMUS `.kairos/knowledge/` as `litmus_sushi_conversion_pattern`
- [ ] Apply the pattern to other LITMUS apps: Statistics, Processing,
      Reporting, FAIR Export
- [ ] Wire `sushi_record_provenance` + `sushi_generate_scaffold` for
      blockchain-recorded reproducibility bundles
