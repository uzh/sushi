# CellBenderApp

Removal of ambient RNA and empty droplets from 10x count matrices with CellBender.
ezRun backend: `EzAppCellBender` (`R/app-cellBender.R`).

## Scope

- **Analysis type:** background correction of droplet-based single-cell counts.
  CellBender `remove-background` models ambient RNA and empty droplets, and outputs a
  corrected count matrix with the droplets it calls as cells.
- **Data:** raw (unfiltered) 10x count matrices from `CellRangerApp`,
  `CellRangerMultiApp` or `CellRangerARCApp`. For multiome data, only the gene
  expression part is used.
- **Output:** corrected count matrices (filtered and unfiltered) per sample, and
  CellBender's report.
- **Not covered here:** clustering and annotation of the corrected counts
  (`ScSeuratApp`).
- **Place in a pipeline:** between CellRanger and `ScSeuratApp`. Its `CountMatrix`
  feeds `ScSeuratApp`.

## What it can do

- **Input.** The raw count matrix of each sample, taken from `UnfilteredCountMatrix`
  or from the CellRanger result folder. When only a matrix folder exists, it is first
  converted to an h5 file with DropletUtils.
- **Multiome inputs.** ATAC peak features are dropped before inference, so only gene
  counts are modelled.
- **`cellbender remove-background`** on a GPU:
  - It estimates, for every droplet, the probability that it contains a cell, and
    removes ambient RNA counts.
  - It outputs the corrected counts for all droplets, and for the droplets called as
    cells at its target false positive rate.
  - Common adjustments go in `cmdOptions`: `--expected-cells`,
    `--total-droplets-included`, `--fpr`, `--epochs`.
- **Output conversion.** The corrected matrices are rewritten in a format Seurat
  reads.

## How to use it

- **Input dataset:** needs `Name`, `Species`, `refBuild`, `refFeatureFile` and
  `CountMatrix`. The unfiltered matrix is found from `UnfilteredCountMatrix` or
  `ResultDir`.
- **Jobs:** one job per sample, always on a GPU node (the partition and GPU are
  fixed by the app).
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `cmdOptions` | Extra `cellbender remove-background` options, for example `--expected-cells`, `--total-droplets-included`, `--fpr`. |
  | `refBuild`, `refFeatureFile` | Metadata carried to downstream apps; CellBender does not use a reference. |
  | `cores`, `ram`, `scratch` | Compute resources around the GPU job. |

- **Outputs (per sample):**
  - `CountMatrix [Link]`: the corrected counts of the called cells
    (`cellbender_filtered_seurat.h5`).
  - `UnfilteredCountMatrix [Link]`: the corrected counts of all droplets
    (`cellbender_raw_seurat.h5`).
  - `Static Report [Link]`: CellBender's report.
- **Gotchas:**
  - The conda environment is named after one CellBender version (for example
    `gi_cellbender_0.3.2`), but the version that runs is the one CellBender prints.
    The two can differ.

## What to describe in the Methods

CellBender's command runs without an ezRun command-line entry. CellBender prints its
own version and command line in the log, as `cellbender:remove-background: CellBender
<version>` and `cellbender:remove-background: Command:`, together with the false
positive rate it used.

**Worth describing**

- **CellBender and its version, as CellBender prints it.** Not the version in the
  conda environment name.
- **What it did:** removal of ambient RNA and calling of cell-containing droplets
  from the raw 10x count matrix, run on a GPU.
- **The settings given,** in words: the expected number of cells and the total number
  of droplets included, when passed (`--expected-cells`,
  `--total-droplets-included`), and the target false positive rate from the log line
  "Using MCKP noise targets computed for FPR <value>".
- **For multiome data:** only gene expression features were used; ATAC peaks were
  excluded.

**Not worth describing**

- **The derived values CellBender reports,** such as the count priors, the number of
  probable cells and empty droplets, and the largest surely-empty droplet. These are
  results of the run, not settings.
- **`--cuda`, the checkpoint interval and the output file names.**
- **The h5 conversions** (DropletUtils, `ptrepack`) and file deletions.
- **`refBuild` and `refFeatureFile`:** CellBender uses no reference.
- **Counts** of cells or UMIs. These are results.
- **R, ezRun and SUSHI.** They run CellBender; they do not analyse the data.

**Example Methods paragraph** (placeholders in angle brackets):

> Ambient RNA was removed and cell-containing droplets were identified with CellBender
> <version> (remove-background) on the raw count matrices, with <n> expected cells and
> <n> total droplets included, at a target false positive rate of <fpr>.
