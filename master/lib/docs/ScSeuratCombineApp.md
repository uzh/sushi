# ScSeuratCombineApp

Joint analysis of several single-cell samples with Seurat: merging, optional
integration, clustering, markers and annotation.
ezRun backend: `EzAppScSeuratCombine` (`R/app-ScSeuratCombine.R`), with helpers in
`R/seuratUtils.R`.

## Scope

- **Analysis type:** combined analysis of the single-sample Seurat objects of a
  dataset:
  - merging, and optional batch integration
  - normalisation, dimensionality reduction and clustering of all cells together
  - cluster markers and cell-type annotation
- **Data:** the outputs of `ScSeuratApp` (one Seurat object per sample, already
  QC-filtered).
- **Output:** a combined Seurat object (`scData.qs2`), a report, and marker and
  cluster tables.
- **Not covered here:** per-sample QC and annotation (`ScSeuratApp`); differential
  expression between conditions within clusters (`ScSeuratCompareApp`).
- **Place in a pipeline:** after `ScSeuratApp`. Its combined object feeds
  `ScSeuratCompareApp` and `ScSeuratCombinedLabelClusters`.

## What it can do

All analysis runs in R, with Seurat and Bioconductor packages.

- **Loading.** Each sample's Seurat object is read, with its per-sample Azimuth and
  CELLxGENE labels when present. Cells keep the QC filtering of `ScSeuratApp`.
  `chosenClusters` (not on the form) can restrict the samples to selected clusters.
- **Uncorrected clustering.** All samples are merged and analysed without
  correction. This is kept for comparison (`umap_noCorrected`, `ident_noCorrected`).
- **Normalisation** (`normalizationMethod`): `SCTransform` (default) or
  `LogNormalize`, with optional regression of the cell cycle
  (`SCT.regress.CellCycle`).
- **Integration** (`integrationMethod`):
  - `Harmony` (default): Harmony on the PCA, grouped by `harmonyGroupBy`
    (`Condition` or `Batch`)
  - `CCA` or `RPCA`: Seurat anchor-based integration (`FindIntegrationAnchors`,
    `IntegrateData`)
  - `none`: the merged data without correction
- **Dimensionality reduction and clustering:** PCA, UMAP and a shared
  nearest-neighbour graph on the first `npcs` components, and Seurat `FindClusters`
  at `resolution`.
- **Cluster markers:** Seurat `FindAllMarkers` on all cells:
  - positive markers only, with the test in `DE.method` (`wilcox` or `LR`)
  - thresholds `min.pct` and `logfc.threshold`
  - covariates in `DE.regress` (`Batch`, `CellCycle`) used as latent variables with
    `LR`
  - with SCTransform, `PrepSCTFindMarkers` is run first
- **Cell-type annotation (optional):** AUCell with CellMarker 2.0 gene sets of the
  tissues in `tissue`; Enrichr (web service) with the libraries in
  `enrichrDatabase`; SingleR with the celldex reference in `SingleR`.
- **TF and pathway activity (`computePathwayTFActivity`):** DoRothEA and PROGENy
  activities with decoupleR.
- **Additional factors** (`additionalFactors`): dataset columns copied onto the cells
  for later comparisons.

## How to use it

- **Input dataset:** the output rows of `ScSeuratApp`: needs `Name`, `Species`,
  `refBuild`, `refFeatureFile` and `Static Report`. The `SC Seurat` column locates
  each object.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `integrationMethod`, `harmonyGroupBy` | Integration method, and the Harmony grouping variable. |
  | `normalizationMethod`, `SCT.regress.CellCycle` | Normalisation, and cell cycle regression. |
  | `npcs`, `pcGenes` | Number of principal components; genes used for PCA. |
  | `resolution` | Clustering resolution. |
  | `DE.method`, `DE.regress`, `min.pct`, `logfc.threshold` | Marker detection. |
  | `tissue`, `enrichrDatabase`, `SingleR` | Annotation methods. |
  | `computePathwayTFActivity` | TF and pathway activity. |
  | `additionalFactors` | Extra dataset columns carried onto cells. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `SC Seurat [Link]` (`scData.qs2`), `Static Report [Link]` and
  `Report [File]`, with `posMarkers.xlsx` and a cluster table.

## What to describe in the Methods

The analysis settings are the job script's `param[['...']]` lines. Package versions
appear in the R session listing at the end of the log. Older runs may lack
parameters added later (for example `normalizationMethod`, `harmonyGroupBy`). Their
method then follows the code of that time; the job script shows which parameters
were set.

**Worth describing**

- **Seurat and its version, in R and its version,** and that the samples processed
  with ScSeurat were combined.
- **The annotation,** by source and release, from `refBuild`, which follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>`. The
  `-<date>` is when FGCZ set the reference up; it is not part of the release.
- **Normalisation** (`SCTransform` or `LogNormalize`), with cell cycle regression
  when set.
- **Integration:** the method, with the Harmony version and its grouping variable
  when Harmony was used; no integration when `none`.
- **PCA, UMAP and clustering:** the number of components (`npcs`) and the resolution.
- **Marker detection:** `FindAllMarkers` with the test (`DE.method`), positive
  markers only, the thresholds, and the latent variables from `DE.regress` when
  used.
- **Each annotation method that ran,** with its reference or database: AUCell with
  CellMarker tissues, Enrichr libraries, SingleR reference.
- **TF and pathway activities,** when they ran.

**Not worth describing**

- **The uncorrected clustering:** it is kept for comparison only.
- **The report, plots and the exploreSC link.**
- **Parallelisation and memory settings, and file conversions.**
- **Counts** of cells, clusters or markers. These are results.
- **ezRun and SUSHI.** ezRun assembles the analysis; its steps are described through
  the points above.

**Example Methods paragraph** (Harmony base case, placeholders in angle brackets):

> The samples were combined and analysed jointly with Seurat <version> in R
> <version>. Data were normalised with SCTransform, and samples were integrated with
> Harmony <version> using <harmonyGroupBy> as the grouping variable. Principal
> component analysis was followed by UMAP and graph-based clustering on the first
> <npcs> components at a resolution of <resolution>. Cluster markers were identified
> with FindAllMarkers (<test>), keeping positive markers with <thresholds>.
