# ScSeuratApp

Single-sample single-cell RNA-seq analysis with Seurat: QC, clustering, markers and
cell-type annotation.
ezRun backend: `EzAppScSeurat` (`R/app-ScSeurat.R`), with helpers in
`R/seuratUtils.R`, `R/scTools.R`, `R/sc-estimateAmbient.R` and
`R/cellxgene_annotation.R`.

## Scope

- **Analysis type:** analysis of one single-cell sample per job:
  - cell and gene QC with empty-droplet and doublet removal
  - normalisation, dimensionality reduction and clustering
  - cluster markers
  - cell-type annotation with several optional methods
- **Data:** a 10x-style count matrix per sample, from CellRanger, CellBender, STARsolo
  or BD Rhapsody. Gene expression only.
- **Output:** a Seurat object (`scData.qs2`), a report, and marker and cluster
  tables.
- **Not covered here:** combining several samples (`ScSeuratCombineApp`); ADT, VDJ
  and ATAC layers (`ScMultiOmicsApp`); background correction (`CellBenderApp`).
- **Place in a pipeline:** after CellRanger or `CellBenderApp`. Its Seurat objects
  feed `ScSeuratCombineApp`, `ScMultiOmicsApp` and `SCEVANApp`.

## What it can do

All analysis runs in R, with Seurat and Bioconductor packages.

- **Cell QC:**
  - Per cell: number of UMIs, number of genes, and the percentages of mitochondrial
    and ribosomal-protein reads.
  - A threshold that is set (`nUMI`, `ngenes`, `perc_mito`, `perc_riboprot`) is
    applied as given. A threshold left empty is replaced by outlier detection: cells
    more than a fixed number of median absolute deviations from the median, on a log
    scale for UMIs and genes.
  - Empty droplets are tested with DropletUtils `emptyDrops` against the raw matrix
    and removed, when a raw matrix with extra barcodes exists (not for CellBender
    input).
  - Doublets are scored with scDblFinder (cluster-based) and removed.
- **Gene filtering:** a gene is kept when at least `cellsFraction` of cells have at
  least `geneMinUMI` UMIs for it. An optional list of genes can be excluded.
- **Cell cycle** phases are assigned with scran `cyclone`.
- **Normalisation:** Seurat `NormalizeData` and `SCTransform` (v2). With
  `SCT.regress.CellCycle`, cell cycle scores are regressed out.
- **Dimensionality reduction and clustering:**
  - PCA, then UMAP and a shared nearest-neighbour graph on the first `npcs`
    components
  - Seurat `FindClusters` at each resolution in `resolution`; the first value sets
    the reported clusters
- **Ambient RNA estimate (`estimateAmbient`).** After clustering, SoupX and DecontX
  each estimate the ambient RNA fraction of every cell, stored as cell metadata. The
  counts used for normalisation, clustering and markers are not corrected.
- **Cluster markers:**
  - Seurat `FindAllMarkers`, positive markers only, with the test in `DE.method`
    (`wilcox` or `LR`), and with `min.pct` and `logfc.threshold`
  - the markers are then filtered by the difference in detection rate
    (`min.diff.pct`) and the adjusted p-value (`pvalue_allMarkers`)
- **Cell-type annotation (each optional):**

  | Method | What it does |
  |---|---|
  | AUCell | Scores CellMarker 2.0 gene sets of the tissues in `tissue`. |
  | Enrichr | Tests the cluster markers against the Enrichr libraries in `enrichrDatabase`, through the Enrichr web service. |
  | Azimuth | Maps cells to the Azimuth reference in `Azimuth`. |
  | Azimuth Pan-Human | `AzimuthPanHuman`, human only. Reuses CellRanger's local annotation when the CellRanger output has it; otherwise sends the expression matrix to the external CloudAzimuth API. |
  | SingleR | Labels cells and clusters with the celldex reference in `SingleR`. |
  | scType | Uses its marker database for the tissue in `sctype.tissue`. |
  | CyteTypeR | Sends cluster information to the external CyteType service (API key needed). |
  | mLLMCelltype | Labels each cluster from its top marker genes with the FGCZ-internal language model, using `mLLMCelltype.tissue` as context. Only gene names are sent. |
  | CELLxGENE label transfer | Transfers labels from the CELLxGENE dataset in `cellxgeneUrl` with Seurat `FindTransferAnchors` / `TransferData`, on a curated reference of up to 10 donors and 3000 cells per label per donor. |

- **TF and pathway activity (`computePathwayTFActivity`, human and mouse):**
  transcription factor activities from DoRothEA regulons and pathway activities from
  PROGENy, with decoupleR (weighted mean).

## How to use it

- **Input dataset:** needs `Name`, `Species`, `refBuild`, `CountMatrix`, `ResultDir`
  and `Condition`. An `UnfilteredCountMatrix` enables the empty-droplet test.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `nUMI`, `ngenes`, `perc_mito`, `perc_riboprot` | Fixed QC thresholds; empty means outlier detection. |
  | `cellsFraction`, `geneMinUMI` | Gene filtering. |
  | `filterByExpression` | Keep only cells matching an expression rule. |
  | `SCT.regress.CellCycle` | Regress cell cycle in SCTransform. |
  | `npcs`, `pcGenes` | Number of principal components; genes used for PCA. |
  | `resolution` | Clustering resolutions; the first is the reported one. |
  | `DE.method`, `min.pct`, `min.diff.pct`, `logfc.threshold`, `pvalue_allMarkers` | Marker detection and filtering. |
  | `estimateAmbient` | SoupX and DecontX ambient RNA estimates. |
  | `tissue`, `enrichrDatabase`, `Azimuth`, `AzimuthPanHuman`, `SingleR`, `sctype.enabled`, `sctype.tissue`, `CyteTypeR`, `mLLMCelltype`, `mLLMCelltype.tissue`, `cellxgeneUrl`, `cellxgeneLabel` | Annotation methods, as above. |
  | `computePathwayTFActivity` | TF and pathway activity. |
  | `refBuild`, `refFeatureFile`, `geneCountModel` | Reference, and the count model for STARsolo input. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):** `SC Seurat [Link]` (`scData.qs2`), `Static Report [Link]`
  and `SC Cluster Report [File]`, with the marker table (`posMarkers.xlsx`) and
  cluster table (`clusterInfos.xlsx`). `CountMatrix` and `ResultDir` are passed on
  for downstream apps.

## What to describe in the Methods

The analysis settings are the job script's `param[['...']]` lines. Package versions
appear in the R session listing at the end of the log. scType, mLLMCelltype, Azimuth
Pan-Human and CyteTypeR each write a log line when they complete, are skipped or fail;
mLLMCelltype also logs the model it used. The other steps write no log line. Their
result files in the report folder show that they ran: `cells.AUC.qs2` (AUCell),
`enrichRout.qs2` (Enrichr), `TFActivity.qs2` and `pathwayActivity.qs2` (decoupleR),
`sctype_results.rds`, `mllmcelltype_results.rds`.

**Worth describing**

- **Seurat and its version, in R and its version.**
- **Cell QC:** the thresholds that were set, and outlier detection by median absolute
  deviation for those left empty; empty-droplet removal with DropletUtils
  `emptyDrops`, when it ran; doublet removal with scDblFinder.
- **Gene filtering** (`cellsFraction`, `geneMinUMI`).
- **Cell cycle assignment** with scran, and its regression when
  `SCT.regress.CellCycle` is true.
- **Normalisation** with SCTransform (v2).
- **PCA, UMAP and clustering:** the number of components (`npcs`) and the reported
  resolution.
- **Marker detection:** `FindAllMarkers` with the test (`DE.method`), positive
  markers only, and the thresholds `min.pct`, `logfc.threshold`, `min.diff.pct` and
  `pvalue_allMarkers`.
- **Each annotation method that ran,** with its reference or database:
  - the CellMarker tissues for AUCell
  - the Enrichr libraries
  - the Azimuth reference
  - the SingleR reference
  - the scType tissue
  - for Azimuth Pan-Human, whether CellRanger's local annotation was reused or
    CloudAzimuth was called
  - for mLLMCelltype, the model (from the log line "mLLMCelltype using model ...")
    and the tissue context; its labels are a result of the analysis. A log warning
    that the tissue is 'auto' means no tissue context was given.
  - for CELLxGENE, the source dataset and label column
- **Ambient RNA estimation** with SoupX and DecontX, when it ran, as a per-cell
  estimate of contamination. The analysis itself used the uncorrected counts.
- **TF and pathway activities** with decoupleR (DoRothEA, PROGENy), when they ran.

**Not worth describing**

- **Annotation methods that were disabled or skipped.**
- **The report, plots and the exploreSC link.**
- **Parallelisation and memory settings,** the internal random seed, and the file
  conversions.
- **Counts** of cells, genes, clusters or markers, the doublet numbers scDblFinder
  reports, and the contamination estimates. These are results.
- **ezRun and SUSHI.** ezRun assembles the analysis; its steps are described through
  the points above.

**Example Methods paragraph** (base case, placeholders in angle brackets):

> Single-cell data were analysed with Seurat <version> in R <version>. Cells were
> filtered on UMI counts, detected genes and mitochondrial and ribosomal-protein
> content, using <fixed thresholds | outlier detection by median absolute
> deviation>; empty droplets were removed with emptyDrops (DropletUtils <version>) and
> doublets with scDblFinder <version>. Genes detected in at least <fraction> of cells
> were kept. Data were normalised with SCTransform (v2), and cell cycle phases were
> assigned with scran. Principal component analysis was followed by UMAP and graph-based
> clustering on the first <npcs> components at a resolution of <resolution>. Cluster
> markers were identified with FindAllMarkers (<test>), keeping positive markers with
> <thresholds>. Cell types were annotated with <methods and references>.
