# CountQCApp

Quality control and exploration of count tables.
ezRun backend: `EzAppCountQC` (`R/app-countQC.R`), report template
`inst/templates/CountQC.qmd`.

## Scope

- **Analysis type:** quality control and unsupervised exploration of a count dataset.
  It covers normalisation, detection of expressed genes, sample correlation and
  clustering, clustering of high-variance genes with GO enrichment, and MDS. There is
  no differential expression testing.
- **Data:** count tables from `FeatureCountsApp` (gene level) or `KallistoApp`
  (transcript level), for any organism with an FGCZ reference annotation.
- **Output:** an interactive HTML report, plus a data file for the exploreDE Shiny
  app.
- **Not covered here:** differential expression between groups (`DESeq2App`,
  `EdgeRApp`, `LimmaApp`).
- **Place in a pipeline:** after counting, before differential expression. No other
  app takes its output as input.

## What it can do

All analysis runs in R, in ezRun's own code and the report template.

- **Loading and annotation.** The counts of all samples are loaded, with gene
  annotation. With `transcriptTypes` set, only features of those types are kept.
- **Optional RUV correction (`runRUV`, not on the form).** RUVSeq `RUVs` estimates
  `kRUVFactors` factors of unwanted variation from the replicate groups, and the
  corrected counts replace the raw counts.
- **Normalisation** (`normMethod`), implemented in ezRun:
  - `logMean` (default): scales each sample so that the mean log signal of the
    expressed genes is the same across samples
  - also available: `quantile` (limma), `median`, `vsn` and `none`
- **Detection of expressed genes:** a gene counts as present in a sample when its
  count is above a signal threshold (`sigThresh`).
- **Count statistics:** read counts per sample, the number of genes above the
  threshold, and dominant genes that take a large share of a sample's reads.
- **Sample correlation and clustering:**
  - Pearson correlation between samples, on log2 signal plus
    `backgroundExpression`, over all genes and over the top genes
  - hierarchical clustering on 1 - correlation (Ward's method)
- **Top genes:** the `topGeneSize` most variable genes, or those with the smallest
  F-test p-value between conditions (`selectByFtest`).
- **Clustering of high-variance genes:**
  - row-centred log2 signal, clustered hierarchically (Ward's method) into groups
  - with `runGO`, GO enrichment of each gene group, with GOstats (hypergeometric
    test) for biological process, molecular function and cellular component
- **MDS:** multidimensional scaling of the samples, on the present genes and on the
  top genes, with at least 4 samples.
- **Scatter plots** of signal between conditions.
- **Enrichr links:** for the top genes per sample, links that open the Enrichr website
  in the reader's browser.

## How to use it

- **Input dataset:** needs `Name`, `Count`, `Species`, `refBuild`, `featureLevel` and
  `refFeatureFile`. A `[Factor]` column such as `Condition` defines the groups.
  `refBuild`, `refFeatureFile` and `transcriptTypes` are taken from the dataset.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `featureLevel` | `gene` or `isoform`, matching the counts. |
  | `normMethod` | Normalisation method (default `logMean`). |
  | `backgroundExpression` | Value added before log2 transformation, which damps noise from low counts. |
  | `topGeneSize` | Number of top genes for clustering and MDS. |
  | `selectByFtest` | Choose the top genes by F-test between conditions instead of by variance. |
  | `runGO` | GO enrichment of the gene clusters. |
  | `transcriptTypes` | Feature types kept. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Static Report [Link]` (the report), `Live Report [Link]` (exploreDE)
  and `Report [File]`.

## What to describe in the Methods

The analysis settings are the job script's `param[['...']]` lines. The report's
"Settings" section also lists the normalisation, the signal threshold and the log2
offset used.

**Worth describing**

- **The software:** ezRun and its version, which carried out the analysis in R, with
  R and its version. Also GOstats and its version when GO enrichment ran. The
  versions appear in the R session listing at the end of the report and log.
- **Normalisation,** in words: for `logMean`, scaling so that the mean log signal of
  expressed genes is equal across samples.
- **The expression threshold:** genes counted as present above a count of
  `sigThresh`.
- **Log2 transformation** with the offset `backgroundExpression`.
- **Sample correlation (Pearson) and hierarchical clustering** (Ward's method, on
  1 - correlation).
- **Selection of the top `topGeneSize` genes,** by variance or by F-test between
  conditions.
- **Clustering of those genes** (Ward's method) and, when `runGO` is true, GO
  enrichment per cluster with GOstats (hypergeometric test; BP, MF, CC).
- **MDS** on present and top genes.
- **RUV correction with RUVSeq, when it was used** (`runRUV`), with the number of
  factors.
- **The annotation and feature types** (`transcriptTypes`), when restricted.

**Not worth describing**

- **The Enrichr links:** they only open the website; no analysis runs.
- **Scatter plots, count statistics tables and the interactive report:**
  presentation of the same data.
- **The exploreDE link.**
- **Counts** of reads, genes, clusters or correlations. These are results.
- **SUSHI:** it only launches the job.

**Example Methods paragraph** (placeholders in angle brackets):

> Count data were analysed with ezRun <version> in R <version>. Counts were
> normalised by <method in words>, genes with more than <sigThresh> counts were
> considered expressed, and signals were log2-transformed after adding <offset>.
> Samples were compared by Pearson correlation and hierarchical clustering (Ward's
> method), and by multidimensional scaling. The <n> most variable genes were
> clustered hierarchically, and each cluster was tested for enrichment of Gene
> Ontology terms with GOstats <version>.
