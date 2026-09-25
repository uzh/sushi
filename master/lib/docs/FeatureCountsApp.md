# FeatureCountsApp

Counting of aligned reads per gene (or other feature) with featureCounts.
ezRun backend: `EzAppFeatureCounts` (`R/app-featureCounts.R`).

## Scope

- **Analysis type:** assignment of aligned reads to annotated features and counting
  them, giving one count table per sample. By default, counts are per gene, from
  exon overlaps.
- **Data:** BAM files from a genome alignment, typically RNA-seq aligned with
  `STARApp`.
- **Output:** a count table and an assignment summary per sample.
- **Not covered here:** alignment itself (`STARApp`); transcript quantification
  without alignment (`KallistoApp`).
- **Place in a pipeline:** after `STARApp`. Its count tables feed `CountQCApp`,
  `DESeq2App`, `EdgeRApp` and `LimmaApp`.

## What it can do

- **Counting with featureCounts** from the Rsubread R package. It runs inside R, not
  as a command-line program.
- **Annotation:**
  - the GTF of the selected reference annotation
  - restricted to the transcripts of the types chosen in `transcriptTypes` (only
    protein-coding by default)
  - transcripts without type information, such as spike-ins, are kept
- **Features and levels.** Reads are matched to GTF features of one type
  (`gtfFeatureType`, default `exon`) and summed per gene, transcript, exon or ORF
  (`featureLevel`).
- **Counting settings:**
  - strandedness
  - paired-end fragments counted as one unit
  - multimapping reads counted or not
  - primary alignments only, or fractional counting of all alignments
  - reads overlapping several features counted or not
  - a minimum overlap with a feature
  - a minimum mapping quality
  - duplicates ignored or counted
- **Settings fixed in the code:**
  - For paired-end data, a pair is counted even when only one end is mapped.
  - Fragment length is not checked.
  - Chimeric fragments are counted.
  - Reads are not extended.
- **Extra sequences (`secondRef`).** Their annotation is added to the GTF, so reads
  on extra sequences are counted too.
- **Counting around transcription start sites (hidden option).** With
  `aroundTSSCounting` set through `specialOptions`, reads are counted in windows
  around gene starts instead.
- **Sorting.** Paired-end BAM files are sorted by read name first.

## How to use it

- **Input dataset:** needs `Name`, `BAM`, `BAI` and `refBuild`. `paired`,
  `strandMode` and `refFeatureFile` are taken from the dataset when present.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile` | Reference annotation. |
  | `transcriptTypes` | Transcript types kept in the annotation (default: protein-coding only). |
  | `strandMode` | `both` (unstranded), `sense` or `antisense`. |
  | `paired` | Count read pairs as fragments. |
  | `featureLevel`, `gtfFeatureType` | Level of the counts, and the GTF feature type that defines it. |
  | `allowMultiOverlap` | Count reads overlapping several features. |
  | `keepMultiHits`, `countPrimaryAlignmentsOnly` | Count multimapping reads; primary alignment only, or fractional counts over all alignments. |
  | `minFeatureOverlap` | Minimum number of overlapping bases. |
  | `minMapQuality` | Minimum mapping quality. |
  | `ignoreDup` | Skip reads marked as duplicates. |
  | `secondRef` | Extra sequences whose features are also counted. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):** `Count [File,Link]` (the count table) and
  `Stats [File,Link]` (featureCounts' assignment summary). The rows carry
  `featureLevel`, `strandMode` and `transcriptTypes` for downstream apps.

## What to describe in the Methods

featureCounts has no command line of its own here. Its settings appear in two places:
- The job script's `param[['...']]` lines hold the settings from the form.
- featureCounts prints a "featureCounts setting" box in the log, followed by the
  strandedness it applied.

**Worth describing**

- **featureCounts, from the Rsubread package, with the Rsubread version.** The version
  appears in the R session listing at the end of the log (`Rsubread_<version>`). R
  and its version are worth naming here, because the counting runs in R.
- **The annotation:** the reference annotation, restricted to the chosen
  `transcriptTypes`, and the counting level (`featureLevel`, from `gtfFeatureType`
  features). The `refBuild` value follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>`. The
  `-<date>` is when FGCZ set the reference up; it is not part of the release.
- **The counting settings, in words:**
  - strandedness (the box's "Strand specific" line)
  - paired-end fragments counted as one unit
  - how multimapping reads were handled (`keepMultiHits`,
    `countPrimaryAlignmentsOnly`)
  - whether reads overlapping several features were counted (`allowMultiOverlap`)
  - the minimum overlap (`minFeatureOverlap`)
  - the minimum mapping quality (`minMapQuality`)
  - whether duplicates were ignored (`ignoreDup`)
- **The fixed settings, for paired-end data:** pairs were counted when only one end
  was mapped, and fragment length was not checked.
- **Extra sequences, when `secondRef` was set.**

**Not worth describing**

- **Sorting by read name** (`samtools sort -n`) and the temporary annotation file
  (`<number>-genes.gtf`).
- **Threads and the temporary-file folder.**
- **Counts** of assigned, unassigned or total reads, and percentages. These are
  results.
- **ezRun and SUSHI.** They run the counting; they do not change it.

**Example Methods paragraph** (paired-end base case, placeholders in angle brackets):

> Reads were counted per <gene> with featureCounts from the Rsubread package
> <version> in R <version>, using the <source> annotation release <release> for
> <organism> <genome build>, restricted to <transcript types> transcripts, with exon
> overlaps summarised per gene. Counting was <reversely stranded>, read pairs were
> counted as fragments, multimapping reads were counted using the primary alignment
> only, reads overlapping several genes were counted, and a minimum overlap of <n>
> bases and a minimum mapping quality of <q> were required.
