# KallistoApp

Transcript quantification of RNA-seq reads by pseudo-alignment with kallisto.
ezRun backend: `EzAppKallisto` (`R/app-kallisto.R`).

## Scope

- **Analysis type:** estimation of transcript abundances (estimated counts and TPM)
  per sample, by pseudo-aligning reads to a transcriptome, without a genome
  alignment.
- **Data:** bulk RNA-seq reads (FASTQ), single-end or paired-end, from any organism
  with an FGCZ reference annotation, or against a user-supplied transcriptome.
- **Output:** a transcript-level abundance table per sample, plus bootstrap
  estimates.
- **Not covered here:** genome alignment (`STARApp`); counting reads from BAM files
  (`FeatureCountsApp`); single-cell data.
- **Place in a pipeline:** after read QC. Its count tables feed `CountQCApp`,
  `DESeq2App`, `EdgeRApp` and `LimmaApp`.

## What it can do

- **Read preprocessing with fastp:**
  - adapter trimming (on by default)
  - front and tail trimming, sliding-window quality trimming, an average-quality
    filter, length capping, poly-X trimming and a minimum read length

  The trimmed reads are used only for quantification.
- **Transcriptome index.** Built from the transcript sequences of the selected
  annotation, restricted to the transcript types chosen in `transcriptTypes` (only
  protein-coding transcripts by default). The index is built once per annotation,
  transcript-type combination and kallisto version, on first use.
- **Custom transcriptomes:**
  - `transcriptFasta` quantifies against a supplied transcript FASTA instead, for
    example a Trinity assembly.
  - `secondRef` adds extra sequences to the transcriptome.

  Both build an index for the run.
- **Quantification with `kallisto quant`:**
  - strand-specific or unstranded, following `strandMode`
  - paired-end or single-end; single-end runs need a fragment length mean and SD
  - bootstrap estimates with a seed
- **GPU mode (`gpu = 1`):** runs on a GPU node and skips the bootstrap estimates.

## How to use it

- **Input dataset:** needs `Name`, `Read1` and `Species`; `Read2` for paired-end
  data. `paired` is set from the presence of `Read2`, and `strandMode` from the
  dataset when present.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile` | Reference annotation whose transcripts are quantified. |
  | `transcriptTypes` | Transcript types included in the index (default: protein-coding only). |
  | `strandMode` | `both` (unstranded), `sense` or `antisense`. |
  | `paired` | Paired-end quantification. |
  | `fragment-length`, `sd` | Fragment length mean and SD for single-end data. When left at 0, ezRun uses 180 and 50. |
  | `bootstrap-samples`, `seed` | Number of bootstrap estimates and their seed. |
  | `transcriptFasta`, `secondRef` | Custom transcriptome, or extra sequences. |
  | `trimAdapter` and the fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `gpu`, `gpu_feature` | Run on a GPU node, without bootstrap estimates. |
  | `markDuplicates` | Shown on the form, but has no effect in this app. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):**
  - `Count [File]`: the abundance table (`abundance.tsv`: transcript, length,
    effective length, estimated counts, TPM).
  - `bootstrappedCount [File]`: the bootstrap estimates (`abundance.h5`); not in GPU
    mode.
  - `runInfo [File]`: kallisto's run summary (`run_info.json`).
  - `PreprocessingLog [File]`: the fastp log.
  - The output rows carry `featureLevel = isoform`.

## What to describe in the Methods

**Worth describing**

- **fastp preprocessing,** with its version and the options it was given that are
  switched on, including whether adapter trimming was applied. On the fastp command
  line, `--adapter_fasta` means adapter trimming was on and
  `--disable_adapter_trimming` means it was off. fastp also applies its built-in
  quality filter, and may trim poly-G tails on data from two-colour Illumina
  instruments.
- **kallisto and its version,** used to pseudo-align reads and quantify transcript
  abundances.
- **The transcriptome.** The index path on the `kallisto quant` line follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>/Genes/genes_<types>_kallistoIndex_v<version>/transcripts.idx`:
  - It gives the organism, the annotation source (e.g. Ensembl, GENCODE), the genome
    build and the annotation release.
  - The `-<date>` is when FGCZ set the reference up; it is not part of the release.
  - `<types>` lists the transcript types included, for example `protein_coding`.
  - A `kallistoIndex.../transcripts` or `Custom_kallistoIndex...` folder in the job's
    working directory means a supplied transcriptome or extra sequences.
- **The quantification settings on the `kallisto quant` line, in words:**
  - strandedness (`--fr-stranded` or `--rf-stranded`; neither means unstranded)
  - single-end mode (`--single`) with the fragment length mean and SD
  - the number of bootstrap samples and the seed
- **Index building, when it happened:** a `kallisto index` line.

**Not worth describing**

- **Housekeeping:** `-o`, `-t` (threads), the `export HDF5_DISABLE_VERSION_CHECK=1`
  prefix, and the file renames afterwards.
- **`markDuplicates`:** it has no effect in this app.
- **samtools:** loaded, but not called.
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **Counts** of reads, pseudo-aligned reads or transcripts. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (paired-end base case, placeholders in angle brackets):

> Reads were preprocessed with fastp <version> (<switched-on options>), including
> adapter trimming. Transcript abundances were quantified with kallisto <version>
> against the <transcript types> transcripts of the <source> annotation release
> <release> for <organism> <genome build>, treating the library as <stranded
> orientation>, with <n> bootstrap samples (seed <seed>).
