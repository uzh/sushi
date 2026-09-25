# Bowtie2App

Alignment of sequencing reads to a reference genome with Bowtie 2.
ezRun backend: `EzAppBowtie2` (`R/app-mapping.R`, function `ezMethodBowtie2`).

## Scope

- **Analysis type:** alignment of reads to a genome without splicing, giving one
  sorted BAM file per sample.
- **Data:** DNA-based assays such as whole-genome or targeted DNA sequencing,
  ChIP-seq, ATAC-seq and CUT&RUN, and reads against small genomes. Single-end or
  paired-end.
- **Not covered here:** RNA-seq, which needs spliced alignment (`STARApp`); BWA
  alignment (`BWAApp`); single-cell data.
- **Place in a pipeline:** after read QC. Its BAM files feed peak calling
  (`MACS3App`), variant calling (`MpileupApp`, the GATK apps) and `DnaBamStatsApp`.

## What it can do

- **Read preprocessing with fastp:**
  - adapter trimming (on by default)
  - front and tail trimming, sliding-window quality trimming, an average-quality
    filter, length capping, poly-X trimming and a minimum read length

  The trimmed reads are used only for the alignment.
- **Alignment with Bowtie 2:**
  - against the genome index of the selected genome build; no annotation is used
  - the index is built once per genome, on first use
  - further alignment options come from `cmdOptions`
- **Extra sequences (optional, `secondRef`).** A FASTA file of extra sequences, for
  example a transgene or a virus, is added to the genome in a custom index built for
  the run. A coverage plot of each extra sequence is also produced.
- **BAM handling:**
  - Alignments are converted to BAM, sorted and indexed.
  - Duplicates are flagged with Picard MarkDuplicates, not removed (on by default,
    `markDuplicates`).
- **Optional BigWig coverage track** (`generateBigWig`), and an IGV session link.

## How to use it

- **Input dataset:** needs `Name`, `Read1` and `Species`; `Read2` for paired-end
  data. `paired` is set from the presence of `Read2`.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild` | Reference genome. Only the genome is used; the annotation part of the selection has no effect here. |
  | `paired` | Paired-end alignment. |
  | `cmdOptions` | Bowtie 2 options. The default `--no-unal` only leaves unaligned reads out of the BAM. |
  | `secondRef` | Extra sequences, see above. |
  | `trimAdapter` and the fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `markDuplicates` | Flag duplicates with Picard (default on). |
  | `generateBigWig` | Write a BigWig coverage track. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs (per sample):**
  - `BAM [File]`, `BAI [File]`: the sorted alignments and their index.
  - `DupMetrics [File,Link]`: Picard duplicate metrics, when `markDuplicates` is on.
  - `PreprocessingLog [File,Link]`: the fastp log.
  - `IGV [File,Link]`: the IGV session link.
  - With the options on, `BigWig [File]` and `SecondRefCoverage [File,Link]`.
- **Gotchas:**
  - A run waits while another job builds the same index, and gives up after a fixed
    timeout.

## What to describe in the Methods

**Worth describing**

- **fastp preprocessing,** with its version and the options it was given that are
  switched on, including whether adapter trimming was applied. On the fastp command
  line, `--adapter_fasta` means adapter trimming was on and
  `--disable_adapter_trimming` means it was off. fastp also applies its built-in
  quality filter, and may trim poly-G tails on data from two-colour Illumina
  instruments.
- **Bowtie 2 and its version, and the genome.**
  - The `-x` path follows
    `<organism>/<source>/<genome build>/Sequence/BOWTIE2Index/genome`, giving the
    organism and the genome build.
  - Bowtie 2 uses no annotation. The `refBuild` parameter names an annotation
    release, but it plays no part in this step.
  - A `BOWTIE2Index/customGenome` index means the genome was extended with extra
    sequences.
  - A `bowtie2-build` line means the index was built during the run.
- **Alignment options on the `bowtie2` line,** other than the housekeeping ones
  below. When there are none, Bowtie 2 aligned with its built-in settings.
- **Paired-end alignment,** when the line carries `-1` and `-2`.
- **Duplicate marking, when a Picard `MarkDuplicates` line is present.** Duplicates
  are flagged, not removed (`REMOVE_DUPLICATES=false`); the optical duplicate
  distance is on that line.

**Not worth describing**

- **Bowtie 2 housekeeping:** `-p` (threads), the `--rg-id` / `--rg` read-group
  labels, `--no-unal` (unaligned reads left out of the file), and the log
  redirection.
- **The annotation release in `refBuild`.**
- **BAM conversion, sorting and indexing** (`samtools view` / `sort` / `index`).
- **The IGV link, the BigWig track and the extra-sequence coverage plot:**
  presentation only.
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **Counts** of reads, alignment rates or duplicate rates. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case, placeholders in angle brackets):

> Reads were preprocessed with fastp <version> (<switched-on options>), including
> adapter trimming, and aligned <as pairs> to the <organism> <genome build> genome
> with Bowtie 2 <version> <with the options ... | with its built-in settings>.
> Duplicate reads were flagged with Picard <version> MarkDuplicates, using an
> optical duplicate distance of <n> pixels, and were not removed.
