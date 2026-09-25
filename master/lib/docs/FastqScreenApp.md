# FastqScreenApp

Contamination and composition screening of sequencing reads.
ezRun backend: `EzAppFastqScreen` (`R/app-fastqscreen.R`).

## Scope

- **Analysis type:** screening of a subsample of reads per sample to see where they
  come from: adapters, host and model-organism genomes, rRNA, common contaminants,
  other species, and (for human samples) pathogenic viruses. Also estimates the
  strandedness of rRNA reads.
- **Data:** bulk sequencing reads (FASTQ), single-end or paired-end, from any
  organism and assay.
- **Output:** a screening report only. The reads themselves go to later apps as they
  were.
- **Not covered here:** 10x single-cell libraries delivered as tar archives
  (`FastqScreen10xApp`); general read-quality control (`FastqcApp`); taxonomic
  profiling of a whole metagenomic sample (`KrakenApp`).
- **Place in a pipeline:** a QC step next to FastQC. No other app takes its output as
  input.

## What it can do

All steps work on the same small subsample of each sample.

1. **Subsampling.** Each sample is reduced to `nReads` reads with ShortRead, using a
   seed fixed in the code. Samples listing several FASTQ files are concatenated
   first.
2. **Choice of read.** `readFileToUse` selects `Read1`, `Read2` or `both`. The
   screens themselves always use one read (the first, or the second when `Read2` is
   chosen). With `both`, the second read is also trimmed and shown in the base
   composition plots.
3. **Preprocessing with fastp.** Adapter trimming is always on. It uses an in-house
   Illumina adapter list plus any adapter given in the dataset. Optional front and
   tail trimming, sliding-window quality trimming, an average-quality filter, length
   capping, poly-X trimming and a minimum read length come from the form.
4. **Adapter screen.** FastQ Screen with Bowtie 2 aligns the subsampled reads,
   before trimming, against an adapter panel.
5. **Genome and rRNA screen.** FastQ Screen with Bowtie 2 aligns the trimmed reads
   against a panel of:
   - the human, mouse and Arabidopsis genomes
   - SILVA rRNA of archaea, bacteria, fungi, plants, chloroplasts, mitochondria,
     metazoa and mammals
   - a tRNA database
   - PhiX, lambda phage, Mycoplasma DNA and UniVec vector sequences
6. **Virus screen (conditional).** For samples whose `Species` starts with "Human" or
   "Homo", reads with no hit in step 5 are aligned with Bowtie 2 against RefSeq
   genomes of human pathogenic viruses.
7. **Species screen against RefSeq mRNA.** The trimmed reads of all samples are
   pooled and aligned with Bowtie 2 against the RefSeq mRNA collection. Per read, the
   best-scoring alignments at or above `minAlignmentScore` are kept, and each read is
   assigned to a species. The top `nTopSpecies` species per sample are reported.
8. **rRNA strandedness.** The trimmed reads are aligned with Bowtie 2 against SILVA
   LSU and SSU rRNA. The best alignments at or above `minAlignmentScore` are counted
   as sense or antisense.
9. **Taxonomic classification.** Kraken 2 classifies the trimmed reads against a
   MiniKraken database. The top species per sample are reported.
10. **Base composition.** Per-position base frequencies of the trimmed reads are
    shown as sequence logos.

## How to use it

- **Input dataset:** needs `Name`, `Read1` and `Read Count`; `Read2` for paired-end
  data.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `nReads` | Reads per sample used for all screens. |
  | `readFileToUse` | Which read is screened: `Read1`, `Read2` or `both`. |
  | `cmdOptions` | Bowtie 2 options for the virus, RefSeq and rRNA alignments (steps 6-8). Not used by FastQ Screen. |
  | `minAlignmentScore` | Minimum Bowtie 2 alignment score kept in steps 6-8. |
  | `nTopSpecies` | Number of species shown per sample. |
  | `trim_front1`, `trim_tail1`, `cut_front`, `cut_tail`, `cut_right` (+ window size, mean quality), `average_qual`, `max_len1`, `max_len2`, `poly_x_min_len`, `length_required`, `cmdOptionsFastp` | fastp preprocessing. A value of 0 switches the option off. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Html [Link]` (the report, `00index.html`) and `Report [File]` (the
  report folder).
- **Gotchas:**
  - The virus screen runs only when the first sample's `Species` starts with "Human"
    or "Homo". It can also be forced with `virusCheck=true` in `specialOptions`.
  - Results describe the subsample, not the full read set.

## What to describe in the Methods

Every screen works on the same subsample, so the description starts there.

**Worth describing**

- **The subsample:** `nReads` reads per sample, and which read was screened
  (`readFileToUse`). The subsampling seed is fixed in the code.
- **fastp preprocessing,** with its version and the options it was given that are
  switched on. Adapter trimming is always on, with an in-house adapter list. fastp
  also applies its built-in quality filter, and may trim poly-G tails on data from
  two-colour Illumina instruments.
- **FastQ Screen, with its version and Bowtie 2 as its aligner:**
  - the adapter screen on the untrimmed reads
  - the genome and rRNA screen on the trimmed reads, against the panel listed under
    "What it can do"

  The panel's SILVA release is part of the configuration file name on the
  `fastq_screen` command line.
- **The direct Bowtie 2 alignments,** with the Bowtie 2 version and their options
  in words:
  - against RefSeq mRNA for the species screen, with the reads of all samples
    aligned together
  - against SILVA LSU/SSU rRNA for strandedness

  In both, the best-scoring alignments with a score of at least
  `minAlignmentScore` were kept. `cmdOptions` applies to these alignments only, not
  to FastQ Screen.
- **The virus screen, when it ran,** for human samples: a `bowtie2` alignment against
  the human pathogenic virus genomes.
- **Kraken 2 and its version,** with the MiniKraken database.
- **Databases by name and release,** not by path. The database folders carry both:

  | Folder | Database |
  |---|---|
  | `minikraken_8GB_<date>` | the MiniKraken 8 GB database of that date |
  | `RefSeq/mRNA/<date>` | the RefSeq mRNA release of that date |
  | `Viruses/ncbi/humanPathogenic_<date>` | the NCBI human pathogenic virus set of that date |
  | `Silva/silva/release_<x>` | SILVA release x |

**Not worth describing**

- **Picard `FastqToSam` / `MergeSamFiles`, `samtools view` / `bam2fq`, `sort` and
  `join`.** They only pool the reads of all samples for the RefSeq alignment and keep
  track of which read came from which sample.
- **Threads and housekeeping options:** `--threads`, `-p`, `--thread`, `-t`,
  `--no-unal`, `--nohits`, `--outdir`, `--report`, `--gzip-compressed`, fastp's
  `--compression`, the `fastp.json` copy, and output redirection.
- **fastp options set to 0** (for example `--max_len1 0`): they are switched off.
- **The base composition logos:** presentation only.
- **Counts** of reads, hits, species or percentages. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case, non-human sample, placeholders in angle
brackets):

> For each sample, <nReads> reads of <read> were subsampled and screened. Reads were
> preprocessed with fastp <version> (<switched-on trimming options>), including
> adapter trimming with an in-house adapter list. Adapter content was assessed on
> the untrimmed reads, and the trimmed reads were screened against the human, mouse
> and Arabidopsis genomes, SILVA release <x> rRNA, tRNA, PhiX, lambda, Mycoplasma and
> UniVec sequences with FastQ Screen <version>, using Bowtie 2 <version> as aligner.
> Trimmed reads were further aligned with Bowtie 2 (<options>) to RefSeq mRNA
> (<release>), pooled across samples, to assign reads to species, and to SILVA
> LSU/SSU rRNA to estimate rRNA strandedness, keeping best-scoring alignments with a
> score of at least <minAlignmentScore>. Taxonomic classification was performed with
> Kraken 2 <version> against the MiniKraken 8 GB database (<date>).
