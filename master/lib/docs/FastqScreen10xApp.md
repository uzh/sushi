# FastqScreen10xApp

Contamination and composition screening of 10x Genomics single-cell libraries.
ezRun backend: `EzAppFastqScreen_10x` (`R/app-fastqscreen_10x.R`), which runs the
same screens as `EzAppFastqScreen` (`R/app-fastqscreen.R`).

## Scope

- **Analysis type:** screening of a subsample of reads per library to see where they
  come from: adapters, host and model-organism genomes, rRNA, common contaminants,
  other species, and (for human samples) pathogenic viruses. Also estimates the
  strandedness of rRNA reads.
- **Data:** 10x libraries delivered as one tar archive per sample (`RawDataDir`).
- **Output:** a screening report only. The reads themselves go to later apps as they
  were.
- **Not covered here:** bulk FASTQ data (`FastqScreenApp`); read-quality control of
  10x libraries (`Fastqc10xApp`).
- **Place in a pipeline:** a QC step before CellRanger. No other app takes its output
  as input.

## What it can do

1. **Choice of read.** From each tar archive, only the first file of the read that
   carries the cDNA insert is extracted: the first `_R3_` file if the archive has
   one, otherwise the first `_R2_` file. The barcode and UMI read (R1) is not
   screened. Libraries whose archive holds several files per read (for example one
   per lane) are screened from that first file only.
2. **Subsampling.** Each library is reduced to `nReads` reads with ShortRead, using a
   seed fixed in the code.
3. **Preprocessing with fastp.** Adapter trimming is always on, with an in-house
   Illumina adapter list. Optional front and tail trimming, sliding-window quality
   trimming, an average-quality filter, length capping, poly-X trimming and a minimum
   read length come from the form.
4. **Adapter screen.** FastQ Screen with Bowtie 2 aligns the subsampled reads, before
   trimming, against an adapter panel.
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
7. **Species screen against RefSeq mRNA.** The trimmed reads of all libraries are
   pooled and aligned with Bowtie 2 against the RefSeq mRNA collection. Per read, the
   best-scoring alignments at or above `minAlignmentScore` are kept, and each read is
   assigned to a species. The top `nTopSpecies` species per library are reported.
8. **rRNA strandedness.** The trimmed reads are aligned with Bowtie 2 against SILVA
   LSU and SSU rRNA. The best alignments at or above `minAlignmentScore` are counted
   as sense or antisense.
9. **Taxonomic classification.** Kraken 2 classifies the trimmed reads against a
   MiniKraken database.
10. **Base composition.** Per-position base frequencies of the trimmed reads are
    shown as sequence logos.

## How to use it

- **Input dataset:** needs `Name`, `RawDataDir` (a `.tar` per sample) and
  `Read Count`.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `nReads` | Reads per library used for all screens. |
  | `readFileToUse` | `Read1` is the value that works here: the 10x step has already placed the insert read there (see Gotchas). |
  | `cmdOptions` | Bowtie 2 options for the virus, RefSeq and rRNA alignments (steps 6-8). Not used by FastQ Screen. |
  | `minAlignmentScore` | Minimum Bowtie 2 alignment score kept in steps 6-8. |
  | `nTopSpecies` | Number of species shown per library. |
  | `trim_front1`, `trim_tail1`, `cut_front`, `cut_tail`, `cut_right` (+ window size, mean quality), `average_qual`, `max_len1`, `max_len2`, `poly_x_min_len`, `length_required`, `cmdOptionsFastp` | fastp preprocessing. A value of 0 switches the option off. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `Html [Link]` (the report, `00index.html`) and `Report [File]` (the
  report folder).
- **Gotchas:**
  - The form's only choice for `readFileToUse` is `Read2`. The 10x step, however,
    puts the insert read into `Read1`, and with `Read2` the shared FastqScreen step
    replaces it with a `Read2` column this dataset does not have. Successful runs
    record `readFileToUse = Read1`.
  - Results describe one file per library, subsampled, not the full library.

## What to describe in the Methods

Every screen works on the same subsample of one read file per library, so the
description starts there.

**Worth describing**

- **The read screened:** the insert read (R2, or R3 where present) of the first read
  file of each library, subsampled to `nReads` reads. The subsampling seed is fixed
  in the code.
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
  - against RefSeq mRNA for the species screen, with the reads of all libraries
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

- **BWA:** loaded, but never called.
- **Picard `FastqToSam` / `MergeSamFiles`, `samtools view` / `bam2fq`, `sort` and
  `join`.** They only pool the reads of all libraries for the RefSeq alignment and
  keep track of which read came from which library.
- **Threads and housekeeping options:** `--threads`, `-p`, `--thread`, `-t`,
  `--no-unal`, `--nohits`, `--outdir`, `--report`, `--gzip-compressed`, fastp's
  `--compression`, the `fastp.json` copy, output redirection, and the
  `export LANG=...` line.
- **fastp options set to 0** (for example `--max_len1 0`): they are switched off.
- **Tar extraction and the base composition logos.**
- **Counts** of reads, hits, species or percentages. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case, non-human sample, placeholders in angle
brackets):

> For each 10x library, <nReads> reads of the insert read from the first read file
> were subsampled and screened. Reads were preprocessed with fastp <version>
> (<switched-on trimming options>), including adapter trimming with an in-house
> adapter list. Adapter content was assessed on the untrimmed reads, and the trimmed
> reads were screened against the human, mouse and Arabidopsis genomes, SILVA release
> <x> rRNA, tRNA, PhiX, lambda, Mycoplasma and UniVec sequences with FastQ Screen
> <version>, using Bowtie 2 <version> as aligner. Trimmed reads were further aligned
> with Bowtie 2 (<options>) to RefSeq mRNA (<release>), pooled across libraries, to
> assign reads to species, and to SILVA LSU/SSU rRNA to estimate rRNA strandedness,
> keeping best-scoring alignments with a score of at least <minAlignmentScore>.
> Taxonomic classification was performed with Kraken 2 <version> against the
> MiniKraken 8 GB database (<date>).
