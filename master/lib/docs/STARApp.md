# STARApp

Spliced alignment of RNA-seq reads to a reference genome with STAR.
ezRun backend: `EzAppSTAR` (`R/app-mapping.R`, function `ezMethodSTAR`).

## Scope

- **Analysis type:** alignment of RNA-seq reads to a genome, splice-aware, giving one
  sorted BAM file per sample, plus splice junctions, chimeric junctions, a
  strandedness estimate and duplication rates.
- **Data:** bulk RNA-seq reads (FASTQ), single-end or paired-end, from any organism
  with an FGCZ reference genome and annotation. Libraries with UMIs are supported.
- **Not covered here:**
  - Gene-level counting (`FeatureCountsApp`); pseudo-alignment and transcript
    quantification (`KallistoApp`).
  - DNA alignment without splicing (`Bowtie2App`, `BWAApp`).
  - Single-cell data (`STARsoloApp`, the CellRanger apps).
- **Place in a pipeline:** after read QC. Its BAM files feed `FeatureCountsApp`,
  `RnaBamStatsApp` and the splicing and variant apps.

## What it can do

- **Read preprocessing with fastp:**
  - adapter trimming (on by default)
  - front and tail trimming, sliding-window quality trimming, an average-quality
    filter, length capping, poly-X trimming and a minimum read length

  The trimmed reads are used only for the alignment.
- **UMI handling (optional).** With `barcodePattern`, and `barcodePattern2` for
  dual-inline UMIs, umi_tools moves the UMI bases into the read name before
  alignment. After alignment, reads are deduplicated by UMI with umi_tools.
- **Alignment with STAR:**
  - against a genome index that includes the annotation's splice junctions
  - the index is built once per genome and annotation, on first use
  - one-pass mapping, or two-pass mapping to pick up novel junctions (`twopassMode`)
  - FGCZ default alignment options: splice-junction-aware read filtering, a minimum
    number of matched bases, mismatch limits, a cap on multimapping reads with
    random ordering, limits on intron size and mate gap, chimeric (fusion)
    junction detection, and strand information for unstranded data
- **Extra sequences (optional, `secondRef`).** A FASTA file of extra sequences, for
  example transgenes or viruses. With a matching GTF next to it, a custom index of
  genome plus extra sequences is built for the run. Without one, the sequences are
  added at mapping time.
- **BAM handling:**
  - Alignments are sorted and indexed.
  - Duplicates can be flagged in the delivered BAM with Picard (`markDuplicates`).
    They are flagged, not removed.
- **QC of the alignments:**
  - The library strandedness is estimated with RSeQC `infer_experiment.py` on a
    sample of reads.
  - Duplication rates are computed with dupRadar on a Picard-marked copy of the BAM,
    for use by `RnaBamStatsApp`.
- **Extra outputs:** splice junctions (`SJ.out.tab`), chimeric junctions, and an IGV
  session link.

## How to use it

- **Input dataset:** needs `Name`, `Read1` and `Species`; `Read2` for paired-end
  data. `paired` is set from the presence of `Read2`. `refBuild` and `strandMode` are
  taken from the dataset when present.
- **Jobs:** one job per sample.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild`, `refFeatureFile` | Reference genome and annotation. |
  | `paired` | Paired-end alignment. |
  | `strandMode` | The library's strandedness (`both`, `sense`, `antisense`). Not passed to STAR; used for the duplication rates and by downstream apps. |
  | `twopassMode` | Two-pass mapping. |
  | `cmdOptions` | STAR alignment options. The default holds the FGCZ settings listed above. |
  | `secondRef` | Extra sequences, see above. |
  | `trimAdapter` and the fastp parameters | Read preprocessing. A value of 0 switches an option off. |
  | `barcodePattern`, `barcodePattern2` | UMI patterns for umi_tools (`N` = UMI base, `X` = discarded base). |
  | `markDuplicates` | Flag duplicates in the delivered BAM. |
  | `cores`, `ram`, `scratch` | Compute resources. Genomes with very many contigs need much more memory. |

- **Outputs (per sample):**
  - `BAM [File]`, `BAI [File]`: the sorted alignments and their index.
  - `Junctions [File]`, `Chimerics [File]`: splice junctions and chimeric junctions.
  - `StrandFile [Link,File]`: the strandedness estimate.
  - `DupRate [File]`: duplication rates.
  - `STARLog [File]`, `PreprocessingLog [File]`: the STAR and fastp logs.
  - `IGV [Link,File]`: the IGV session link.
- **Gotchas:**
  - The chimeric junction file is empty when the options carry no `--chimSegmentMin`.
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
- **UMI extraction and deduplication, when they ran** (`umi_tools extract` and
  `umi_tools dedup`). Describe the UMI structure in words, for example "an 8-base
  UMI followed by 6 discarded bases on read 2".
- **STAR and its version, and the reference.** The `--genomeDir` path follows
  `<organism>/<source>/<genome build>/Annotation/Release<_><release>-<date>/Genes/genes_STARIndex`:
  - It gives the organism, the annotation source (e.g. Ensembl, GENCODE), the genome
    build and the annotation release. The annotation supplied the splice junctions of
    the index.
  - The `-<date>` is when FGCZ set the reference up; it is not part of the release.
  - A `Custom_STARIndex` folder, or `--genomeFastaFiles` on the STAR line, means the
    genome was extended with extra sequences.
  - A `STAR --runMode genomeGenerate` line means the index was built during the run.
- **Two-pass mapping,** when the STAR line carries `--twopassMode Basic`.
- **The alignment settings on the STAR line, in words:**
  - read filtering: `--outFilterType`, `--outFilterMatchNmin`,
    `--outFilterMismatchNmax`, `--outFilterMismatchNoverLmax`
  - multimapping: `--outFilterMultimapNmax`, `--outSAMmultNmax`,
    `--outMultimapperOrder`
  - splice and gap limits: `--alignSJoverhangMin`, `--alignSJDBoverhangMin`,
    `--alignIntronMax`, `--alignMatesGapMax`
  - read ends: `--alignEndsProtrude`
  - chimeric junction detection: the `--chim...` options
  - strand information: `--outSAMstrandField intronMotif`
- **Duplicate marking of the delivered BAM, when `markDuplicates` is true.** It shows
  as a Picard `MarkDuplicates` line with `O=<sample>.bam`. Duplicates are flagged,
  not removed; the optical duplicate distance is on that line.
- **QC of the alignments,** briefly: strandedness estimated with RSeQC
  `infer_experiment.py` on a sample of reads (its `-s` value), and duplication rates
  computed with dupRadar.

**Not worth describing**

- **STAR housekeeping options:** `--runThreadN`, `--outStd`, `--outSAMtype`,
  `--outSAMattrRGline`, `--outSAMattributes`, `--readFilesCommand` and
  `--genomeLoad`. Also `--sjdbOverhang`, which only restates how the index was built.
- **The dupRadar step's own `MarkDuplicates` line.** It works on a temporary copy
  (`<number>-<sample>.bam` to `..._duprm.bam`) and does not mark the delivered BAM.
- **`strandMode`:** it is not passed to STAR.
- **Sorting and indexing** (`samtools sort` / `index`), the IGV link, file renames,
  and the environment activation for umi_tools.
- **fastp options set to 0,** which are switched off, and fastp's `--thread` and
  `--compression`.
- **Dev/jdk and Dev/Python:** support modules only.
- **Counts** of reads, mapped reads, junctions or percentages. These are results.
- **R, ezRun and SUSHI.** They run the tools; they do not analyse the data.

**Example Methods paragraph** (base case without UMIs, placeholders in angle
brackets):

> Reads were preprocessed with fastp <version> (<switched-on options>), including
> adapter trimming, and aligned to the <organism> <genome build> genome with STAR
> <version>, using a genome index with splice junctions from the <source> annotation
> release <release>. Alignment <in one or two passes> required at least <n> matched
> bases, allowed at most <n> mismatches and <fraction> mismatches per aligned length,
> filtered alignments by annotated and novel splice junctions, and retained reads
> mapping to up to <n> loci. Chimeric junctions were detected with a minimum segment
> length of <n>.
