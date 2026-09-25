# MpileupApp

Joint variant calling on aligned reads with bcftools.
ezRun backend: `EzAppMpileup` (`R/app-mpileup.R`), report template
`inst/templates/Mpileup.Rmd`.

## Scope

- **Analysis type:** calling of small variants (by default SNVs only) jointly across
  all samples of a dataset, with bcftools `mpileup`, `call` and `filter`. It gives one
  multi-sample VCF and a report.
- **Data:** BAM files from a genome alignment (DNA-seq, or RNA-seq for expressed
  variants), for any organism with an FGCZ reference genome.
- **Output:** a filtered, bgzipped and indexed VCF, and a report with sample
  clustering by genotype.
- **Not covered here:** GATK-based calling (`GatkDnaHaplotyperApp`,
  `GatkRnaHaplotyperApp` and related apps); somatic calling (`Mutect2App`);
  structural variants (`DellyApp`, `PbsvApp`).
- **Place in a pipeline:** after `STARApp`, `Bowtie2App` or `BWAApp`.

## What it can do

- **BAM preparation per sample:**
  - optional restriction to a genomic region (`region`)
  - read groups added and reads coordinate-sorted with Picard
    `AddOrReplaceReadGroups`
  - duplicates marked (not removed) with Picard `MarkDuplicates`

  bcftools `mpileup` skips reads flagged as duplicates, so marked duplicates do not
  contribute to the calls.
- **Joint calling across all samples:**
  - `bcftools mpileup` against the reference genome, with the options in
    `mpileupOptions`. The default skips indels and adds per-sample and per-strand
    allele depth annotations.
  - `bcftools call` with the options in `callOptions`. The default uses the
    multiallelic caller, keeps alternative alleles, and outputs variant sites only.
  - `bcftools filter` with the expression in `filterOptions`. The default keeps sites
    where every sample has a minimum read depth.
- **Report:**
  - hierarchical clustering of the samples by genotype (Ward's method), with
    genotypes below a fixed minimum read depth treated as missing
  - the variant positions per chromosome

## How to use it

- **Input dataset:** needs `Name`, `BAM`, `BAI`, `refBuild` and `Species`. `paired` is
  taken from the dataset when present.
- **Jobs:** one job for the whole dataset.
- **Parameters:**

  | Parameter | Effect |
  |---|---|
  | `refBuild` | Reference genome. |
  | `region` | Limit calling to a chromosome or a region like `chr1:1000-2000`. Empty: whole genome. |
  | `mpileupOptions` | Options for `bcftools mpileup`. |
  | `callOptions` | Options for `bcftools call`. |
  | `filterOptions` | Options for `bcftools filter`. |
  | `cores`, `ram`, `scratch` | Compute resources. |

- **Outputs:** `VCF [File]` and `TBI [File]` (the filtered VCF and its index),
  `Html [Link]` (the report) and `Report [File]`.

## What to describe in the Methods

**Worth describing**

- **Duplicate marking** with Picard `MarkDuplicates` and its version, and the
  resulting exclusion of duplicates from calling.
- **bcftools and its version** (from the job script's `module add` line), with the
  three steps and their options in words:
  - `mpileup` against the reference genome: for example "indels skipped; allele depth
    annotations added"
  - `call`: for example "multiallelic caller, variant sites only"
  - `filter`: for example "sites kept where every sample had a read depth above 10"

  All samples were called jointly.
- **The reference genome.** The `-f` path follows
  `<organism>/<source>/<genome build>/Sequence/WholeGenomeFasta/genome.fa`, giving the
  organism and the genome build. No annotation is used in calling.
- **The region, when `region` was set.**
- **Sample clustering by genotype** (Ward's method), briefly, as part of the report.

**Not worth describing**

- **Picard `AddOrReplaceReadGroups`:** it only labels and sorts the reads.
- **Housekeeping:** `-Ou`, `--output-type`, `--output`, the `java` memory and
  temporary-folder settings, `TMP_DIR`, `MAX_RECORDS_IN_RAM`, `VERBOSITY`,
  `VALIDATION_STRINGENCY`, the file copies and indexing (`samtools index`, tabix).
- **The plots of variant positions:** presentation only.
- **Counts** of variants or genotypes. These are results.
- **R, ezRun and SUSHI.** They run the tools; the calling itself is done by bcftools.

**Example Methods paragraph** (placeholders in angle brackets):

> Duplicate reads were marked with Picard <version> MarkDuplicates and excluded from
> variant calling. Variants were called jointly across all samples with bcftools
> <version> against the <organism> <genome build> genome: bcftools mpileup
> (<options in words>), bcftools call (<options in words>) and bcftools filter
> (<filter in words>). Samples were clustered by genotype with Ward's method.
