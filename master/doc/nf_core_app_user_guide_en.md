# nf-core SUSHI App User Guide

Last updated: 2026-02-06

---

## Table of Contents

1. [Overview](#1-overview)
2. [Basic Usage](#2-basic-usage)
3. [Available Pipelines](#3-available-pipelines)
4. [Parameter Configuration](#4-parameter-configuration)
5. [Reference Genome Selection](#5-reference-genome-selection)
6. [Project Defaults](#6-project-defaults)
7. [Pipeline-Specific Usage](#7-pipeline-specific-usage)
8. [YAML Configuration for Customization](#8-yaml-configuration-for-customization)
9. [Troubleshooting](#9-troubleshooting)
10. [Known Limitations](#10-known-limitations)

---

## 1. Overview

SUSHI supports nf-core pipelines (https://nf-co.re/) as native SUSHI Apps. nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow.

### Key Features

- **Automatic Registration**: All nf-core pipelines are automatically registered as SUSHI Apps
- **Zero Configuration**: Most pipelines work out of the box (Convention over Configuration)
- **Reference Genome Integration**: Automatic integration with SUSHI's reference genome selector
- **Parameter Auto-Detection**: Required and optional parameters are automatically extracted from pipeline schemas
- **Nextflow Version Auto-Compatibility**: Pipeline versions compatible with the installed Nextflow are automatically selected

### Requirements

- Nextflow installed on the compute cluster
- Singularity/Apptainer available for container execution
- Network access to nf-core API

---

## 2. Basic Usage

### Step 1: Select the App

From the SUSHI Application list, choose an application that starts with `NfCore`.

Examples:
- `NfCoreRnaseqApp` - RNA-seq analysis
- `NfCoreChipseqApp` - ChIP-seq analysis
- `NfCoreSarekApp` - Variant calling
- `NfCoreFetchngsApp` - Download data from public databases

### Step 2: Select a Dataset

Select your input dataset. Most nf-core pipelines require datasets with these columns:

| Column | Description | Required |
|--------|-------------|----------|
| Name | Sample name | Yes |
| Read1 [File] | FASTQ R1 file path | Pipeline-dependent |
| Read2 [File] | FASTQ R2 file path | For paired-end |

### Step 3: Configure Parameters

Set pipeline-specific parameters. See "4. Parameter Configuration" for details.

### Step 4: Submit the Job

Review parameters and click "Submit" to run the job.

---

## 3. Available Pipelines

All pipelines registered on nf-core are automatically available as SUSHI Apps. Main categories:

| Category | Representative Pipelines | Description |
|----------|-------------------------|-------------|
| Transcriptomics | rnaseq, smrnaseq | RNA-seq analysis |
| Genomics | sarek | Variant calling |
| Epigenetics | chipseq, atacseq, methylseq | Epigenomic analysis |
| SingleCell | scrnaseq | Single-cell analysis |
| Metagenomics | ampliseq, taxprofiler | Metagenomic analysis |
| DataFetch | fetchngs | Download from public DBs |

See the full list at https://nf-co.re/pipelines.

---

## 4. Parameter Configuration

### Auto-Managed Parameters (handled by SUSHI)

The following parameters are automatically managed by SUSHI and do not need to be set:

- `--input` : Samplesheet path
- `--outdir` : Output directory
- `--email` : Notification email
- `--publish_dir_mode` : Output file publishing mode

### Resource Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| cores | CPU cores | Varies by category (8–16) |
| ram | Memory (GB) | Varies by category (30–128) |
| scratch | Temp disk (GB) | Varies by category (100–500) |

### Category-Based Default Resources

| Category | cores | ram (GB) | scratch (GB) |
|----------|-------|----------|--------------|
| Transcriptomics | 8 | 64 | 200 |
| Genomics | 8 | 64 | 200 |
| SingleCell | 16 | 128 | 500 |
| Epigenetics | 8 | 30 | 100 |
| Metagenomics | 8 | 64 | 200 |
| QC | 4 | 16 | 100 |
| DataFetch | 4 | 16 | 500 |

### Pipeline-Specific Parameters

Parameters from the pipeline's `nextflow_schema.json` are automatically displayed in the UI:

- Parameters with `[REQUIRED from API]` prefix: **Required** parameters extracted from the nf-core schema
- Parameters without prefix: **Optional** parameters from the nf-core schema. These can be omitted but are easily configurable from the UI

> **Note**: Parameters from `institutional` or `generic` sections and those with the `hidden` attribute are not shown in the UI.

---

## 5. Reference Genome Selection

For pipelines that require a reference genome (rnaseq, chipseq, sarek, etc.), SUSHI's standard reference selector (`refBuild`) is automatically displayed.

### How It Works

1. The system checks the pipeline's `nextflow_schema.json` for a `reference_genome_options` section
2. If fasta, gtf, or genome parameters are found, the `refBuild` selector appears
3. FASTA and GTF paths are automatically resolved from the selected reference

### Notes

- Only references installed on the SUSHI server are available for selection
- Even if a pipeline expects a `--genome` parameter, SUSHI passes `--fasta` and `--gtf` directly

---

## 6. Project Defaults

### Saving Project Defaults

1. Configure parameters on the parameter screen
2. Check "Save these parameters as project defaults" before submitting
3. Next time you use the same app with the same project, saved parameters are applied as defaults

### Notes

- Dropdown (SELECT) components remain as dropdowns after loading project defaults
- Boolean values (true/false) are correctly converted
- Project defaults are saved per app name and project combination

---

## 7. Pipeline-Specific Usage

### NfCoreRnaseqApp (RNA-seq)

**Required Dataset Columns:**
- `Name` - Sample name
- `Read1 [File]` - FASTQ R1
- `Read2 [File]` - FASTQ R2 (for paired-end)

**Key Parameters:**
- `refBuild` - Reference genome (auto-displayed)
- `pseudo_aligner` - "salmon" recommended (workaround for STAR crash)
- `skip_alignment` - Set to "true" when using Salmon pseudo-aligner

**About strandedness:**
- If your dataset doesn't have a `StrandMode` column, `auto` is set automatically
- nf-core/rnaseq will auto-detect strand information

**Version:** The system automatically selects a version compatible with the installed Nextflow. For example, on Nextflow 24.10.x, it will auto-select a compatible version instead of 3.22.2 (which requires Nextflow 25.04+).

### NfCoreFetchngsApp (Data Retrieval)

**Input Method:** Accession IDs (not samplesheet)

**Key Parameters:**
- `accession_ids` - Enter SRA/ENA/GEO/DDBJ accession IDs separated by commas or spaces

**Note:** No MultiQC report is generated.

### NfCoreChipseqApp (ChIP-seq)

**Required Dataset Columns:**
- `Name` - Sample name
- `Read1 [File]` - FASTQ R1
- `Read2 [File]` - FASTQ R2
- `Antibody` - Antibody name
- `Control` - Control sample name

### NfCoreSarekApp (Variant Calling)

**Required Dataset Columns:**
- `Name` - Sample name
- `Read1 [File]` - FASTQ R1
- `Read2 [File]` - FASTQ R2
- `Patient` - Patient ID
- `Lane` - Sequencing lane
- `Sex` - Sex (optional)
- `Status` - Tumor/normal status (optional)

---

## 8. YAML Configuration for Customization

### Overview

Most nf-core pipelines work **without any configuration**, but for special cases you can customize behavior through `config/nf_core_pipelines.yml`.

### When Do You Need YAML?

| Scenario | YAML Needed? | Description |
|----------|--------------|-------------|
| Standard pipeline usage | No | Auto-configuration works |
| Pin a specific version | Yes | Set `default_version` |
| Add custom parameters | Yes | Define `custom_params` |
| Create dropdown selectors | Yes | Use `type: select` |
| Override resources | Yes | Override in `params` |
| Disable MultiQC | Yes | Set `skip_multiqc: true` |

### Basic Syntax

```yaml
pipelines:
  # Pipeline name (use the nf-core pipeline name)
  pipeline_name:
    # Pin version (optional)
    default_version: "3.14.0"
    
    # Override category (optional, usually auto-detected)
    analysis_category: Transcriptomics
    
    # Override input type (optional, usually auto-detected)
    input_type: samplesheet  # samplesheet | id_list | fasta_list
    
    # Explicitly set reference selector (optional, usually auto-detected)
    use_ref_selector: true
    
    # Skip MultiQC report (optional)
    skip_multiqc: false
    
    # Override resources (optional)
    params:
      cores: 16
      ram: 128
      scratch: 500
    
    # Add custom parameters (optional)
    custom_params:
      - name: param_name
        type: select       # select | boolean | text_area | string
        required: true
        options: ["opt1", "opt2", "opt3"]
        default: "opt1"
        description: "Description of the parameter"
```

### custom_params Type Reference

| type | UI Component | Description |
|------|-------------|-------------|
| `select` | Dropdown menu | Specify choices with `options` |
| `boolean` | Checkbox | true/false |
| `text_area` | Multi-line text box | For multi-line text input |
| `string` (default) | Text field | Single-line text input |

### Example 1: Add a Dropdown Parameter

```yaml
pipelines:
  rnaseq:
    custom_params:
      - name: pseudo_aligner
        type: select
        required: false
        options: ["", "salmon"]
        default: "salmon"
        description: "Pseudo-aligner to use"
```

This renders `pseudo_aligner` as a dropdown in the UI, allowing users to select either empty (disabled) or "salmon".

### Example 2: Pin a Pipeline Version

```yaml
pipelines:
  smrnaseq:
    default_version: "2.4.0"
```

### Example 3: Increase Resources

```yaml
pipelines:
  myspecial_pipeline:
    params:
      cores: 32
      ram: 256
      scratch: 1000
```

### Example 4: Add a Required Custom Parameter

```yaml
pipelines:
  smrnaseq:
    custom_params:
      - name: mirtraceSpecies
        type: select
        required: true
        options: ["hsa", "mmu", "rno", "ath"]
        default: "hsa"
        description: "3-letter species code for miRBase"
```

### Example 5: Use a Text Area (for fetchngs accession IDs)

```yaml
pipelines:
  fetchngs:
    input_type: id_list
    skip_multiqc: true
    custom_params:
      - name: accession_ids
        type: text_area
        required: true
        default: ""
        description: "Enter SRA/ENA/GEO/DDBJ accession IDs"
```

### Overriding Samplesheet Mapping

Usually auto-detected, but can be specified manually:

```yaml
pipelines:
  custom_pipeline:
    samplesheet_mapping:
      sample: Name
      fastq_1: "Read1 [File]"
      fastq_2: "Read2 [File]"
      custom_column: "MyCustomColumn"
```

Left side: nf-core samplesheet column name. Right side: SUSHI dataset column name.

---

## 9. Troubleshooting

### Q: Parameter turns into a text box instead of dropdown

**Cause:** In older versions, saving Project Defaults converted dropdowns to text boxes.

**Fix:** This has been fixed in the current version. Clear your browser cache if the issue persists.

### Q: Next DataSet Name is empty

**Cause:** In older versions, the `@name` attribute was not properly set for nf-core apps.

**Fix:** This has been fixed in the current version.

### Q: Strandedness error in RNA-seq

**Error:** `Strandedness must be provided and be one of 'auto', 'forward', 'reverse' or 'unstranded' (NA)`

**Cause:** The dataset was missing a `StrandMode` column.

**Fix:** The current version automatically defaults to `auto` when `StrandMode` is not present.

### Q: How to change the pipeline version

**Method 1:** Change `pipelineVersion` in the SUSHI parameter UI

**Method 2:** Set `default_version` in `config/nf_core_pipelines.yml`

### Q: Pipeline not showing in the list

**Cause:** Possible connection issues with the nf-core API.

**Check:**
1. Inspect the cache file at `tmp/nfcore_cache/pipelines.json`
2. If the cache is stale (>24 hours), delete it and restart the server

### Q: Nextflow version error

**Cause:** Some pipelines (newer versions) require Nextflow 25.04+.

**Fix:** The auto-compatibility feature normally selects a compatible version automatically. If automatic selection causes issues, explicitly set `default_version` in `config/nf_core_pipelines.yml`.

### Q: Will apps auto-update when the server's Nextflow is updated?

**Answer:** Yes. The Nextflow version auto-compatibility feature detects the installed Nextflow version at startup and checks each pipeline's releases from newest to oldest for compatibility. When Nextflow is updated, the next time apps are loaded, they will automatically select newer compatible pipeline versions.

---

## 10. Known Limitations

1. **Nextflow Version Compatibility**: The auto-compatibility feature handles most cases automatically, but if all releases of a pipeline require a newer Nextflow than installed, you may need to set `default_version` manually.

2. **Container Cache**: First-time runs may take longer due to Singularity/Apptainer container downloads.

3. **API Cache**: Pipeline information is cached for 24 hours. Newly published pipelines may take up to 24 hours to appear.

4. **Custom Pipelines**: Non-nf-core custom Nextflow pipelines are not directly supported by this framework.

5. **Dataset Columns**: If required columns for a pipeline are missing from the SUSHI dataset, samplesheet generation will fail. Check the pipeline documentation for required columns.

---

## Appendix: Auto-Mapped Columns

The following nf-core samplesheet columns are automatically mapped to SUSHI dataset columns:

| nf-core Column | SUSHI Column | Notes |
|---------------|-------------|-------|
| sample | Name | Common to all pipelines |
| fastq_1 | Read1 [File] | FASTQ R1 |
| fastq_2 | Read2 [File] | FASTQ R2 |
| fasta | FASTA [File] | FASTA file |
| strandedness | StrandMode | RNA-seq (default: auto) |
| antibody | Antibody | ChIP-seq |
| control | Control | ChIP-seq |
| patient | Patient | Sarek |
| lane | Lane | Sarek |
| sex | Sex | Sarek |
| status | Status | Sarek |
| bam | BAM [File] | BAM file |
| bai | BAI [File] | BAI file |
| condition | Condition | Metadata |
| replicate | Replicate | Metadata |
