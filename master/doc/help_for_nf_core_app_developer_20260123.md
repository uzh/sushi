# nf-core App Developer Guide for SUSHI

This guide is for developers who maintain and customize nf-core applications within the SUSHI workflow management system.

## Table of Contents

1. [Overview](#overview)
2. [Automatic Pipeline Loading](#automatic-pipeline-loading)
3. [Customizing Pipeline Configuration](#customizing-pipeline-configuration)
4. [Creating a Manual SUSHI App](#creating-a-manual-sushi-app)
5. [Input Types](#input-types)
6. [Parameter Configuration](#parameter-configuration)
7. [Samplesheet Mapping](#samplesheet-mapping)
8. [Nextflow & Apptainer Tips](#nextflow--apptainer-tips)
9. [Troubleshooting](#troubleshooting)

---

## Overview

SUSHI automatically loads all nf-core pipelines from the nf-core API at startup. Each pipeline is converted to a dynamic `SushiApp` class (e.g., `NfCoreRnaseqApp`, `NfCoreFetchngsApp`).

**Key Features:**
- 142+ nf-core pipelines available by default
- Automatic parameter detection from `nextflow_schema.json`
- Customizable via YAML configuration
- Support for different input types (samplesheet, ID list, etc.)

---

## Automatic Pipeline Loading

### How It Works

1. At Rails startup, `NfCoreAppFactory` fetches pipeline metadata from `https://nf-co.re/pipelines.json`
2. For each pipeline, a dynamic Ruby class is created inheriting from `SushiFabric::SushiApp`
3. Pipeline parameters are loaded from `nextflow_schema.json` when the app is selected

### Default Behavior

| Attribute | Default Value |
|-----------|---------------|
| `required_columns` | `['Name']` |
| `analysis_category` | Mapped from nf-core topics or `'nf-core'` |
| `input_type` | Auto-detected from `schema_input.json` |
| `params.cores` | 8 |
| `params.ram` | 30 |
| `params.scratch` | 100 |

---

## Customizing Pipeline Configuration

To customize a pipeline, edit `config/nf_core_pipelines.yml`:

### Basic Configuration

```yaml
pipelines:
  rnaseq:
    nf_core_name: rnaseq
    default_version: "3.14.0"
    analysis_category: Transcriptomics
    params:
      cores: 16
      ram: 64
      scratch: 200
```

### Full Configuration Options

```yaml
pipelines:
  my_pipeline:
    # Required
    nf_core_name: my_pipeline          # nf-core pipeline name
    
    # Optional - Override defaults
    default_version: "1.0.0"           # Pipeline version (default: latest)
    analysis_category: MyCategory       # SUSHI category for grouping
    required_columns:                   # Required dataset columns
      - Name
      - Read1
      - Read2
    
    # Input handling
    input_type: samplesheet            # samplesheet, id_list, fasta_list, auto
    samplesheet_mapping:               # Column mapping for samplesheet
      sample: Name
      fastq_1: "Read1 [File]"
      fastq_2: "Read2 [File]"
    
    # Custom parameters (shown in UI)
    custom_params:
      - name: my_param
        type: text_area                # text_area, boolean, integer, select
        required: true
        default: ""
        description: "Description shown in UI"
    
    # Output options
    skip_multiqc: false                # Set true if pipeline doesn't generate MultiQC
    
    # Resource defaults
    params:
      cores: 8
      ram: 30
      scratch: 100
```

---

## Creating a Manual SUSHI App

For pipelines requiring special handling, create a static Ruby class:

### File Location

Create `lib/NfCore<PipelineName>App.rb`:

```ruby
#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreMyPipelineApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'NfCoreMyPipeline'
    @analysis_category = 'MyCategory'
    @description = <<-EOS
My pipeline description.
<a href='https://nf-co.re/mypipeline'>nf-core/mypipeline</a>
EOS
    @required_columns = ['Name', 'Read1', 'Read2']
    @required_params = ['my_required_param']
    
    # Parameters
    @params['process_mode'] = 'DATASET'
    @params['nfcorePipeline'] = 'mypipeline'
    @params['pipelineVersion'] = '1.0.0'
    @params['cores'] = 8
    @params['ram'] = 30
    @params['scratch'] = 100
    
    # Custom parameters
    @params['my_param'] = ''
    @params['my_param', 'description'] = 'Enter value here'
    
    @modules = ["Dev/jdk", "Tools/Nextflow"]
  end
  
  def next_dataset
    result_dir = File.join(@result_dir, "#{@params['name']}_result")
    {
      'Name' => @params['name'],
      'Result [File]' => result_dir,
      'MultiQC [Link]' => File.join(result_dir, 'multiqc', 'multiqc_report.html')
    }
  end
  
  def commands
    r_app_path = '/srv/GT/analysis/masaomi/2026/FGCZ/nf_core_sushi_app_20260109/EzAppNfCoreGeneric.R'
    
    cmd = run_RApp('EzAppNfCoreGeneric')
    
    # Add Apptainer cache settings
    cache_settings = <<~SHELL
      export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
      export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
    SHELL
    
    cmd = cache_settings + cmd
    cmd.sub!("}\nparam = list()", "}\nsource('#{r_app_path}')\nparam = list()")
    cmd
  end
end
```

### Priority

Static `.rb` files take precedence over dynamically generated classes.

---

## Input Types

### samplesheet (Default)

Standard CSV format for most nf-core pipelines:

```yaml
input_type: samplesheet
samplesheet_mapping:
  sample: Name
  fastq_1: "Read1 [File]"
  fastq_2: "Read2 [File]"
```

### id_list

For pipelines like `fetchngs` that accept accession IDs:

```yaml
input_type: id_list
custom_params:
  - name: accession_ids
    type: text_area
    required: true
    description: "Enter SRA/ENA/GEO IDs, comma-separated"
```

### fasta_list

For pipelines that accept FASTA file lists:

```yaml
input_type: fasta_list
custom_params:
  - name: fasta_files
    type: text_area
    required: true
    description: "Enter FASTA file paths"
```

---

## Parameter Configuration

### Parameter Types in YAML

| Type | UI Display | Example |
|------|------------|---------|
| `text_area` | Multi-line text field | Accession IDs, file lists |
| `boolean` | Checkbox | `skip_trimming` |
| `integer` | Number input | `max_cpus` |
| `select` | Dropdown | Predefined options |

### Example

```yaml
custom_params:
  - name: genome
    type: select
    options: ["GRCh38", "GRCm39", "BDGP6"]
    default: "GRCh38"
    description: "Reference genome"
  
  - name: skip_qc
    type: boolean
    default: false
    description: "Skip quality control steps"
```

### Required Parameters

Add `required: true` to enforce input validation:

```yaml
custom_params:
  - name: my_param
    required: true  # User must fill this before proceeding
```

---

## Samplesheet Mapping

Map SUSHI dataset columns to nf-core samplesheet columns:

```yaml
samplesheet_mapping:
  # nf-core column: SUSHI column
  sample: Name
  fastq_1: "Read1 [File]"
  fastq_2: "Read2 [File]"
  strandedness: StrandMode
  antibody: Antibody
  control: Control
```

The `[File]` suffix is handled automatically - both `Read1` and `Read1 [File]` will match.

---

## Nextflow & Apptainer Tips

### Apptainer/Singularity Cache

All nf-core container images are cached at:
```
/misc/fgcz01/nextflow_apptainer_cache/
```

This is set automatically via environment variables:
```bash
export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
```

### Resource Configuration

A custom Nextflow config is generated for each job:

```groovy
executor {
  name = 'local'
  cpus = 8
  memory = '64.GB'
}

process {
  cpus = 8
  memory = '30.GB'
}
```

### Version Pinning

Always specify a version in `default_version`:

```yaml
rnaseq:
  default_version: "3.14.0"  # Don't use "master" in production
```

### Checking Available Versions

```bash
# List releases
curl -s https://api.github.com/repos/nf-core/rnaseq/releases | jq '.[].tag_name'
```

---

## Troubleshooting

### Pipeline Not Appearing

1. Check if pipeline exists at `https://nf-co.re/pipelines`
2. Clear cache: `rm -rf master/tmp/nfcore_cache/`
3. Restart Rails server

### Parameters Not Loading

1. Verify `nextflow_schema.json` exists in the pipeline repo
2. Check Rails logs for API errors
3. Override in `nf_core_pipelines.yml`

### Input File Format Error

nf-core pipelines require specific file extensions:
- Samplesheet: `.csv`, `.tsv`
- ID list: `.csv` (not `.txt`)

### Container Pull Failures

```bash
# Check cache directory permissions
ls -la /misc/fgcz01/nextflow_apptainer_cache/

# Manually pull container
singularity pull docker://nfcore/rnaseq:3.14.0
```

### Job Fails Immediately

Check the job log files:
```
/srv/gstore/projects/pXXXX/oXXXXX_NfCore.../scripts/*_o.log
/srv/gstore/projects/pXXXX/oXXXXX_NfCore.../scripts/*_e.log
```

---

## Quick Reference

### Configuration File Locations

| File | Purpose |
|------|---------|
| `config/nf_core_pipelines.yml` | Pipeline customization |
| `lib/NfCore*App.rb` | Static app definitions |
| `tmp/nfcore_cache/` | API response cache (24h TTL) |

### Common Settings

```yaml
# Minimal configuration
rnaseq:
  default_version: "3.14.0"

# Full configuration
fetchngs:
  input_type: id_list
  skip_multiqc: true
  custom_params:
    - name: accession_ids
      type: text_area
      required: true
```

### Useful Commands

```bash
# Clear nf-core cache
rm -rf master/tmp/nfcore_cache/

# Check registered apps (Rails console)
NfCoreAppFactory.registered_class_names

# Force reload
NfCoreAppFactory.register_dynamic_apps
```
