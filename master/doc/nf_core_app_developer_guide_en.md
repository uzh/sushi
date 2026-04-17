# nf-core SUSHI App Developer Guide

Last updated: 2026-02-06

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Dynamic Class Generation (Metaprogramming)](#2-dynamic-class-generation-metaprogramming)
3. [Convention over Configuration](#3-convention-over-configuration)
4. [NfCoreAppFactory in Detail](#4-nfcoreappfactory-in-detail)
5. [NfCoreInfoFetcher in Detail](#5-nfcoreinfofetcher-in-detail)
6. [NfCoreConventions in Detail](#6-nfcoreconventions-in-detail)
7. [Nextflow Version Auto-Compatibility](#7-nextflow-version-auto-compatibility)
8. [R Script (EzAppNfCoreGeneric.R)](#8-r-script-ezappnfcoregenericr)
9. [End-to-End Execution Flow](#9-end-to-end-execution-flow)
10. [Adding New Pipelines](#10-adding-new-pipelines)
11. [Debugging & Troubleshooting](#11-debugging--troubleshooting)

---

## 1. Architecture Overview

The nf-core SUSHI App system enables automatic execution of nf-core pipelines within the SUSHI framework.

### Component Architecture

```
┌─────────────────────────────────────────────┐
│              SUSHI Web UI                    │
│   (set_parameters.html.erb)                 │
│   (run_application_controller.rb)           │
├─────────────────────────────────────────────┤
│          NfCoreAppFactory                   │
│   ┌───────────────────────┐                 │
│   │ Dynamic Class Creation│ (Class.new)     │
│   │ (Metaprogramming)     │                 │
│   └───────────┬───────────┘                 │
│               ↓                             │
│   ┌─────────────────────┐                   │
│   │ NfCoreInfoFetcher   │ ← nf-core API    │
│   │ (API + Caching)     │ ← GitHub Raw     │
│   └─────────────────────┘                   │
│   ┌─────────────────────┐                   │
│   │ NfCoreConventions   │ (Static defs)     │
│   │ (Conventions)       │                   │
│   └─────────────────────┘                   │
│   ┌─────────────────────┐                   │
│   │ nf_core_pipelines   │ (Exceptions only) │
│   │ .yml                │                   │
│   └─────────────────────┘                   │
├─────────────────────────────────────────────┤
│         SushiFabric::SushiApp               │
│         (Base class)                        │
├─────────────────────────────────────────────┤
│    EzAppNfCoreGeneric.R                     │
│    (Samplesheet generation / NF execution)  │
├─────────────────────────────────────────────┤
│         Nextflow + Singularity              │
│         (Pipeline execution)                │
└─────────────────────────────────────────────┘
```

### File Overview

| File | Role |
|------|------|
| `lib/nf_core_app_factory.rb` | Main module for dynamic class generation |
| `lib/nf_core_info_fetcher.rb` | nf-core API communication and caching |
| `lib/nf_core_conventions.rb` | Static definitions for Convention over Configuration |
| `config/nf_core_pipelines.yml` | Exception pipeline configuration (minimal) |
| `app/controllers/run_application_controller.rb` | Dynamic UI parameter injection |
| `app/views/run_application/set_parameters.html.erb` | Parameter setting view template |
| `lib/sushi_fabric/lib/sushi_fabric/sushiApp.rb` | SUSHI App base class |
| `EzAppNfCoreGeneric.R` | R script (samplesheet generation / execution) |

---

## 2. Dynamic Class Generation (Metaprogramming)

### Overview

With 100+ nf-core pipelines, it is impractical to maintain individual Ruby class files for each one. Instead, Ruby's metaprogramming capabilities are used to dynamically generate classes at SUSHI startup.

### Class.new for Dynamic Class Creation

Ruby's `Class.new` method creates a new class inheriting from `SushiFabric::SushiApp`:

```ruby
# Conceptual overview of NfCoreAppFactory#create_dynamic_class
klass = Class.new(SushiFabric::SushiApp) do
  # Closure captures config variables
  define_method(:initialize) do
    super()
    @name = app_config['name']
    @analysis_category = app_config['analysis_category']
    @description = app_config['description']
    @params['nfcorePipeline'] = app_config['nf_core_name']
    @params['pipelineVersion'] = app_config['default_version']
    # ... parameter setup
  end

  define_method(:next_dataset) do
    # Define output dataset
  end

  define_method(:commands) do
    # Generate R script execution command
    cmd = run_RApp('EzAppNfCoreGeneric')
    cmd.sub!("}\nparam = list()", "}\nsource('#{r_app_path}')\nparam = list()")
    cmd
  end
end

# Register as global constant
Object.const_set("NfCoreRnaseqApp", klass)
```

### Closure Usage

Methods defined with `define_method` capture outer-scope variables (`app_config`, `r_app_path`) as closures. This ensures each pipeline's configuration is correctly bound to its class instances.

### Object.const_set for Constant Registration

```ruby
Object.const_set(class_name, klass)
```

This registers the dynamically created class as a top-level constant like `NfCoreRnaseqApp`, making it visible in SUSHI's App list.

### Metaprogrammatic Parameter Processing

Parameters fetched from the API are automatically mapped to SUSHI UI components based on their `type`:

| API Parameter Type | SUSHI Display | Implementation |
|-------------------|---------------|----------------|
| `enum` | Dropdown | `@params[name] = enum_values` (Array) |
| `boolean` | true/false dropdown | `@params[name] = true/false` |
| `integer` | Text field (numeric) | `@params[name] = value.to_i` |
| `string` | Text field | `@params[name] = value.to_s` |

---

## 3. Convention over Configuration

### Design Philosophy

This system follows the same design philosophy as Rails' "Convention over Configuration":

- **Follow conventions, no config needed**: Standard nf-core pipelines work without YAML
- **Configure exceptions only**: Only special pipelines need `nf_core_pipelines.yml` entries
- **Auto-detect from schemas**: Extract as much information as possible from API schemas

### Configuration Resolution Priority

All configuration items are resolved in this order:

```
YAML Explicit Config > API Auto-Detection > Convention Default > Global Default
```

Specific examples:

| Config Item | YAML | Auto-Detection | Convention | Default |
|-------------|------|---------------|-----------|---------|
| category | `analysis_category:` | Inferred from nf-core topic | `CATEGORY_MAPPING` | `'nf-core'` |
| resources | `params: cores/ram/scratch` | - | `RESOURCE_DEFAULTS[category]` | `cores:8, ram:30, scratch:100` |
| input_type | `input_type:` | Parsed from `schema_input.json` | - | `'samplesheet'` |
| samplesheet_mapping | `samplesheet_mapping:` | `schema_input.json` + `COLUMN_MAPPING` | - | `{}` |
| use_ref_selector | `use_ref_selector:` | Parsed from `nextflow_schema.json` | - | `false` |
| version | `default_version:` | Nextflow-compatible version | - | `latest_version` |

### When YAML Configuration is Needed

YAML is only required for:

1. **Custom UI components** - Dropdowns, text areas, or other special inputs
2. **Version pinning** - When you want a specific version instead of auto-selection
3. **Resource overrides** - When the pipeline needs resources different from the category standard
4. **Special input formats** - `id_list` (fetchngs) or `fasta_list` that aren't standard
5. **MultiQC skipping** - For pipelines that don't generate MultiQC reports

---

## 4. NfCoreAppFactory in Detail

### Module Structure

```ruby
module NfCoreAppFactory
  PIPELINE_CONFIG_PATH  # YAML config file path
  LIB_DIR               # Generated file output directory
  DEFAULTS              # Global default values

  # Key methods
  def self.load_pipelines          # Load pipeline list
  def self.build_config(key, info) # Build config (apply conventions)
  def self.register_dynamic_apps   # Register dynamic classes
  def self.create_dynamic_class    # Create class via Class.new
end
```

### load_pipelines Flow

```
1. NfCoreInfoFetcher.fetch_all_pipelines_with_info
   → Get all pipeline info from pipelines.json

2. YAML.load_file(nf_core_pipelines.yml)
   → Load manual configuration

3. deep_merge(api_info, yaml_config)
   → YAML overrides API data

4. Apply defaults
   → Apply default values to unset items
```

### build_config Flow

```
Input: pipeline_info (API + YAML merged)
  ↓
1. Category:  YAML > NfCoreConventions.category_for > 'nf-core'
  ↓
2. Resources: Convention(category) + YAML merge
  ↓
3. InputType: YAML > NfCoreInfoFetcher.detect_input_type > 'samplesheet'
  ↓
4. Mapping:   YAML > NfCoreInfoFetcher.auto_samplesheet_mapping > {}
  ↓
5. RefSelector: YAML > NfCoreInfoFetcher.needs_ref_selector? > false
  ↓
6. Version:   YAML > find_compatible_version > latest > 'master'
  ↓
Output: config Hash (all settings finalized)
```

### register_dynamic_apps Flow

```
1. load_pipelines() → all_pipelines
  ↓
2. for each pipeline:
   a. build_config → config
   b. class_name = "NfCore#{camelize(name)}App"
   c. skip if Object.const_defined?(class_name)  ← Don't overwrite existing classes
   d. create_dynamic_class(config) → klass
   e. Object.const_set(class_name, klass)
  ↓
3. Add to @@registered_classes
```

---

## 5. NfCoreInfoFetcher in Detail

### API Endpoints

| URL Template | Data Retrieved | Cache TTL |
|-------------|---------------|-----------|
| `https://nf-co.re/pipelines.json` | All pipeline listings | 24 hours |
| `https://raw.githubusercontent.com/nf-core/{name}/master/assets/schema_input.json` | Samplesheet column definitions | 24 hours |
| `https://raw.githubusercontent.com/nf-core/{name}/{version}/nextflow_schema.json` | Parameter schema | 24 hours |
| `https://raw.githubusercontent.com/nf-core/{name}/{version}/nextflow.config` | Nextflow requirements | 24 hours |

### Caching Mechanism

```ruby
def self.fetch_with_cache(url, cache_path)
  # 1. Cache is valid (< 24 hours) → Use cache
  # 2. Connect to API → Success: Update cache + return data
  # 3. API failure → Use stale cache if available
  # 4. No cache → Return nil
end
```

Cache files are stored in `tmp/nfcore_cache/`:
- `pipelines.json` - All pipeline info
- `{name}_schema_input.json` - Samplesheet schema
- `{name}_{version}_nextflow_schema.json` - Parameter schema
- `{name}_{version}_nextflow.config` - Nextflow config

### Auto-Detection Methods

| Method | Input | Output | Description |
|--------|-------|--------|-------------|
| `detect_input_type` | pipeline_name | `'samplesheet'/'id_list'/'fasta_list'` | Infer from schema_input.json properties |
| `auto_samplesheet_mapping` | pipeline_name | `{nf_col => sushi_col}` | Auto-generate column mapping |
| `needs_ref_selector?` | pipeline_name, version | `true/false` | Detect reference genome section |
| `fetch_all_params` | pipeline_name, version | `{required: [...], optional: [...]}` | Get all params categorized |

---

## 6. NfCoreConventions in Detail

### COLUMN_MAPPING

Maps nf-core samplesheet column names to SUSHI dataset column names. Special cases can include default values:

```ruby
# Standard mapping
'sample'       => 'Name'
'fastq_1'      => 'Read1 [File]'
'fastq_2'      => 'Read2 [File]'

# Mapping with default value
'strandedness' => { sushi_col: 'StrandMode', default: 'auto' }
```

### CATEGORY_MAPPING

Infers SUSHI categories from nf-core topic tags:

```ruby
'rna-seq'     => 'Transcriptomics'
'chip-seq'    => 'Epigenetics'
'dna-seq'     => 'Genomics'
'single-cell' => 'SingleCell'
```

### RESOURCE_DEFAULTS

Provides resource defaults based on category:

```ruby
'Transcriptomics' => { 'cores' => 8,  'ram' => 64,  'scratch' => 200 }
'SingleCell'       => { 'cores' => 16, 'ram' => 128, 'scratch' => 500 }
'QC'               => { 'cores' => 4,  'ram' => 16,  'scratch' => 100 }
```

---

## 7. Nextflow Version Auto-Compatibility

### Background

Newer nf-core pipeline releases often require newer Nextflow versions. When the server's installed Nextflow is older, the latest pipeline version won't work.

### How Auto-Compatibility Works

```
1. installed_nextflow_version
   Execute `nextflow -v` to get version string like "24.10.3"
   (Cached per process lifetime)

2. find_compatible_version(pipeline_name, releases)
   Iterate releases from newest to oldest:
     a. fetch_required_nextflow_version(pipeline, version)
        → Parse "nextflowVersion = '!>=24.10.5'" from nextflow.config
     b. nextflow_version_compatible?(required)
        → installed >= required? → Yes: Return this version
     c. No: Try next older release

3. Version resolution priority:
   YAML default_version > compatible latest > API latest > 'master'
```

### Version Comparison

```ruby
def self.compare_versions(v1, v2)
  # "24.10.3" vs "24.10.5" → [-1, 0, 1]
  # Semantic versioning comparison
  parts1 = v1.split('.').map(&:to_i)
  parts2 = v2.split('.').map(&:to_i)
  parts1 <=> parts2
end
```

### Automatic Updates

When Nextflow is updated on the server:
- Next SUSHI startup detects new version via `installed_nextflow_version`
- `find_compatible_version` selects newer pipeline versions
- **No manual changes required**

---

## 8. R Script (EzAppNfCoreGeneric.R)

### Overview

`EzAppNfCoreGeneric.R` is a generic R script that receives SUSHI parameters and executes nf-core pipelines.

### Main Processing Flow

```
1. Receive Parameters
   - nfcorePipeline, pipelineVersion, inputType
   - samplesheetMapping (JSON string)
   - refBuild (reference genome)

2. Build Samplesheet (buildSamplesheetFromMapping)
   - Convert SUSHI dataset → nf-core samplesheet format
   - Map columns based on samplesheetMapping
   - Handle strandedness defaults (RNA-seq: 'auto')

3. Resolve References
   - Resolve FASTA, GTF paths from refBuild
   - Pass as --fasta, --gtf options to Nextflow

4. Execute Nextflow
   nextflow run nf-core/{pipeline} -r {version} \
     --input samplesheet.csv \
     --outdir _result \
     -profile singularity \
     -c custom_resources.config \
     [--fasta ... --gtf ...]

5. Process Results
   - Organize output files
   - Generate MultiQC report links
```

### Strandedness Default Logic

```r
# Only for RNA-seq related pipelines
if (pipeline %in% c("rnaseq", "smrnaseq", "rnafusion")) {
  if (!"StrandMode" %in% colnames(dataset)) {
    samplesheet$strandedness <- "auto"
  }
}
```

---

## 9. End-to-End Execution Flow

### At SUSHI Startup

```
SUSHI Server Start
  ↓
run_application_controller.rb: load nf-core apps
  ↓
NfCoreAppFactory.register_dynamic_apps
  ↓
  ├── NfCoreInfoFetcher.fetch_all_pipelines_with_info
  │   └── fetch pipelines.json (cache or API)
  │       → 100+ pipeline info records
  │
  ├── YAML.load_file(nf_core_pipelines.yml)
  │   → Exception config (rnaseq, fetchngs, smrnaseq)
  │
  ├── for each pipeline:
  │   ├── build_config (apply conventions)
  │   │   ├── category_for(topic)
  │   │   ├── resources_for(category)
  │   │   ├── detect_input_type(name)
  │   │   ├── auto_samplesheet_mapping(name)
  │   │   ├── needs_ref_selector?(name, version)
  │   │   └── find_compatible_version(name, releases)
  │   │
  │   ├── create_dynamic_class(config)
  │   │   ├── fetch_all_params(name, version)
  │   │   ├── Class.new(SushiFabric::SushiApp) { ... }
  │   │   └── define_method(:initialize, :next_dataset, :commands)
  │   │
  │   └── Object.const_set("NfCore#{Name}App", klass)
  │
  └── Result: 100+ SUSHI App classes registered
```

### At Job Execution

```
User: Select NfCoreRnaseqApp
  ↓
run_application_controller#set_parameters
  ├── Build parameter form
  ├── Inject refBuild selector (if use_ref_selector == true)
  └── Display to user
  ↓
User: Configure parameters + Submit
  ↓
SushiFabric::SushiApp#run
  ↓
NfCoreRnaseqApp#commands
  ├── run_RApp('EzAppNfCoreGeneric')
  ├── source('EzAppNfCoreGeneric.R')
  └── Generate shell script
  ↓
Cluster Job Execution (shell script)
  ↓
EzAppNfCoreGeneric.R
  ├── Build samplesheet from dataset
  ├── Resolve reference genome
  ├── Build Nextflow command
  └── Execute: nextflow run nf-core/rnaseq -r 3.19.0 ...
  ↓
Nextflow Pipeline Execution
  ├── Pull containers (Singularity)
  ├── Run pipeline steps
  └── Generate output + MultiQC report
  ↓
SUSHI: Register output dataset
```

---

## 10. Adding New Pipelines

### Three Approaches

There are three approaches to making an nf-core pipeline available as a SUSHI App:

| Approach | How | Use Case |
|----------|-----|----------|
| **Automatic (Convention)** | Do nothing | Standard pipelines |
| **YAML Exception Config** | Add to `nf_core_pipelines.yml` | Custom params / version pinning |
| **Static .rb File** | Place `lib/NfCore*App.rb` | Custom R script / complex logic |

### Static .rb File Definition

For pipelines that require logic beyond what dynamic generation supports, you can place a static Ruby class file in `lib/` for full customization.

**How it works:**

`NfCoreAppFactory.register_dynamic_apps` checks before generating each pipeline class:

```ruby
# Skip creating if class is already defined (static .rb file or previous registration)
if Object.const_defined?(class_name)
  next
end
```

Since SUSHI loads `.rb` files from `lib/` before calling `register_dynamic_apps`, any class defined by a static file takes precedence. **Static files always win** over dynamic generation.

**Existing static files:**

| File | Pipeline | Features |
|------|----------|----------|
| `lib/NfCoreCutNRunApp.rb` | nf-core/cutandrun | Custom R script (`EzAppNfCoreCutAndRun`), spike-in genome config, IGV session |
| `lib/NfCoreSmRnaSeqApp.rb` | nf-core/smrnaseq | Custom R script (`EzAppNfCoreSmRnaSeq`), `grandchild_datasets` for per-sample count files |
| `lib/NfCoreAtacSeqApp.rb` | nf-core/atacseq | Custom R script (`EzAppNfCoreAtacSeq`), `grandchild_datasets`, two-group analysis option |

**When static files are appropriate:**

1. **Custom R scripts** - Pipeline-specific pre/post-processing not handled by `EzAppNfCoreGeneric.R`
2. **grandchild_datasets** - Need to produce per-sample output datasets
3. **set_default_parameters** - Dynamic defaults that depend on other params (e.g., paired → add Read2 to required)
4. **Complex next_dataset** - Non-standard output directory structure
5. **Special modules** - Additional tools like `Tools/BEDTools` required

**Static file template:**

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
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Genomics'
    @description = <<-EOS
      My custom pipeline description.
      <a href='https://nf-co.re/mypipeline'>nf-core/mypipeline</a>
    EOS
    @required_columns = ['Name', 'Read1', 'Read2']
    @required_params = ['refBuild', 'pipelineVersion']
    
    # Resource parameters
    @params['cores'] = '8'
    @params['ram'] = '100'
    @params['scratch'] = '200'
    
    # Reference selector
    @params['refBuild'] = ref_selector
    
    # Pipeline version
    @params['pipelineVersion'] = '1.0.0'
    
    # Custom parameters
    @params['myOption'] = ['optionA', 'optionB']
    @params['myOption', 'description'] = 'Select an option'
    
    @modules = ["Dev/jdk", "Tools/Nextflow"]
  end
  
  def next_dataset
    report_file = File.join(@result_dir, "#{@params['name']}_result")
    {
      'Name' => @params['name'],
      'Result [File]' => report_file,
      'MultiQC [Link]' => File.join(report_file, 'multiqc', 'multiqc_report.html')
    }
  end
  
  def commands
    run_RApp('EzAppMyCustomPipeline')  # Use a custom R script
  end
end
```

> **Note**: Static files do not benefit from `NfCoreAppFactory` automatic features (version auto-compatibility, API parameter fetching, convention application). All configuration must be managed manually.

### Standard Pipelines (Automatic)

**No action required.** When a new pipeline is added to nf-core:

1. After `pipelines.json` cache refresh (24 hours), it is auto-detected
2. `NfCoreAppFactory` automatically generates a class
3. Conventions are applied for configuration

### Special Pipelines (YAML Needed)

Add to `config/nf_core_pipelines.yml`:

```yaml
pipelines:
  my_special_pipeline:
    # Pin version
    default_version: "1.0.0"

    # Custom UI parameters
    custom_params:
      - name: special_option
        type: select
        required: true
        options: ["A", "B", "C"]
        default: "A"
        description: "Special option for this pipeline"

    # Resource overrides
    params:
      cores: 32
      ram: 256
      scratch: 1000
```

### Adding New Column Mappings

Add to `COLUMN_MAPPING` in `lib/nf_core_conventions.rb`:

```ruby
COLUMN_MAPPING = {
  # ... existing mappings ...
  'new_nfcore_column' => 'NewSushiColumn',
  # With default value:
  'another_column' => { sushi_col: 'SushiCol', default: 'default_val' }
}
```

### Adding New Categories

Add to `lib/nf_core_conventions.rb`:

```ruby
CATEGORY_MAPPING['new-topic'] = 'NewCategory'
RESOURCE_DEFAULTS['NewCategory'] = { 'cores' => 8, 'ram' => 64, 'scratch' => 200 }
```

---

## 11. Debugging & Troubleshooting

### Clear All Caches

```bash
rm -rf master/tmp/nfcore_cache/
```

Deletes all caches; data will be re-fetched from the API on next startup.

### Clear Cache for a Specific Pipeline

```bash
rm master/tmp/nfcore_cache/rnaseq_*
```

### Check Logs

The SUSHI server log contains messages like:

```
NfCoreInfoFetcher: Detected installed Nextflow version: 24.10.3
NfCoreInfoFetcher: Pipeline rnaseq v3.19.0 is compatible with Nextflow 24.10.3
NfCoreInfoFetcher: Pipeline genomeassembler latest (1.1.0) requires Nextflow >= 24.10.5,
  but installed is 24.10.3. Falling back to compatible version: 1.0.1
NfCoreAppFactory: Total nf-core apps tracked: 105, newly created: 103
```

### Inspect Dynamic Classes

```ruby
# From Rails console
NfCoreAppFactory.registered_class_names
# => ["NfCoreRnaseqApp", "NfCoreChipseqApp", ...]

NfCoreRnaseqApp.new.params
# => Parameter listing
```

### Directly Check API Data

```bash
# Pipeline count
curl -s https://nf-co.re/pipelines.json | jq '.remote_workflows | length'

# Samplesheet schema
curl -s https://raw.githubusercontent.com/nf-core/rnaseq/master/assets/schema_input.json | jq .

# Parameter schema
curl -s https://raw.githubusercontent.com/nf-core/rnaseq/3.19.0/nextflow_schema.json | jq '.definitions | keys'

# Nextflow requirement
curl -s https://raw.githubusercontent.com/nf-core/rnaseq/3.22.2/nextflow.config | grep nextflowVersion
```

### Common Issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| Pipeline not showing | API connection error | Check cache / network |
| Version error | Nextflow incompatible | Verify auto-compatibility is working |
| Parameters not displayed | nextflow_schema.json fetch failed | Clear cache |
| Samplesheet error | Column mapping mismatch | Check `COLUMN_MAPPING` / override in YAML |
| UI shows text box | Project Defaults loading bug | Clear browser cache (fixed) |
