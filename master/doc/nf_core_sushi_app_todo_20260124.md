# nf-core SUSHI App Integration - TODO

**Date**: 2026-01-24  
**Branch**: nf-core-sushi-app

## Current Status

- [x] Dynamic nf-core app loading from API (142 pipelines)
- [x] fetchngs pipeline test - SUCCESS
- [ ] rnaseq pipeline test - IN PROGRESS (resource configuration issues)
- [x] refBuild selector integration for FGCZ references
- [x] Strandedness value conversion (SUSHI → nf-core)
- [x] Developer documentation created

## Next Steps

### 1. NfCoreRnaseqApp Test Completion
- [ ] Debug current resource/CPU configuration issue
- [ ] Verify STAR index generation with correct max_cpus setting
- [ ] Complete end-to-end rnaseq pipeline test
- [ ] Validate output dataset structure

### 2. Converter Improvement for Automatic SUSHI App Conversion
- [ ] Analyze nf-core API JSON structure for automatic mapping
- [ ] Implement automatic `samplesheet_mapping` inference from `schema_input.json`
  - `sample` → `Name`
  - `fastq_1` → `Read1 [File]`
  - `fastq_2` → `Read2 [File]`
  - `strandedness` → `StrandMode`
- [ ] Implement automatic `use_ref_selector` detection from `nextflow_schema.json`
  - Detect if pipeline has `--genome`, `--fasta`, `--gtf` parameters
- [ ] Implement automatic value conversion rules
  - Strandedness: `both` → `unstranded`, `sense` → `forward`, etc.
- [ ] Reduce YAML configuration to exceptions only

### 3. Pre-Release Tasks
- [ ] Test additional pipelines (chipseq, atacseq, sarek)
- [ ] Add more pipeline-specific configurations as needed
- [ ] Remove debug output from production code
- [ ] Update documentation with final implementation details
- [ ] Code review and cleanup

## Known Issues

1. **Nextflow version compatibility**: nf-core/rnaseq 3.22.2 requires Nextflow 25.04+, using 3.14.0 instead
2. **Resource configuration**: STAR genome generation requires proper max_cpus limiting
3. **Parameter visibility**: Internal parameters (inputType, samplesheetMapping) visible in UI

## Files Modified

- `lib/nf_core_app_factory.rb` - Dynamic class generation
- `lib/nf_core_info_fetcher.rb` - API integration and caching
- `config/nf_core_pipelines.yml` - Pipeline configurations
- `app/controllers/run_application_controller.rb` - Parameter loading
- `EzAppNfCoreGeneric.R` - Nextflow execution and input generation
