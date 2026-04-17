#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class NfCoreDenovotranscriptApp < SushiFabric::SushiApp
  include GlobalVariables

  # Override run to call SushiFabric::SushiApp's run method directly
  # (skipping GlobalVariables' run which requires an argument)
  def run
    SushiFabric::SushiApp.instance_method(:run).bind(self).call
  end

  def initialize
    super
    @name = 'NfCoreDenovotranscript'
    @params['process_mode'] = 'DATASET'
    @analysis_category = 'Transcriptomics'
    @description = <<-EOS
      De novo transcriptome assembly pipeline from nf-core.<br/>
      Assembles transcriptomes from bulk RNA-seq short reads without a reference genome.<br/>
      Supports Trinity and rnaSPAdes assemblers with EvidentialGene integration.<br/>
      <a href='https://nf-co.re/denovotranscript/dev/'>nf-core/denovotranscript</a>
    EOS
    @required_columns = ['Name', 'Read1']
    @required_params = ['pipelineVersion']

    # nf-core pipeline identifiers
    @params['nfcorePipeline'] = 'denovotranscript'
    @params['pipelineVersion'] = 'dev'
    @params['pipelineVersion', 'description'] = 'Pipeline version (use "dev" for development branch)'

    # Input configuration for EzAppNfCoreGeneric
    @params['inputType'] = 'samplesheet'
    @params['samplesheetMapping'] = {
      'sample' => 'Name',
      'fastq_1' => 'Read1 [File]',
      'fastq_2' => 'Read2 [File]'
    }.to_json

    # Slurm resource parameters
    @params['cores'] = '8'
    @params['cores', 'context'] = 'slurm'
    @params['ram'] = '100'
    @params['ram', 'context'] = 'slurm'
    @params['scratch'] = '500'
    @params['scratch', 'context'] = 'slurm'

    # Nextflow module version
    @params['nextflowModuleVersion'] = ['25.10', '25.04', '24.10', '24.04']
    @params['nextflowModuleVersion', 'description'] = 'Nextflow version to load via module'

    # Assembly options
    @params['assemblers'] = ['trinity,rnaspades', 'trinity', 'rnaspades', 'trinity_no_norm', 'trinity,rnaspades,trinity_no_norm']
    @params['assemblers', 'description'] = 'Comma-separated assemblers to use: trinity, rnaspades, trinity_no_norm'

    # Strand-specificity for rnaSPAdes
    @params['ss'] = ['', 'rf', 'fr']
    @params['ss', 'description'] = 'Strand-specificity for rnaSPAdes (leave empty for unstranded)'

    # rRNA removal
    @params['remove_ribo_rna'] = false
    @params['remove_ribo_rna', 'description'] = 'Enable ribosomal RNA removal with SortMeRNA'

    # BUSCO options
    @params['busco_lineage'] = 'auto'
    @params['busco_lineage', 'description'] = 'BUSCO lineage for assembly QC (e.g., mammalia_odb10, auto)'

    # Workflow mode options
    @params['qc_only'] = false
    @params['qc_only', 'description'] = 'Only perform QC, skip assembly and quantification'

    # Additional command options
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Additional Nextflow command-line options'

    @params['mail'] = ''
    @params['name'] = 'NfCoreDenovotranscript'

    @modules = ["Dev/jdk", "Tools/Nextflow"]
  end

  def next_dataset
    result_dir = File.join(@result_dir, "#{@params['name']}_result")

    dataset = {
      'Name' => @params['name'],
      'Result [File]' => result_dir,
      'MultiQC [Link]' => File.join(result_dir, 'multiqc', 'multiqc_report.html')
    }

    if @dataset && @dataset.first
      inherit_cols = @dataset.first.keys.reject do |col|
        col == 'Name' || col =~ /^Read[12]/ || col == 'Species'
      end
      inherit_cols.each do |col|
        dataset[col] = @dataset.first[col]
      end
    end

    dataset
  end

  def commands
    cmd = run_RApp('EzAppNfCoreGeneric')

    # Apptainer cache settings
    cache_settings = <<~SHELL
      export NXF_SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
      export SINGULARITY_CACHEDIR=/misc/fgcz01/nextflow_apptainer_cache/
    SHELL

    cmd = cache_settings + cmd

    # Source external R file (EzAppNfCoreGeneric is not yet in ezRun)
    r_app_source = '/srv/GT/analysis/masaomi/2026/FGCZ/nf_core_sushi_app_20260109/EzAppNfCoreGeneric.R'
    cmd.sub!("}\nparam = list()", "}\nsource('#{r_app_source}')\nparam = list()")

    cmd
  end
end

if __FILE__ == $0
  usecase = NfCoreDenovotranscriptApp.new
  usecase.project = "p1001"
  usecase.user = "masaomi"
end
