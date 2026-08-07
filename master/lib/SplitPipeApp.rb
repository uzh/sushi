#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class SplitPipeApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'ParseSplitPipe'
    @analysis_category = 'SingleCell'
    @description =<<-EOS
This wrapper runs the <a href='https://support.parsebiosciences.com/'>Parse Biosciences split-pipe</a> pipeline (Evercode WT split-seq).<br/>
Every selected row is treated as one sublibrary of the same experiment: each is processed with <tt>split-pipe --mode all</tt> and, when more than one sublibrary is selected, they are merged with <tt>split-pipe --mode comb</tt> into a single per-sample result.<br/>
Read1 must be the mRNA/cDNA read (--fq1) and Read2 the barcode read (--fq2).
    EOS
    ## general
    @required_columns = ['Name', 'Read1', 'Read2', 'Species']
    @required_params = ['name', 'refBuild', 'kit', 'chemistry']
    @params['process_mode'] = 'DATASET'
    @params['cores'] = ['16', '8', '32']
    @params['cores', "context"] = "slurm"
    @params['ram'] = ['60', '40', '120']
    @params['ram', "context"] = "slurm"
    @params['scratch'] = '300'
    @params['scratch', "context"] = "slurm"
    @params['name'] = 'ParseSplitPipe'
    ## reference
    @params['refBuild'] = ref_selector
    @params['refBuild', "context"] = "referfence genome assembly"
    @params['refFeatureFile'] = 'genes.gtf'
    @params['refFeatureFile', "context"] = "SplitPipe"
    @params['featureLevel'] = 'gene'
    @params['featureLevel', "context"] = "SplitPipe"
    @params['transcriptTypes'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'multi_selection'] = true
    @params['transcriptTypes', 'selected'] = ['protein_coding', 'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA', 'long_noncoding', 'short_noncoding', 'pseudogene']
    @params['transcriptTypes', 'description'] = 'Transcript types kept when building the reference on the fly. The default (all types) matches the pre-built shared Parse indices, so an existing index is reused; a different selection triggers a fresh split-pipe --mode mkref build.'
    ## split-pipe
    @params['kit'] = ['WT', 'WT_mini', 'WT_mega', 'WT_mega_384', 'WT_penta', 'WT_penta_384']
    @params['kit', 'description'] = 'Parse Evercode WT kit. WT_mega_384/WT_penta/WT_penta_384 require chemistry v3 or v4.'
    @params['chemistry'] = ['v3', 'v2', 'v1', 'v4']
    @params['chemistry', 'description'] = 'Parse chemistry version.'
    @params['sampleWells'] = 'all-well A1-A12'
    @params['sampleWells', 'description'] = "Sample-to-well mapping passed as --sample. One or more '<name> <wells>' specs separated by ';', e.g. 'all-well A1-A12' or 'sampleA A1-A6; sampleB A7-A12'. The well range must match the kit (WT_mini=12/A1-A12, WT=48/A1-D12, WT_mega=96/A1-H12). Ignored when a sampleLoadingTable is given."
    @params['sampleLoadingTable'] = ''
    @params['sampleLoadingTable', 'description'] = 'Optional full path to a Parse SampleLoadingTable (.xlsx) or a plain "<name> <wells>" per-line file. When set, it overrides sampleWells (--samp_sltab / --samp_list).'
    @params['saveAnndata'] = false
    @params['saveAnndata', 'description'] = 'Also write an AnnData (.h5ad) output (--save_anndata).'
    @params['cmdOptions'] = ''
    @params['cmdOptions', 'description'] = 'Extra split-pipe --mode all flags (e.g. --tcr_analysis, --bcr_analysis, --crispr). Do not repeat options already set by the app.'
    @params['cmdOptions', "context"] = "SplitPipe"
    @params['specialOptions'] = ''
    @params['mail'] = ""
    @modules = ["Dev/R"]
    @inherit_tags = ["Factor", "B-Fabric"]
  end
  def set_default_parameters
  end
  def next_dataset
    # Experiment-level row: report + result dir. The per-biological-sample count
    # matrices are emitted as grandchild datasets (one row per --sample), so
    # downstream ScSeurat / CellBender can run per sample.
    report_dir = File.join(@result_dir, @params['name'])
    {
      'Name' => @params['name'],
      'Species' => (dataset = @dataset.first and dataset['Species']),
      'refBuild' => @params['refBuild'],
      'refFeatureFile' => @params['refFeatureFile'],
      'featureLevel' => @params['featureLevel'],
      'transcriptTypes' => @params['transcriptTypes'],
      'SCDataOrigin' => 'ParseBio',
      'ResultDir [File]' => report_dir,
      'Report [Link]' => File.join(report_dir, '00index.html')
    }.merge(extract_columns(@inherit_tags))
  end
  def grandchild_datasets
    # One row per biological sample, pointing at the 10x-format matrices the R
    # method writes (filtered_feature_bc_matrix / raw_feature_bc_matrix). Sample
    # names are derived from the sampleWells specs (';'-separated '<name> <wells>').
    # When a sampleLoadingTable is used or sampleWells is empty, the sample names
    # are not known here, so no grandchildren are emitted (the experiment row and
    # the on-disk 10x matrices are still produced).
    return [] if @params['sampleLoadingTable'].to_s.strip != ''
    return [] if @params['sampleWells'].to_s.strip == ''

    sample_names = @params['sampleWells'].to_s.split(';').map { |spec|
      spec.strip.split(/\s+/).first
    }.compact.reject(&:empty?).uniq

    report_dir = File.join(@result_dir, @params['name'])
    species = (dataset = @dataset.first and dataset['Species'])
    sample_names.map do |sample_name|
      sample_dir = File.join(report_dir, sample_name)
      {
        'Name' => sample_name,
        'Species' => species,
        'refBuild' => @params['refBuild'],
        'refFeatureFile' => @params['refFeatureFile'],
        'featureLevel' => @params['featureLevel'],
        'transcriptTypes' => @params['transcriptTypes'],
        'SCDataOrigin' => 'ParseBio',
        'ResultDir [File]' => sample_dir,
        'Report [Link]' => File.join(report_dir, '00index.html'),
        'CountMatrix [Link]' => File.join(sample_dir, 'filtered_feature_bc_matrix'),
        'UnfilteredCountMatrix [Link]' => File.join(sample_dir, 'raw_feature_bc_matrix')
      }.merge(extract_columns(@inherit_tags))
    end
  end
  def commands
    run_RApp("EzAppSplitPipe", conda_env: "gi_parse_v1.8.2")
  end
end

if __FILE__ == $0

end
