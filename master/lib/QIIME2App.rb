#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class QIIME2App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'QIIME2App'
    @analysis_category = 'Metagenomics'
    @description =<<-EOS
    Data processing with QIIME2. For short reads/Illumina data only.
    <a href='https://qiime2.org/'>QIIME2 main website with tutorials and manuals.</a>

EOS
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name', 'Read1']
    @required_params = ['paired','trim_left','truncate_len','sampling_depth','max_rarefaction_depth','min_freq','min_samples','group']
    @params['cores'] = '2'
    @params['ram'] = '7'
    @params['scratch'] = '10'
    @params['paired'] = false
    @params['paired', 'description'] = 'whether the reads are paired end; if false then only Read1 is considered even if Read2 is available.'
    @params['trim_left'] = '0'
    @params['trim_left', 'description'] = 'Position at which sequences should be trimmed due to low quality. Default:0, assuming good quality reads.'
    @params['truncate_len'] = '150'
    @params['truncate_len', 'description'] = 'Position at which sequences should be truncated due to decrease in quality. Assuming good quality reads and 150bp illumina sequencing.'
    @params['sampling_depth'] = '1000'
    @params['sampling_depth', 'description'] = 'Total frequency that each sample should be rarefied to prior to computing diversity metrics.'
    @params['max_rarefaction_depth'] = '4000'
    @params['max_rarefaction_depth', 'description'] = 'Choose maximum depth for generating alpha rarefaction curves for exploring explore alpha diversity as a function of sampling depth.'
    @params['min_freq'] = '100'
    @params['min_freq', 'description'] = 'The minimum total frequency that a feature must have to be retained for differential abundance calculation.'
    @params['min_samples'] = '15'
    @params['min_samples', 'description'] = 'The minimum number of samples that a feature must be observed in to be retained for differential abundance calculation.' 
    @params['group'] = true
    @params['group', 'description'] = 'There needs to be a group assignment column. Ensure the column name is in the format "NAME [Factor]"'
    @params['grouping'] = ''
    @params['grouping', 'description'] = 'Type in the same of the group assignment column. If the name is in the format "NAME [Factor]" type in the NAME'
    @params['database'] = ["silva", "greengenes"]
    @params['database', 'description'] = 'Choose marker gene reference database'
    @params['primer'] = ["V1-V3(1)", "V1-V3(2)", "V3", "V4", "V3-V4", "V4-V5", "V3-V5", "V4-V6"]
    @params['primer', 'description'] = 'Region of 16S rRNA genes to be used for training Naive Bayes classifier. It has been shown taxonomic classification improves when Naive Bayes classifier is trained on only the region of the target sequences that was sequenced.'
    @params['forward_primer'] = [ "DAGAGTTTGATCMTGGCTCAG", "GAGAGTTTGATYMTGGCTCAG", "GATCCTACGGGAGGCAGCA", "GTGCCAGCMGCCGCGGTAA", "CCTACGGGNGGCWGCAG", "GTGCCAGCMGCCGCGGTAA", "CCTACGGGAGGCAGCAG", "GTGCCAGCMGCNGCGG3"]
    @params['forward_primer', 'description'] = "Choose accoding to the same order as the primer parameter"
    @params['reverse_primer'] = [ "TMTTACCGCGGCNGCTGGCAC", "ACCGCGGCTGCTGGCAC", "CTTACCGCGGCTGCTGGC", "GGACTACHVGGGTWTCTAAT", "GACTACHVGGGTATCTAATCC", "CCGTCAATTCMTTTRAGTTT", "CCGTCAATTCMTTTRAGT", "GGGTTNCGNTCGTTG"]
    @params['reverse_primer', 'description'] = "Choose accoding to the same order as the primer parameter"
    @params['name'] = 'QIIME2'
    @params['mail'] = ""
    @inherit_tags = ['B-Fabric']
    @modules = ['Dev/R']
  end
   def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
 def set_default_parameters
    @params['paired'] = dataset_has_column?('Read2')
  end
  def next_dataset
     nds = {'Name'=>@params['name']}
     nds['ResultDir [File]'] = File.join(@result_dir, 'Results_Folder/')
     nds['Demux Report [Link]'] = File.join(@result_dir, 'Results_Folder/demux_seqs.qzv.zip.folder/data/index.html')
     nds['Denoising stats [Link]'] = File.join(@result_dir, 'Results_Folder/dada2_denoising_stats.qzv.zip.folder/data/index.html')
     nds['Feature Table [Link]'] = File.join(@result_dir, 'Results_Folder/table.qzv.zip.folder/data/index.html')
     nds['Rep Seqs Report [Link]'] = File.join(@result_dir, 'Results_Folder/dada2_rep_set.qzv.zip.folder/data/index.html')
     nds['Taxonomy Barplot [Link]'] = File.join(@result_dir, 'Results_Folder/taxa-bar-plots.qzv.zip.folder/data/index.html')
     nds['Taxonomy List [Link]'] = File.join(@result_dir, 'Results_Folder/taxonomy.qzv.zip.folder/data/index.html')
     nds['Shannon Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/shannon_group_significance.qzv.zip.folder/data/index.html')
     nds['Jaccard Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/jaccard_group_significance.qzv.zip.folder/data/index.html')
     nds['Bray Curtis Diversity [Link]'] = File.join(@result_dir, 'Results_Folder/bray_curtis_group_significance.qzv.zip.folder/data/index.html')
     nds['Jaccard Emperor Plot [Link]'] = File.join(@result_dir, 'Results_Folder/jaccard_emperor_plot.qzv.zip.folder/data/index.html')
     nds['Bray Curtis Emperor Plot [Link]'] = File.join(@result_dir, 'Results_Folder/bray_curtis_emperor_plot.qzv.zip.folder/data/index.html')
     nds['Alpha rarefaction [Link]'] = File.join(@result_dir, 'Results_Folder/alpha-rarefaction.qzv.zip.folder/data/index.html')
     nds['Differential abundace [Link]'] = File.join(@result_dir, 'Results_Folder/ancom_group.qzv.zip.folder/data/index.html')
     nds
  end
  def commands
     run_RApp("EzAppQIIME2", conda_env: "qiime2-2022.2")
  end
end

if __FILE__ == $0

end
