#!/usr/bin/env ruby
# encoding: utf-8
Version = '20220210-150606'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class MACS3App < SushiFabric::SushiApp
  def initialize
    super
    @name = 'MACS3'
    @analysis_category = 'Peaks'
    @description =<<-EOS
Capturing the influence of genome complexity to evaluate the significance of enriched ChIP regions<br/>
<a href='https://macs3-project.github.io/MACS/'>https://macs3-project.github.io/MACS/</a>
    EOS
    @required_columns = ['Name','BAM','BAI', 'refBuild']
    @required_params = ['refBuild','paired']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '40'
    @params['scratch'] = '100'
    @params['refBuild'] = ref_selector
    @params['paired'] = false
    @params['useControl'] = false
    @params['useControl', 'description'] = 'Input or other controls needs to be specified in the input_dataset with a separate Control [File] column linking to the bam file'
    @params['shiftATAC'] = false
    @params['shiftATAC', 'description'] = 'should all reads aligning to + strand were offset by +4bp, all reads aligning to the - strand are offset -5 bp'
    @params['annotatePeaks'] = true
    @params['refFeatureFile'] = 'genes.gtf'
    @params['mode'] = ['ChIP-seq', 'ATAC-seq']
    @params['mode', 'description'] = 'Call MACS3 for ChIP-seq or ATAC-seq data.'
    @params['cmdOptions'] = '--nomodel --bw 200 --extsize 147'
    @params['cmdOptions', 'description'] = 'default options deactivate model building and extend reads in 5p->3p direction to fix-sized fragments.'
    @params['genomeSize'] = '0'
    @params['genomeSize', 'description'] = 'if you keep 0 the genome size is calculated from the given reference genome'
    @params['removeDuplicates'] = true
    @params['removeDuplicates', 'description'] = 'remove duplicates before running MACS3. The resulting processed.bam does not contain duplicates anymore.'
    @params['specialOptions'] = ''
    @params['mail'] = ''
    @modules = ["Tools/samtools", "Tools/UCSC", "Tools/BEDTools", "Dev/R", "Tools/Picard"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    @params['refBuild'] = @dataset[0]['refBuild']
    if dataset_has_column?('refFeatureFile')
      @params['refFeatureFile'] = @dataset[0]['refFeatureFile']
    end
    if dataset_has_column?('paired')
      @params['paired'] = @dataset[0]['paired'].to_s.downcase
    end
 end

  def next_dataset
    bw_link = File.join(@result_dir, "#{@dataset['Name']}_processed.bw")
    bam_link = File.join(@result_dir, "#{@dataset['Name']}_processed.bam")
    bai_link = File.join(@result_dir, "#{@dataset['Name']}_processed.bam.bai")
    
    peakfile_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.xls")
    bedfile_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.bed")
    peakseq_link = File.join(@result_dir, "#{@dataset['Name']}_peaks.fa")

    {'Name'=>@dataset['Name'],
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'refFeatureFile'=>@params['refFeatureFile'],
     'paired'=>@params['paired'],
     'CalledPeaks [File]'=>peakfile_link,
     'BED [File]'=>bedfile_link,
     'PeakSequences [File]'=>peakseq_link,
     'BigWigFile [File]'=>bw_link,
     'BAM [File]'=>bam_link,
     'BAI [File]'=>bai_link
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    '. "/usr/local/ngseq/miniforge3/etc/profile.d/conda.sh"' +
    "\nconda activate gi_macs3\n" +
    run_RApp("EzAppMacs3")
  end
end

if __FILE__ == $0

end
