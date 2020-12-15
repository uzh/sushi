#!/usr/bin/env ruby
# encoding: utf-8
Version = '20171109-094054'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class DiscovarDenovoApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DiscovarDenovo'
    @analysis_category = 'Assemble'
    @params['process_mode'] = 'DATASET'
    @description =<<-EOS
DISCOVAR de novo is a large (and small) de novo genome assembler. It quickly generates highly accurate and complete assemblies using the same single library data as used by DISCOVAR. It currently doesn’t support variant calling – for that, please use DISCOVAR instead.
Refer to <a href='https://software.broadinstitute.org/software/discovar/blog/'>https://software.broadinstitute.org/software/discovar/blog/</a>

As a preprocess, Trimmomatic can run.
Refer to <a href='http://www.usadellab.org/cms/?page=trimmomatic'>http://www.usadellab.org/cms/?page=trimmomatic</a>
    EOS
    @required_columns = ['Name','Read1']
    #@required_params = ['paired', 'quality_type', 'leading', 'trailing', 'slidingwindow', 'avgqual', 'headcrop', 'minlen']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'

    # for Trimmomatic
    @params['trimming'] = true
    @params['trimming', 'description'] = 'if you set this false, the following options will be ignored'
    @params['paired'] = false
    @params['paired', 'description'] = 'either the reads are paired-ends or single-end'
    @params['quality_type'] = ['phred33', 'phred64']
    @params['quality_type', 'description'] = 'Fastq quality score type, if you use Illumina HiSeq or MySeq, chose phred33'
    @params['illuminaclip'] = {
      'allIllumina-forTrimmomatic-20160202.fa' => '/srv/GT/databases/contaminants/allIllumina-forTrimmomatic-20160202.fa',
      'All Illumina Adapter' => '/srv/GT/databases/contaminants/illuminaContaminants.fa',
      'FastQC checking Adapter' => '/srv/GT/databases/adapter/adapter_list.fa',
      #'Original Adapter' => 'under construction'
    }
    @params['seed_mismatchs'] = '1'
    @params['seed_mismatchs', 'description'] = 'Specifies the maximum mismatch count which will still allow a full match to be performed'
    @params['palindrome_clip_threshold'] = '30'
    @params['palindrome_clip_threshold', 'description'] = 'THIS IS NOT USED. In order to use this method, adapter name should start with Prefix in the fast file. Specifies how accurate the match between the two _adapter ligated_ reads must be for PE palindrome read alignment'
    @params['simple_clip_threshold'] = '10'
    @params['simple_clip_threshold', 'description'] = 'Specifies how accurate the match between any adapter etc. sequence must be against a read. It will change depending on the adapter length. For FastQC checking adapter, 7 may be proper.'
    @params['leading'] = '5'
    @params['leading', 'description'] = 'Cut bases off the start of a read, if below a threshold quality'
    @params['trailing'] = '5'
    @params['trailing', 'description'] = 'Cut bases off the end of a read, if below a threshold quality'
    @params['slidingwindow'] = '5:15'
    @params['slidingwindow', 'description'] = 'Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold'
    @params['avgqual'] = '20'
    @params['avgqual', 'description'] = 'Drop the read if the average quality is below the specified level'
    @params['headcrop'] = 0
    @params['headcrop', 'description'] = 'The number of bases to remove from the start of the read'
    @params['minlen'] = '30'
    @params['minlen', 'description'] = 'Drop the read if it is below a specified length'

    @modules = ["QC/Trimmomatic", "Assembly/DiscovarDenovo"]
    #@inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
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
    {'Name'=>"Assembly", 
     'Results [File]'=>File.join(@result_dir, "discovar_out"),	
    }
  end
  def se_pe
    @params['paired'] ? 'PE' : 'SE'
  end
  def commands
    command = ""
    # Trimmomatic
    if @params['trimming']
      @dataset_hash.each do |row|
        dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
        p dataset
        adapters_fa = "#{dataset['Name']}_adapters.fa"
        if dataset['Adapter1']
          command << "echo '>Adapter1' > #{adapters_fa}\n"
          command << "echo '#{dataset["Adapter1"]}' >> #{adapters_fa}\n"
          if dataset['Adapter2']
            command << "echo '>Adapter2' >> #{adapters_fa}\n"
            command << "echo #{dataset['Adapter2']} >> #{adapters_fa}\n"
          end
        end
        unless @params['illuminaclip'].to_s.empty?
            command << "cat #{@params['illuminaclip']} >> #{adapters_fa}\n"
        end
        command << "java -jar $Trimmomatic_jar #{se_pe} -threads #{@params['cores']} -#{@params['quality_type']} #{File.join(SushiFabric::GSTORE_DIR, dataset['Read1'])}"
        if @params['paired']
          command << " #{File.join(SushiFabric::GSTORE_DIR, dataset['Read2'])}"
        end
        output_R1 = File.basename(dataset['Read1']).gsub('fastq.gz', 'trimmed.fastq.gz')
        command << " #{output_R1}"
        if @params['paired']
          output_unpared_R1 = File.basename(dataset['Read1']).gsub('fastq.gz', 'unpaired.fastq.gz')
          command << " #{output_unpared_R1}"
        end
        if @params['paired']
          output_R2 = File.basename(dataset['Read2']).gsub('fastq.gz', 'trimmed.fastq.gz')
          output_unpared_R2 = File.basename(dataset['Read2']).gsub('fastq.gz', 'unpaired.fastq.gz')
          command << " #{output_R2} #{output_unpared_R2}"
        end
        command << " ILLUMINACLIP:#{adapters_fa}:#{@params['seed_mismatchs']}:#{@params['palindrome_clip_threshold']}:#{@params['simple_clip_threshold']}"
        command << " LEADING:#{@params['leading']} TRAILING:#{@params['trailing']} SLIDINGWINDOW:#{@params['slidingwindow']} AVGQUAL:#{@params['avgqual']} HEADCROP:#{@params['headcrop']} MINLEN:#{@params['minlen']}\n"
      end
    end

    # DiscvarDenovo
    command << "mkdir discovar_out\n"
    file_cols = get_columns_with_tag("File")
    reads = if @params['trimming']
              file_cols.map{|data_set| data_set.values.map{|file| File.basename(file).gsub(/.fastq.gz/, '.trimmed.fastq.gz')}}.flatten.join(",")
            else
              file_cols.map{|data_set| data_set.values.map{|file| File.join(@gstore_dir, file)}}.flatten.join(",")
            end
    command << "DiscovarDeNovo READS=#{reads} OUT_DIR=discovar_out NUM_THREADS=#{@params['cores']} MAX_MEM_GB=#{@params['ram']}\n"

    command
  end
end

if __FILE__ == $0

end

