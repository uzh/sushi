#!/usr/bin/env ruby
# encoding: utf-8
Version = '20210414-121343'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class TrimmomaticApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'Trimmomatic'
    @analysis_category = 'Prep'
    @description =<<-EOS
Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.
Refer to <a href='http://www.usadellab.org/cms/?page=trimmomatic'>http://www.usadellab.org/cms/?page=trimmomatic</a>
    EOS
    @required_columns = ['Name','Read1']
    @required_params = ['paired', 'quality_type']
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '15'
    @params['scratch'] = '100'
    @params['paired'] = false
    @params['paired', 'description'] = 'either the reads are paired-ends or single-end'
    @params['quality_type'] = ['phred33', 'phred64']
    @params['quality_type', 'description'] = 'Fastq quality score type, if you use Illumina HiSeq or MySeq, chose phred33'
    @params['illuminaclip'] = {
      'none' => '',
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
    @params['crop'] = ''
    @params['crop', 'description'] = 'The number of bases to keep, from the start of the read'
    @params['headcrop'] = '0'
    @params['headcrop', 'description'] = 'The number of bases to remove from the start of the read'
    @params['minlen'] = '30'
    @params['minlen', 'description'] = 'Drop the read if it is below a specified length'
    @params['mail'] = ""
    @modules = ["Dev/jdk", "QC/Trimmomatic"]
    @inherit_tags = ["Factor"]
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
   dataset =  {'Name'=>@dataset['Name'],
    'Read1 [File]' => File.join(@result_dir, "#{File.basename(@dataset['Read1'].to_s).gsub('fastq.gz','trimmed.fastq.gz')}"),
    'Adapters [File]' => File.join(@result_dir, "#{@dataset['Name']}_adapters.fa"),
    'Read Count' => 0
    }.merge(extract_columns(@inherit_tags))
  if @params['paired'] 
      dataset['Read2 [File]'] = File.join(@result_dir, "#{File.basename(@dataset['Read2'].to_s).gsub('fastq.gz','trimmed.fastq.gz')}")
  end
  dataset
  end
  def se_pe
    @params['paired'] ? 'PE' : 'SE'
  end
  def commands
    command = ""
    adapters_fa = "#{@dataset['Name']}_adapters.fa"
    if @dataset['Adapter1']
      command << "echo '>Adapter1' > #{adapters_fa}\n"
      command << "echo '#{@dataset["Adapter1"]}' >> #{adapters_fa}\n"
      if @dataset['Adapter2']
        command << "echo '>Adapter2' >> #{adapters_fa}\n"
        command << "echo #{@dataset['Adapter2']} >> #{adapters_fa}\n"
      end
    end
    unless @params['illuminaclip'].to_s.empty?
        command << "cat #{@params['illuminaclip']} >> #{adapters_fa}\n"
    end
    command << "java -jar $Trimmomatic_jar #{se_pe} -threads #{@params['cores']} -#{@params['quality_type']} #{File.join(SushiFabric::GSTORE_DIR, @dataset['Read1'])}"
    if @params['paired']
      command << " #{File.join(SushiFabric::GSTORE_DIR, @dataset['Read2'])}"
    end
    output_R1 = File.basename(@dataset['Read1']).gsub('fastq.gz', 'trimmed.fastq.gz')
    command << " #{output_R1}"
    if @params['paired']
      output_unpared_R1 = File.basename(@dataset['Read1']).gsub('fastq.gz', 'unpaired.fastq.gz')
      command << " #{output_unpared_R1}"
    end
    if @params['paired']
      output_R2 = File.basename(@dataset['Read2']).gsub('fastq.gz', 'trimmed.fastq.gz')
      output_unpared_R2 = File.basename(@dataset['Read2']).gsub('fastq.gz', 'unpaired.fastq.gz')
      command << " #{output_R2} #{output_unpared_R2}"
    end
    command << " ILLUMINACLIP:#{adapters_fa}:#{@params['seed_mismatchs']}:#{@params['palindrome_clip_threshold']}:#{@params['simple_clip_threshold']}"
    unless @params['leading'].strip.empty?
      command << " LEADING:#{@params['leading']}"
    end
    unless @params['trailing'].strip.empty?
      command << " TRAILING:#{@params['trailing']}"
    end
    unless @params['slidingwindow'].strip.empty?
      command << " SLIDINGWINDOW:#{@params['slidingwindow']}"
    end
    unless @params['avgqual'].strip.empty?
      command << " AVGQUAL:#{@params['avgqual']}"
    end
    unless @params['crop'].strip.empty?
      command << " CROP:#{@params['crop']}"
    end
    unless @params['headcrop'].strip.empty?
      command << " HEADCROP:#{@params['headcrop']}"
    end
    unless @params['minlen'].strip.empty?
      command << " MINLEN:#{@params['minlen']}"
    end
    command
  end
end

if __FILE__ == $0

end

