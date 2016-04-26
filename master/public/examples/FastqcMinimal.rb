require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class FastqcMinimal < SushiFabric::SushiApp
  def initialize
    super
    @name = 'FastqcMinimal'
    @analysis_category = 'QC'
    @required_columns = ['Name','Read1']
    @required_params = ['cores']
    @params['cores'] = 1
  end
  def next_dataset
    {
     'Name'=>@dataset['Name'], 
     'Report [Link, File]'=>File.join(@result_dir, File.basename(@dataset['Read1'].to_s).gsub('.fastq.gz', '_fastqc.zip'))
    }
  end
  def commands
    <<-EOS
which fastqc
fastqc -v
fastqc --extract -o . -t #{@params['cores']} #{@gstore_dir}/#{@dataset['Read1']}
    EOS
  end
end
