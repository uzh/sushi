#!/usr/bin/env ruby
# encoding: utf-8

GENOME_REF_DIRS = ['/srv/GT/reference', '/srv/GT/assembly']
GENOME_REF_DIRS_FAVORITE = ['/srv/GT/reference-favorite']
EZ_GLOBAL_VARIABLES = '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
#EZ_GLOBAL_VARIABLES = '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES_DEMO.txt'

# compatible method to R vector, c(x,x,x) 
def c(*list)
  list
end
if File.exist?(EZ_GLOBAL_VARIABLES)
  FALSE = false
  load EZ_GLOBAL_VARIABLES
end

module GlobalVariables
  SUSHI = 'Supercalifragilisticexpialidocious!!'
  def fetch_genome_root(reference)
    selector = ref_selector({}, 'fetch_genome_root')
    if full_path = selector[reference]
      genome_root = full_path.gsub(/#{reference}/, '')
    else
      raise "Reference is not found: #{reference}"
    end
    genome_root.delete_suffix("/")
  end
  def refBuilder_selector(base_dirs, shown_pattern=nil, value_pattern=nil)
    selector = {}
    Dir[*base_dirs].sort.select{|dir| File.directory?(dir)}.each do |dir|
      key = if shown_pattern
              dir.gsub(/#{shown_pattern.keys.join("|")}/, shown_pattern)
            else
              dir
            end
      value = if value_pattern
                dir.gsub(/#{value_pattern.keys.join("|")}/, value_pattern)
              else
                File.basename(dir)
              end
      selector[key] = value
    end
    selector
  end
  def ref_selector(value_replace_pattern=GENOME_REF_DIRS.inject({}){|hash, dir| hash[dir+"/"]=''; hash}, cache_name='ref_selector')
    selector = {'select'=>''}
    selector = Rails.cache.fetch(cache_name, expired_in: 1.hour) do
      base_pattern_favorite = GENOME_REF_DIRS_FAVORITE.map{|dir| "#{dir}/*/*/*/Annotation/Release*"}.flatten
      shown_replace_pattern_favorite = GENOME_REF_DIRS_FAVORITE.inject({}){|hash, dir| hash[dir+"/"]=''; hash}
      value_replace_pattern_favorite =GENOME_REF_DIRS_FAVORITE.inject({}){|hash, dir| hash[dir+"/"]=''; hash}
      refBuilds_favorites = refBuilder_selector(base_pattern_favorite, shown_replace_pattern_fovarite, value_replace_pattern_forvarite)
      separator = {'-'*50 => '-'}
      refBuilds_favorites.merge!(separator)

      base_pattern = GENOME_REF_DIRS.map{|dir| "#{dir}/*/*/*"}
      shown_replace_pattern = GENOME_REF_DIRS.inject({}){|hash, dir| hash[dir+"/"]=''; hash}
      refBuilds = refBuilder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

      base_pattern = GENOME_REF_DIRS.map{|dir| ["#{dir}/*/*/*/Annotation/Version*", "#{dir}/*/*/*/Annotation/Release*"]}.flatten
      versions = refBuilder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

      base_pattern = GENOME_REF_DIRS.map{|dir| "#{dir}/*/*/*/Sequence"}
      sequences = refBuilder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

      refBuilds.keys.each do |refBuild_key|
        if versions.keys.find{|version_key| version_key=~/#{refBuild_key}/}
          refBuilds.delete(refBuild_key)
        end
        unless sequences.keys.find{|sequence_key| sequence_key=~/#{refBuild_key}/}
          refBuilds.delete(refBuild_key)
        end
      end
      refBuilds.merge!(versions)
      refBuilds = refBuilds.sort.to_h
      refBuilds = refBuilds_favorites.merge(refBuilds)

      selector.merge!(refBuilds)
      selector
    end
    selector
  end
  def extract_column(type)
    factors = get_columns_with_tag(type)
    dataset = {}
    factors.first.keys.each do |colname|
      dataset[colname+" [#{type}]"] = @dataset[colname]
    end
    dataset
  rescue
    {}
  end
  def extract_columns(tags)
    dataset = {}
    tags.each do |tag|
      additional = extract_column(tag)
      dataset.merge!(additional)
    end
    dataset
  end
  def factor_dataset
    # deprecated
    # please use extract_column() method
    extract_column("Factor")
  end
  def run(class_name)
    opt = OptionParser.new do |o|
      o.banner = "Usage: ruby #{__FILE__} [options]"
      o.on(:dataset_id, '-i dataset_id', '--dataset_id', Integer, 'DataSet ID in Sushi DB')
      o.on(:dataset_tsv, '-d dataset_tsv', '--dataset', String, 'DataSet file (.tsv) (This option is prior to dataset_id option)')
      o.on(:parameterset_tsv, '-m parameterset_tsv', '--parameterset', String, 'Parameterset file (.tsv)')
      o.on(:run_mode, '-r', '--run', 'Real run mode. without this option, it runs with test run mode which checks only DataSet and Parameters and no submittion')
      o.on(:project, '1001', '-p', '--project', String, 'Project Number (default: 1001)')
      o.on(:user, 'sushi_lover', '-u', '--user', String, 'Submit user (default: sushi_lover)')
      o.parse!(ARGV)
    end
    opt.project = 'p' + opt.project
    unless opt.dataset_id or opt.dataset_tsv
      puts
      warn "\e[31mERROR\e[0m: Either dataset_id or dataset_tsv is required"
      print opt.help
      exit
    end

    usecase = class_name.new

    usecase.project = opt.project
    usecase.user = opt.user
    usecase.dataset_sushi_id = opt.dataset_id
    usecase.dataset_tsv_file = opt.dataset_tsv
    usecase.parameterset_tsv_file = opt.parameterset_tsv
    if opt.run_mode
      usecase.run
    else
      usecase.test_run
    end
  end
  def run_RApp(app_name = self.class.to_s[0].downcase+self.class.to_s[1,20], lib_path: nil, conda_env: nil)
    command = ''
    if conda_env
      command << ". '/usr/local/ngseq/miniconda3/etc/profile.d/conda.sh'\n"
      command << "conda activate #{conda_env}\n"
    end
    command << "R --vanilla --slave<<  EOT\n"
    command << "EZ_GLOBAL_VARIABLES <<- '#{EZ_GLOBAL_VARIABLES}'\n"
    if lib_path
      command << ".libPaths('#{lib_path}')\n"
    end
    command<<  "if (!library(ezRun, logical.return = TRUE)){\n"
    command<<  "message('retry loading ezRun')\n"
    command<<  "Sys.sleep(120)\n"
    command<<  "library(ezRun)\n"
    command<<  "}\n"
    command << "param = list()\n"
    param = @params
    param.keys.each do |key|
      command << "param[['#{key}']] = '#{param[key]}'\n"
    end
    command << "param[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "param[['resultDir']] = '#{@result_dir}'\n"
    command << "param[['isLastJob']] = #{@last_job.to_s.upcase}\n"
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n"
    end
    if @params['process_mode'] == 'DATASET'
      command <<  "input = '#{@input_dataset_tsv_path}'\n"
    else # sample mode
      command << "input = list()\n" 
      input = @dataset
      input.keys.each do |key|
        command << "input[['#{key}']] = '#{input[key]}'\n" 
      end
    end
    command<<  "#{app_name}\\$new()\\$run(input=input, output=output, param=param)\n"
    command<<  "EOT\n"
    command
  end

end
