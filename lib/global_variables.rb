#!/usr/bin/env ruby
# encoding: utf-8

def c(*list)
  list
end
GLOBAL_VARIABLES = '/usr/local/ngseq/opt/sushi_scripts/GLOBAL_VARIABLES.txt'
load GLOBAL_VARIABLES
#R_SCRIPT_DIR="/usr/local/ngseq/opt/sushi_scripts"

module GlobalVariables
  SUSHI = 'Supercalifragilisticexpialidocious!!'
  def builder_selector(base_dir, shown_pattern=nil, value_pattern=nil)
    selector = {}
    Dir[base_dir].sort.select{|dir| File.directory?(dir)}.each do |dir|
      key = if shown_pattern
              dir.gsub(shown_pattern.keys.first, shown_pattern)
            else
              dir
            end
      value = if value_pattern
                dir.gsub(value_pattern.keys.first, value_pattern)
              else
                File.basename(dir)
              end
      selector[key] = value
    end
    selector
  end
  def ref_selector
    selector = {'select'=>''}

    base_pattern = "/srv/GT/reference/*/*/*"
    shown_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    value_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    builds = builder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

    base_pattern = "/srv/GT/reference/*/*/*/Annotation/Version*"
    versions = builder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

    builds.keys.each do |build_key|
      if versions.keys.find{|version_key| version_key=~/#{build_key}/}
        builds.delete(build_key)
      end
    end
    builds.merge!(versions)
    builds = builds.sort.to_h

    selector.merge!(builds)
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
  def run_RApp(app_name = self.class.to_s[0].downcase+self.class.to_s[1,20])
    command = "/usr/local/ngseq/bin/R --vanilla --slave<<  EOT\n"
    command << "GLOBAL_VARIABLES <<- '#{GLOBAL_VARIABLES}'\n"
    command << "R_SCRIPT_DIR <<- '#{R_SCRIPT_DIR}'\n"
    command<<  "source(file.path(R_SCRIPT_DIR, 'init.R'))\n"
    command << "config = list()\n"
    config = @params
    config.keys.each do |key|
      command << "config[['#{key}']] = '#{config[key]}'\n"
    end
    command << "config[['dataRoot']] = '#{@gstore_dir}'\n"
    command << "config[['resultDir']] = '#{@result_dir}'\n"
    command << "output = list()\n"
    output = next_dataset
    output.keys.each do |key|
      command << "output[['#{key}']] = '#{output[key]}'\n"
    end
    command<<  "inputDatasetFile = '#{@input_dataset_tsv_path}'\n"
    command<<  "runApp('#{app_name}', input=inputDatasetFile, output=output, config=config)\n"
    command<<  "EOT\n"
    command
  end

end
