#!/usr/bin/env ruby
# encoding: utf-8
Version = '20131011-153152'

require 'csv'
require 'fileutils'
require 'active_record'
require 'yaml'

CONFIG = 'sushi_configure.yml'
current_dir = File.dirname(File.expand_path(__FILE__))
config_yml = File.join(current_dir, CONFIG)
config = if File.exist?(config_yml)
           YAML.load(File.read(config_yml))
         else
           {}
         end
WORKFLOW_MANAGER = config[:workflow_manager]||'druby://localhost:3000'
GSTORE_DIR = config[:gstore_dir]||'gstore'
sushi_app_dir = File.expand_path('../..', __FILE__)
SUSHI_APP_DIR = config[:sushi_app_dir]||sushi_app_dir

begin
  ::Project
rescue
  ActiveRecord::Base.establish_connection(
              :adapter  => 'sqlite3',
              :database => "#{SUSHI_APP_DIR}/db/development.sqlite3" 
          )
  require "#{SUSHI_APP_DIR}/app/models/project"
  require "#{SUSHI_APP_DIR}/app/models/data_set"
  require "#{SUSHI_APP_DIR}/app/models/sample"
end

def save_data_set(data_set_arr, headers, rows)
  data_set_hash = Hash[*data_set_arr]
  if project = Project.find_by_number(data_set_hash['ProjectNumber'].to_i)
    data_set = DataSet.new
    data_set.name = data_set_hash['DataSetName']
    data_set.project = project
    if parent_id = data_set_hash['ParentID'] and parent_data_set = DataSet.find_by_id(parent_id.to_i)
      data_set.data_set = parent_data_set
    end
    if comment = data_set_hash['Comment'] and !comment.to_s.empty?
      data_set.comment = comment
    end

    sample_hash = {}
    rows.each do |row|
      headers.each_with_index do |header, i|
       sample_hash[header]=row[i]
      end
      sample = Sample.new
      sample.key_value = sample_hash.to_s
      sample.save unless sample.saved?
      data_set.samples << sample
    end

    data_set.md5 = data_set.md5hexdigest
    unless data_set.saved?
      project.data_sets << data_set
      parent_data_set.data_sets << data_set if parent_data_set
      data_set.save
    end
    data_set.id
  end
end

class Hash
  attr_reader :defaults
  alias :set :[]=
  alias :get :[]
  def []=(k1,k2,v=nil)
    if v
      @desc ||= {}
      @desc.set([k1,k2].join('_'),v)
    else
      @defaults ||= {}
      if !@defaults[k1] and k2
        if k2.instance_of?(Array)
          @defaults.set(k1,k2.first)
        elsif k2.instance_of?(Hash) and k2.first
          @defaults.set(k1,k2.first.last)
        else
          @defaults.set(k1,k2)
        end
      end
      set(k1,k2)
    end
  end
  def default_value(k,v=nil)
    if v
      @defaults[k] = v
    else
      @defaults[k]
    end
  end
  def data_type(k)
    @defaults[k].class
  end
  def data_types
    Hash[@defaults.map{|k,v| [k, v.class]}]
  end
  def [](k1, k2=nil)
    if k2
      if @desc
        @desc.get([k1,k2].join('_'))
      else
        nil
      end
    else
      get(k1)
    end
  end
end
class String
  def tag?(tag)
    scan(/\[(.*)\]/).flatten.join =~ /#{tag}/
  end
end
class SushiApp
  attr_reader :params
  attr_reader :job_ids
  attr_reader :next_dataset_id
  attr_reader :required_columns
  attr_reader :required_params
  attr_reader :dataset_hash
  attr_reader :analysis_category
  attr_reader :description
  attr_reader :name
  attr_accessor :dataset_tsv_file
  attr_accessor :parameterset_tsv_file
  attr_accessor :dataset_sushi_id
  attr_accessor :project
  attr_accessor :user
  attr_accessor :next_dataset_name
  attr_accessor :next_dataset_comment
  def initialize
    @gstore_dir = GSTORE_DIR
    @project = nil
    @name = nil
    @params = {}
    @params['cores'] = nil
    @params['ram'] = nil
    @params['scratch'] = nil
    @params['node'] = ''
    @params['process_mode'] = 'SAMPLE'
    @job_ids = []
    @required_columns = []
  end
  def set_input_dataset
    if @dataset_tsv_file
      dataset_tsv = CSV.readlines(@dataset_tsv_file, {:headers=>true, :col_sep=>"\t"})
      @dataset_hash = []
      dataset_tsv.each do |row|
        @dataset_hash << row.to_hash
      end
    elsif @dataset_sushi_id
      @dataset_hash = []
      if dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
        dataset.samples.each do |sample|
          @dataset_hash << sample.to_hash
        end
      end
    end
    @dataset_hash
  end
  def get_columns_with_tag(tag)
    #@factor_cols = @dataset_hash.first.keys.select{|header| header =~ /\[#{tag}\]/}.map{|header| header.gsub(/\[.+\]/,'').strip}
    @dataset_hash.map{|row| 
      Hash[*row.select{|k,v| k=~/\[#{tag}\]/}.map{|k,v| [k.gsub(/\[.+\]/,'').strip,v]}.flatten]
    }
  end
  def set_output_files
    @dataset = {}
    next_dataset.keys.select{|header| header.tag?('File')}.each do |header|
      @output_files ||= []
      @output_files << header
    end
    @output_files = @output_files.uniq
  end
  def check_required_columns
    if @dataset_hash and @required_columns and (@required_columns-@dataset_hash.map{|row| row.keys}.flatten.uniq.map{|colname| colname.gsub(/\[.+\]/,'').strip}).empty?
      true
    else
      false
    end
  end
  def check_application_parameters
    if @required_params and (@required_params - @params.keys).empty?
      @output_params = @params.clone
    end
  end
  def set_user_parameters
    # this should be done in an instance of applicaiton subclass
    if @parameterset_tsv_file
      parameterset_tsv = CSV.readlines(@parameterset_tsv_file, :col_sep=>"\t")
      headers = []
      parameterset_tsv.each do |row|
        header, value = row
        headers << header
        @params[header] = if @params.data_type(header) == String
                            value
                          else
                            eval(value)
                          end
      end
      (@params.keys - headers).each do |key|
        @params[key] = @params.default_value(key)
      end
    end
    @params
  end
  def set_dir_paths
    ## sushi figures out where to put the resulting dataset
    unless @name and @project
      raise "should set #name and #project"
    end
    @name.gsub!(/\s/,'_')
    result_dir_base = if @next_dataset_name
                        [@next_dataset_name, Time.now.strftime("%Y-%m-%d--%H-%M-%S")].join("_")
                      else
                        [@analysis_category, @name, @dataset_sushi_id.to_s, Time.now.strftime("%Y-%m-%d--%H-%M-%S")].join("_")
                      end
    @result_dir = File.join(@project, result_dir_base)
    @scratch_result_dir = File.join("/scratch", result_dir_base)
    @job_script_dir = File.join(@scratch_result_dir, 'scripts')
    @gstore_result_dir = File.join(@gstore_dir, @result_dir)
    @gstore_script_dir = File.join(@gstore_result_dir, 'scripts')
    @gstore_project_dir = File.join(@gstore_dir, @project)
    set_file_paths
  end
  def prepare_result_dir
    FileUtils.mkdir_p(@scratch_result_dir)
    FileUtils.mkdir_p(@job_script_dir)
  end
  def job_header
    @scratch_dir = if @params['process_mode'] == 'SAMPLE'
                     @scratch_result_dir + "_" + @dataset['Name']
                   else
                     @scratch_result_dir
                   end
    @out.print <<-EOF
#!/bin/bash
set -e
set -o pipefail

#### SET THE STAGE
SCRATCH_DIR=#{@scratch_dir}
GSTORE_DIR=#{@gstore_dir}
mkdir $SCRATCH_DIR || exit 1
cd $SCRATCH_DIR || exit 1
echo "Job runs on `hostname`"
echo "at $SCRATCH_DIR"

    EOF
  end
  def job_footer
    @out.print "#### JOB IS DONE WE PUT THINGS IN PLACE AND CLEAN AUP\n"
    if @output_files
      @output_files.map{|header| next_dataset[header]}.each do |file|
        # in actual case, to save under /srv/gstore/
        src_file = File.basename(file)
        dest_dir = File.dirname(File.join(@gstore_dir, file))
        @out.print copy_commands(src_file, dest_dir).join("\n"), "\n"
      end
    end
    @out.print <<-EOF
cd ~
rm -rf #{@scratch_dir} || exit 1
    EOF

  end
  def job_main
    @out.print "#### NOW THE ACTUAL JOBS STARTS\n"
    @out.print commands, "\n\n"
  end
  def next_dataset
    # this should be overwritten in a subclass
  end
  def commands
    # this should be overwritten in a subclass
  end
  def submit_command(job_script)
    gsub_options = []
    gsub_options << "-c #{@params['cores']}" unless @params['cores'].to_s.empty?
    gsub_options << "-n #{@params['node']}" unless @params['node'].to_s.empty?
    gsub_options << "-r #{@params['ram']}" unless @params['ram'].to_s.empty?
    gsub_options << "-s #{@params['scratch']}" unless @params['scratch'].to_s.empty?
    gsub_options << "-u #{@user}" if @user
    command = "wfm_monitoring --server #{WORKFLOW_MANAGER} --project #{@project.gsub(/p/,'')} --logdir #{@gstore_script_dir} #{job_script} #{gsub_options.join(' ')}"
  end
  def submit(job_script)
    command = submit_command(job_script)
    puts "submit: #{command}"
    job_id = `#{command}`
    job_id = job_id.to_i
    unless job_id.to_i > 1
      raise 'failed in job submitting'
    end
    job_id
  end
  def preprocess
    # this should be overwritten in a subclass
  end
  def set_file_paths
    @parameter_file = 'parameters.tsv'
    @input_dataset_file = 'input_dataset.tsv'
    @next_dataset_file = 'dataset.tsv'
    @input_dataset_tsv_path = File.join(@gstore_result_dir, @input_dataset_file)
    @parameters_tsv_path = File.join(@gstore_result_dir, @input_dataset_file)
    @next_dataset_tsv_path = File.join(@gstore_result_dir, @next_dataset_file)
  end
  def save_parameters_as_tsv
    file_path = File.join(@scratch_result_dir, @parameter_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      @output_params.each do |key, value|
        out << [key, value]
      end
    end
    file_path
  end
  def save_input_dataset_as_tsv
    file_path = File.join(@scratch_result_dir, @input_dataset_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      headers = @dataset_hash.map{|row| row.keys}.flatten.uniq
      out << headers
      @dataset_hash.each do |row|
        out << headers.map{|header| row[header]}
      end
    end
    file_path
  end
  def save_next_dataset_as_tsv
    headers = @result_dataset.map{|row| row.keys}.flatten.uniq
    file_path = File.join(@scratch_result_dir, @next_dataset_file)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << headers
      @result_dataset.each do |row_hash|
        out << headers.map{|header| row_hash[header]}
      end
    end
    file_path
  end
  def copy_commands(org_dir, dest_parent_dir)
    commands = []
    if @gstore_dir =~ /srv\/gstore/
      commands << "g-req -w copy #{org_dir} #{dest_parent_dir}"
    else
      commands << "mkdir -p #{dest_parent_dir}"
      commands << "cp -r #{org_dir} #{dest_parent_dir}"
    end
    commands
  end
  def copy_dataset_parameter_jobscripts
    org = @scratch_result_dir
    dest = @gstore_project_dir
    copy_commands(org, dest).each do |command|
      puts command
      unless system command
        raise "fails in copying next_dataset files from /scratch to /gstore"
      end
    end
    sleep 1
    command = "rm -rf #{@scratch_result_dir}"
    `#{command}`
  end
  def make_job_script
    @out = open(@job_script, 'w')
    job_header
    job_main
    job_footer
    @out.close
  end
  def sample_mode
    @dataset_hash.each do |row|
      @dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
      ## WRITE THE JOB SCRIPT
      sample_name = @dataset['Name']||@dataset.first
      @job_script = if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
                      File.join(@job_script_dir, @analysis_category + '_' + sample_name) + '_' + dataset.name.gsub(/\s+/,'_') + '.sh'
                    else 
                      File.join(@job_script_dir, @analysis_category + '_' + sample_name) + '.sh'
                    end
      make_job_script
      @job_scripts << @job_script
      @result_dataset << next_dataset
    end
  end
  def dataset_mode
    @job_script = if @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
                    File.join(@job_script_dir, @analysis_category + '_' + dataset.name.gsub(/\s+/,'_') + '.sh')
                  else 
                    File.join(@job_script_dir, @analysis_category + '_' + 'job_script.sh')
                  end
    make_job_script
    @job_scripts << @job_script
    @result_dataset << next_dataset
  end
  def run
    test_run

    ## the user presses RUN
    prepare_result_dir

    ## copy application data to gstore 
    save_parameters_as_tsv
    save_input_dataset_as_tsv


    ## sushi writes creates the job scripts and builds the result data set that is to be generated
    @result_dataset = []
    @job_scripts = []
    if @params['process_mode'] == 'SAMPLE'
      sample_mode
    elsif @params['process_mode'] == 'DATASET'
      dataset_mode
    else 
      #stop
      warn "the process mode (#{@params['process_mode']}) is not defined"
      raise "stop job submitting"
    end

    # job submittion
    @job_scripts.each_with_index do |job_script, i|
      job_id = submit(job_script)
      @job_ids << job_id
      print "Submit job #{File.basename(job_script)} job_id=#{job_id}"
    end

		puts
    print 'job scripts: '
    p @job_scripts
    print 'result dataset: '
    p @result_dataset

    # copy application data to gstore 
    next_dataset_tsv_path = save_next_dataset_as_tsv

    if !@job_ids.empty? and @dataset_sushi_id and dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
      data_set_arr = []
      headers = []
      rows = []
      next_dataset_name = if name = @next_dataset_name
                            name.to_s
                          else
                            "#{@analysis_category}_#{@name.gsub(/\s/,'').gsub(/_/,'')}_#{dataset.id}"
                          end
      data_set_arr = {'DataSetName'=>next_dataset_name, 'ProjectNumber'=>@project.gsub(/p/,''), 'ParentID'=>@dataset_sushi_id, 'Comment'=>@next_dataset_comment.to_s}
      csv = CSV.readlines(next_dataset_tsv_path, :col_sep=>"\t")
      csv.each do |row|
        if headers.empty?
          headers = row
        else
          rows << row
        end
      end
      @next_dataset_id = save_data_set(data_set_arr.to_a.flatten, headers, rows)
    end
    Thread.new do
      copy_dataset_parameter_jobscripts
    end
  end
  def test_run
    set_dir_paths
    set_input_dataset
    preprocess
    set_output_files
    set_user_parameters

    failures = 0
    print 'check project name: '
    unless @project
      puts "\e[31mFAILURE\e[0m: project number is required but not found. you should set it in usecase."
      puts "\tex.)"
      puts "\tapp = #{self.class}.new"
      puts "\tapp.project = 'p1001'"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@project=#{@project}"
    end

    print 'check user name: '
    unless @user
      puts "\e[31mWARNING\e[0m: user number is ought to be added but not found. you should set it in usecase. Default will be 'sushi lover'"
      puts "\tex.)"
      puts "\tapp = #{self.class}.new"
      puts "\tapp.user = 'masa'"
    else
      puts "\e[32mPASSED\e[0m:\n\t@user=#{@user}"
    end

    print 'check application name: '
    if @name.to_s.empty?
      puts "\e[31mFAILURE\e[0m: application name is required but not found. you should set it in application class."
      puts "\tex.)"
      puts "\tclass #{self.class}"
      puts "\t def initialize"
      puts "\t  @name = '#{self.class}'"
      puts "\t end"
      puts "\tend"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@name=#{@name}"
    end

    print 'check analysis_category: '
    if @analysis_category.to_s.empty?
      puts "\e[31mFAILURE\e[0m: analysis_category is required but not found. you should set it in application class."
      puts "\tex.)"
      puts "\tclass #{self.class}"
      puts "\t def initialize"
      puts "\t  @analysis_category = 'Mapping'"
      puts "\t end"
      puts "\tend"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@analysis_category=#{@analysis_category}"
    end

    print 'check dataset: '
    if !@dataset_hash or @dataset_hash.empty?
      puts "\e[31mFAILURE\e[0m: dataset is not found. you should set it by using #{self.class}#dataset_sushi_id or #{self.class}#dataset_tsv_file properties"
      puts "\tex.)"
      puts "\tusecase = #{self.class}.new"
      puts "\tusecase.dataset_tsv_file = \"dataset.tsv\""
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@dataset_hash.length = #{@dataset_hash.length}"
    end

    print 'check required columns: '
    unless check_required_columns
      puts "\e[31mFAILURE\e[0m: required_column(s) is not found in dataset. you should set it in application class."
      puts "\tex.)"
      puts "\tdef initialize"
      puts "\t  @required_columns = ['Name', 'Read1']"
      puts
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end
    puts "\trequired columns: #{@required_columns}"
    puts "\tdataset  columns: #{@dataset_hash.map{|row| row.keys}.flatten.uniq}" if @dataset_hash

    print 'check required parameters: '
    unless check_application_parameters
      puts "\e[31mFAILURE\e[0m: required_param(s) is not set yet. you should set it in usecase"
      puts "\tmissing params: #{@required_params-@params.keys}" if @required_params
      puts "\tex.)"
      puts "\tusecase = #{self.class}.new"
      if @required_params
        puts "\tusecase.params['#{(@required_params-@params.keys)[0]}'] = parameter"
      else
        puts "\tusecase.params['parameter name'] = default_parameter"
      end
      puts
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end
    puts "\tparameters: #{@params.keys}"
    puts "\trequired  : #{@required_params}"

    print 'check next dataset: '
    @dataset={}
    unless self.next_dataset
      puts "\e[31mFAILURE\e[0m: next dataset is not set yet. you should overwrite SushiApp#next_dataset method in #{self.class}"
      puts "\tnote: the return value should be Hash (key: column title, value: value in a tsv table)"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end

    print 'check output files: '
    if !@output_files or @output_files.empty?
      puts "\e[31mWARNING\e[0m: no output files. you will not get any output files after the job running. you can set @output_files (array) in #{self.class}"
      puts "\tnote: usually it should be define in initialize method"
      puts "\t      the elements of @output_files should be chosen from #{self.class}#next_dataset.keys"
      puts "\t      #{self.class}#next_dataset.keys: #{self.next_dataset.keys}" if self.next_dataset
    else
      puts "\e[32mPASSED\e[0m:"
    end

    print 'check commands: '
    if @params['process_mode'] == 'SAMPLE'
      @dataset_hash.each do |row|
        @dataset = Hash[*row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
        unless com = commands
          puts "\e[31mFAILURE\e[0m: any commands is not defined yet. you should overwrite SushiApp#commands method in #{self.class}"
          puts "\tnote: the return value should be String (this will be in the main body of submitted job script)"
          failures += 1
        else
          puts "\e[32mPASSED\e[0m:"
          puts "generated command will be:"
          puts "\t"+com.split(/\n/).join("\n\t")+"\n"
        end
      end
    elsif @params['process_mode'] == 'DATASET'
      unless com = commands
        puts "\e[31mFAILURE\e[0m: any commands is not defined yet. you should overwrite SushiApp#commands method in #{self.class}"
        puts "\tnote: the return value should be String (this will be in the main body of submitted job script)"
        failures += 1
      else
        puts "\e[32mPASSED\e[0m:"
        puts "generated command will be:"
        puts "\t"+com.split(/\n/).join("\n\t")+"\n"
      end
    end

    print 'check workflow manager: '
    begin
      hello = `wfm_hello #{WORKFLOW_MANAGER}`
    rescue
    end
    unless hello =~ /hello/
      puts "\e[31mFAILURE\e[0m: workflow_manager does not reply. check if workflow_manager is working"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end

    if failures > 0
      puts
      puts "\e[31mFailures (#{failures})\e[0m: All failures should be solved"
      raise "test run fails"
    else
      puts "All checks \e[32mPASSED\e[0m"
    end
  end
end



