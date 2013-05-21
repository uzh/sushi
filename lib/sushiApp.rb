#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130521-170105'

require 'csv'
require 'fileutils'
require 'active_record'

SUSHI_APP_DIR='/srv/GT/analysis/masaomi/sushi/work_party'
SUSHI_DB_TYPE='sqlite3'
WORKSPACE_DIR='/srv/GT/analysis/masaomi/sushi/work_lunch/gstore/sushi'
#WORKSPACE_DIR='/srv/gstore/projects'


class Hash
  attr_reader :defaults
  alias :set :[]=
  def []=(k,v)
    @defaults ||= {}
    unless @defaults[k]
      @defaults.set(k,v)
    end
    set(k,v)
  end
  def default_value(k)
    @defaults[k]
  end
  def data_type(k)
    @defaults[k].class
  end
  def data_types
    Hash[@defaults.map{|k,v| [k, v.class]}]
  end
end
class SushiApp
  attr_reader :params
  attr_reader :job_ids
  attr_accessor :dataset_tsv_file
  attr_accessor :parameterset_tsv_file
  attr_accessor :dataset_sushi_id
  attr_accessor :project
  def initialize
    @workspace = WORKSPACE_DIR ## will be on g-store; writable by sushi
    @project = ''
    @params = {}
    @params['cores'] = 1
    @params['ram'] = 1
    @params['scratch'] = 1
    @params['node'] = ''
    @params['process_mode'] = 'SAMPLE'
    @job_ids = []
  end
  def save_params_as_tsv
    file_name = 'parameters.tsv'
    file_path = File.join(@result_dir_extended, file_name)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << @output_params.keys
      out << @output_params.values
    end
  end
  def set_input_dataset
    if @dataset_tsv_file
      dataset_tsv = CSV.readlines(@dataset_tsv_file, {:headers=>true, :col_sep=>"\t"})
      @dataset_hash = []
      dataset_tsv.each do |row|
        @dataset_hash << row.to_hash
      end
    elsif @dataset_sushi_id
      require "#{SUSHI_APP_DIR}/app/models/data_set"
      require "#{SUSHI_APP_DIR}/app/models/sample"
      ActiveRecord::Base.establish_connection(
            :adapter  => SUSHI_DB_TYPE,
            :database => "#{SUSHI_APP_DIR}/db/development.#{SUSHI_DB_TYPE}"
        )
      @dataset_hash = []
      if dataset = DataSet.find_by_id(@dataset_sushi_id.to_i)
        dataset.samples.each do |sample|
          @dataset_hash << sample.to_hash
        end
      end
    end
  end
  def check_required_columns
    if (@required_columns-@dataset_hash.map{|row| row.keys}.flatten.uniq).empty?
      true
    else
      false
    end
  end
  def check_application_paramters
    if (@required_params - @params.keys).empty?
      @output_params = @params.clone
    end
  end
  def output_parameters
    @output_params.keys.each do |key|
      @output_params[key] = @params[key]
    end
  end
  def set_user_parameters
    # this should be done in an instance of applicaiton subclass
    if @parameterset_tsv_file
      parameterset_tsv = CSV.readlines(@parameterset_tsv_file, :headers=>true, :col_sep=>"\t")
      parameterset_tsv.each do |row|
        row.headers.each do |header|
          @params[header] ||= row[header]
          @params[header] = if @params.data_type(header) == String
                                    row[header]
                                  else
                                    eval(row[header])
                                  end
        end
      end
    end
  end
  def prepare_result_dir
    ## sushi figures out where to put the resulting dataset
    @result_dir = File.join(@project, [@analysis_category, @name, Time.now.strftime("%Y-%m-%d--%H-%M-%S")].join("_"))
    @result_dir_extended = File.join(@workspace, @result_dir)
    FileUtils.mkdir_p(@result_dir_extended)
  end
  def job_header
    @out.print <<-EOF
#!/bin/sh

#### SET THE STAGE
#{File.join("SCRATCH_DIR=/scratch/", [File.basename(@result_dir), @dataset['Sample']].join("_"))}
WORKSPACE_DIR=#{@workspace}
mkdir $SCRATCH_DIR
cd $SCRATCH_DIR

    EOF
  end
  def job_footer
    @out.print "#### JOB IS DONE WE PUT THINGS IN PLACE AND CLEAN AUP\n"
    @output_files.map{|header| next_dataset[header]}.each do |file|
      # in actual case, to save under /srv/gstore/
      if WORKSPACE_DIR =~ /srv\/gstore/
        @out.print "gstore-request copy ", File.basename(file), " ", File.join(@workspace, file), "\n"
      else
        @out.print "cp ", File.basename(file), " ", File.join(@workspace, file), "\n"
      end
    end
    @out.print <<-EOF
cd ~
rm -rf $SCRATCH_DIR
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
  def submit(job_script)
    gsub_options = []
    gsub_options << "-c #{@params['cores']}" unless @params['cores'].to_s.empty?
    gsub_options << "-n #{@params['node']}" unless @params['node'].to_s.empty?
    gsub_options << "-r #{@params['ram']}" unless @params['ram'].to_s.empty?
    gsub_options << "-s #{@params['scratch']}" unless @params['scratch'].to_s.empty?
    job_id = `wfm_monitoring #{job_script} #{gsub_options.join(' ')}`
    job_id = job_id.to_i
    unless job_id.to_i > 1
      raise 'failed in job submitting'
    end
    job_id
  end
  def save_input_dataset_as_tsv
    file_name = 'input_dataset.tsv'
    file_path = File.join(@result_dir_extended, file_name)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      headers = @dataset_hash.map{|row| row.keys}.flatten.uniq
      out << headers
      @dataset_hash.each do |row|
        out << headers.map{|header| row[header]}
      end
    end
  end
  def save_next_dataset_as_tsv
    headers = @result_dataset.map{|row| row.keys}.flatten.uniq
    file_name = 'next_dataset.tsv'
    file_path = File.join(@result_dir_extended, file_name)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << headers
      @result_dataset.each do |row_hash|
        out << headers.map{|header| row_hash[header]}
      end
    end
  end
  def preprocess
    # this should be overwritten in a subclass
  end
  def run
    test_run

    ## the user presses RUN
    prepare_result_dir

    ## sushi writes creates the job scripts and builds the result data set that is to be generated
    @result_dataset = []
    @job_scripts = []
    if @params['process_mode'] == 'SAMPLE'
      @dataset_hash.each do |row|
        @dataset = row
        ## WRITE THE JOB SCRIPT
        @job_script = File.join(@result_dir_extended, row['Sample']) + '.sh'
        @out = open(@job_script, 'w')
        job_header
        job_main
        job_footer
        @out.close
        @job_scripts << @job_script
        @result_dataset << next_dataset
      end
      @job_scripts.each do |job_script|
        job_id = submit(job_script)
        @job_ids << job_id
        print "Submit job #{File.basename(job_script)} job_id=#{job_id}"
        open(File.join(@result_dir_extended, 'get_log.sh'), 'w') do |out|
          out.print "#!/bin/sh\n\n"
          out.print 'CURDIR=`dirname $0`', "\n"
          out.print "wfm_get_log #{job_id} with_err > $CURDIR/log.dat\n"
        end
      end
    elsif @params['process_mode'] == 'DATASET'
      # This should be implemented in a subclass (Application class)
    else 
      #stop
      warn "the process mode (#{@params['process_mode']}) is not defined"
      raise "stop job submitting"
    end

		puts
    print 'job scripts: '
    p @job_scripts
    print 'result dataset: '
    p @result_dataset

    save_params_as_tsv
    save_input_dataset_as_tsv
    save_next_dataset_as_tsv
  end
  def test_run
    preprocess
    set_input_dataset
    set_user_parameters

    failures = 0
    print 'check project name: '
    if @project.empty?
      puts "\e[31mFAILURE\e[0m: project number is required but not found. you should set it in usecase."
      puts "\tex.)"
      puts "\tapp = #{self.class}.new"
      puts "\tapp.project = 'p1001'"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@project=#{@project}"
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
      ptus "\tusecase = #{self.class}.new"
      puts "\tusecase.dataset_tsv_file = \"dataset.tsv\""
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:\n\t@dataset_hash.length = #{@dataset_hash.length}"
    end

    print 'check required columns: '
    unless check_required_columns
      puts "\e[31mFAILURE\e[0m: required_column(s) is not found in dataset"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end
    puts "\trequired columns: #{@required_columns}"
    puts "\tdataset  columns: #{@dataset_hash.map{|row| row.keys}.flatten.uniq}"

    print 'check required parameters: '
    unless check_application_paramters
      puts "\e[31mFAILURE\e[0m: required_param(s) is not set yet. you should set it in usecase"
      puts "\tmissing params: #{@required_params-@params.keys}"
      puts "\tex.)"
      puts "\tusecase = #{self.class}.new"
      puts "\tusecase.params['#{(@required_params-@params.keys)[0]}'] = parameter"
      puts
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end
    puts "\tparameters: #{@params.keys}"
    puts "\trequired  : #{@required_params}"

    print 'check next dataset: '
    @dataset={}
    @result_dir='.'
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
      puts "\t      #{self.class}#next_dataset.keys: #{self.next_dataset.keys}"
    else
      puts "\e[32mPASSED\e[0m:"
    end

    print 'check commands: '
    unless commands
      puts "\e[31mFAILURE\e[0m: any commands is not defined yet. you should overwrite SushiApp#commands method in #{self.class}"
      puts "\tnote: the return value should be String (this will be in the main body of submitted job script)"
      failures += 1
    else
      puts "\e[32mPASSED\e[0m:"
    end

    print 'check workflow manager: '
    begin
      hello = `wfm_hello`
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



