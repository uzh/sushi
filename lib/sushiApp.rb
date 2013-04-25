#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130425-085436'

require 'csv'
require 'fileutils'
require 'active_record'

SUSHI_APP_DIR='/srv/GT/analysis/masaomi/sushi/test_new_data_set'
SUSHI_DB_TYPE='sqlite3'
WORKSPACE_DIR='/srv/GT/analysis/masaomi/sushi/sushi_launch_box/gstore/sushi'
#WORKSPACE_DIR='/srv/gstore/projects'

class SushiApp
  class Param
    attr_accessor :value
    attr_reader :type
    def initialize(type, default)
      @value = default
      @type = type
    end
  end
  attr_reader :params
  attr_accessor :dataset_tsv_file
  attr_accessor :parameterset_tsv_file
  attr_accessor :dataset_sushi_id
  def initialize
    @workspace = WORKSPACE_DIR ## will be on g-store; writable by sushi
    @project = "p999"
    @params = {}
    @params['process_mode'] = Param.new(String, 'SAMPLE')
  end
  def save_params_as_tsv
    file_name = 'parameters.tsv'
    file_path = File.join(@result_dir_extended, file_name)
    CSV.open(file_path, 'w', :col_sep=>"\t") do |out|
      out << @output_params.keys
      out << @output_params.values.map{|param| param.value}
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
    @output_params = @params.clone
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
          @params[header] ||= Param.new(String, row[header])
          @params[header].value = if @params[header].type == String
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
    @output_files.map{|header| output[header]}.each do |file|
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
  def output
    # this should be overwritten in a subclass
  end
  def commands
    # this should be overwritten in a subclass
  end
  def submit(job_script)
    job_id = `wfm_monitoring #{job_script}`
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
  def run
    set_input_dataset
    exit unless check_required_columns
    check_application_paramters
    set_user_parameters

    ## the user presses RUN
    prepare_result_dir

    ## sushi writes creates the job scripts and builds the result data set that is to be generated
    @result_dataset = []
    @job_scripts = []
    if @params['process_mode'].value == 'SAMPLE'
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
        @result_dataset << output
      end
      @job_scripts.each do |job_script|
        print "Submit job #{File.basename(job_script)} job_id=#{submit(job_script)}"
      end
    elsif @params['process_mode'].value == 'DATASET'
      # This should be implemented in a subclass (Application class)
    else 
      #stop
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
end



