#!/usr/bin/env ruby
# encoding: utf-8

require 'spec_helper'
require 'csv'

describe SushiApp do
  let(:sushi) {SushiApp.new}
  context 'when new' do
    it {should be_an_instance_of SushiApp}
    its(:params) {should be_an_instance_of Hash}
    it "params.keys.sort should have some elements" do 
      subject.params.keys.sort.should == ["cores", "node", "process_mode", "ram", "scratch"]
    end
    its(:job_ids) {should be_an_instance_of Array}
    its(:dataset_tsv_file) {should be_nil} 
    its(:parameterset_tsv_file) {should be_nil}
    its(:dataset_sushi_id) {should be_nil}
    its(:project) {should be_nil}
    its(:user) {should be_nil}
    its(:next_dataset) {should be_nil}
    its(:commands) {should be_nil}
  end
  describe '#copy_commands' do
    subject {sushi.copy_commands('file1','path')}
    context 'when @gstore includes "gstore"' do
      before do
        sushi.instance_variable_set(:@gstore_dir, '/srv/gstore')
      end
      it {should eq ["g-req -w copy file1 path"]}
    end
    before do
      sushi.instance_variable_set(:@gstore_dir, '/srv/dummy')
    end
    context 'when @gstore does not include "gstore"' do
      it {should eq ["mkdir -p path","cp -r file1 path"]}
    end
  end
  describe '#set_input_dataset' do
    let(:result) do
      [{"Sample"=>"mysample1",
        "Read1"=>"p1001/data/short-ama_E1_R1.fastq.gz",
        "Read2"=>"p1001/data/short-ama_E1_R2.fastq.gz",
        "Species"=>"Arabidopsis"}]
    end
    before do 
      @sushi = SushiApp.new
      headers = ['Sample', 'Read1', 'Read2', 'Species']
      fields = ['mysample1', 'p1001/data/short-ama_E1_R1.fastq.gz', 'p1001/data/short-ama_E1_R2.fastq.gz', 'Arabidopsis']
      row = CSV::Row.new(headers, fields)
      @csv = [row]
    end
    subject {@sushi.set_input_dataset}
    context 'when @dataset_tsv_file is used' do
      before do
        @sushi.dataset_tsv_file = 'tsv'
        CSV.should_receive(:readlines).and_return(@csv)
      end
      it {should eq result}
    end
    context 'when @dataset_sushi_id is used' do
      before do
        @sushi.dataset_sushi_id = 1
        dataset = double :dataset, :samples=>@csv
        DataSet.should_receive(:find_by_id).and_return(dataset)
      end
      it {should eq result}
    end
  end
  describe '#save_parameters_as_tsv' do
    let(:out) {[]}
    before do
      sushi.instance_variable_set(:@scratch_result_dir, '/scratch_result_dir')
      sushi.instance_variable_set(:@parameter_file, 'parameters.tsv')
      output_params = {'cores' => '4', 'ram' => '16'}
      output_params.each do |key, value|
        out << [key, value]
      end
      CSV.stub(:open).and_yield(out)
      sushi.instance_variable_set(:@output_params, output_params)
    end
    context 'return value' do
      subject {sushi.save_parameters_as_tsv}
      it {should eq "/scratch_result_dir/parameters.tsv"}
    end
    context 'output tsv file' do
      subject {out}
      it {should eq [['cores','4'], ['ram','16']]}
    end
  end
  describe '#set_user_parameters' do
    before do
      parameterset_tsv_file = [['cores','4'], ['ram','16'], ['process_mode','SAMPLE']]
      CSV.stub(:readlines).and_return(parameterset_tsv_file)
    end
    context 'when @parameterset_tsv_file is set' do
      before do
        sushi.parameterset_tsv_file = true
      end
      let(:result) {{"cores"=>4, "ram"=>16, "scratch"=>nil, "node"=>"", "process_mode"=>"SAMPLE"}}
      subject {sushi.set_user_parameters}
      it {should eq result}
    end
    context 'when @parameterset_tsv_file is not set' do
      before do
        sushi.parameterset_tsv_file = false
      end
      let(:result) {{"cores"=>nil, "ram"=>nil, "scratch"=>nil, "node"=>"", "process_mode"=>"SAMPLE"}}
      subject {sushi.set_user_parameters}
      it {should eq result}
    end
  end
  describe '#submit_command' do
    before do
      sushi.params['cores'] = '4'
      sushi.params['ram']   = '7'
      sushi.params['scratch'] = '16'
      sushi.params['node']  = 'fgcz-c-046'
      sushi.params['user']  = 'sushi'
      sushi.project = 'p1001'
      sushi.instance_variable_set(:@gstore_script_dir, '/gstore_script_dir')
    end
    let(:result) {"wfm_monitoring --server #{WORKFLOW_MANAGER} --project 1001 --logdir /gstore_script_dir job_script -c 4 -n fgcz-c-046 -r 8 -s 16"}
    subject {sushi.submit_command('job_script')}
    it {should eq result}
  end
  describe '#prepare_result_dir' do
    before do
      FileUtils.stub(:mkdir_p).twice.and_return('mkdir_p done')
    end
    subject {sushi.prepare_result_dir}
    it {should eq 'mkdir_p done'}
  end
  describe '#set_file_paths' do
    before do
      sushi.instance_variable_set(:@gstore_result_dir, '/gstore_result_dir')
    end
    subject {sushi.set_file_paths}
    it {should eq "/gstore_result_dir/dataset.tsv"}
  end
  describe '#set_dir_paths' do
    context 'when @name and @project are not set' do
      it {lambda{sushi.set_dir_paths}.should raise_error(RuntimeError)}
    end
    context 'when @name and @project are set' do
      before do
        sushi.instance_variable_set(:@name, 'name')
        sushi.project = 'p1001'
        sushi.should_receive(:set_file_paths).and_return('set_file_paths')
      end
      subject {sushi.set_dir_paths}
      it {should eq 'set_file_paths'}
    end
  end
  describe '#set_output_files' do
    before do
      next_dataset = {
        'Name' => 'name',
        'Bam [File]' => 'test.bam'
      }
      sushi.stub(:next_dataset).and_return(next_dataset)
    end
    subject {sushi.set_output_files}
    it {should eq ["Bam [File]"]}
  end
  describe '#check_required_columns' do
    let(:required_columns) {['Name', 'Species']}
    before do
      sushi.instance_variable_set(:@required_columns, required_columns)
    end
    context 'when required_columns is not satified' do
      let(:dataset_hash) {
        [{'Name' => 'name1'},
         {'Name' => 'name2'}]
      }
      before do
        sushi.instance_variable_set(:@dataset_hash, dataset_hash)
      end
      subject {sushi.check_required_columns}
      it {should be_false}
    end
    context 'when required_columns is statisfied' do
      let(:dataset_hash) {
        [{'Name' => 'name1', 'Species [File]' => 'A.thaliana'},
         {'Name' => 'name2', 'Species [File]' => 'A.thaliana'}]
      }
      before do
        sushi.instance_variable_set(:@dataset_hash, dataset_hash)
      end
      subject {sushi.check_required_columns}
      it {should be_true}
    end
  end
  describe '#check_application_parameters' do
    let(:required_params) {['Name', 'Bam', 'Species']}
    before do
      sushi.instance_variable_set(:@required_params, required_params)
    end
    context 'when required_params is not satisfied' do
      before do
        sushi.params['Name'] = 'name1'
        sushi.params['Species'] = 'A.thaliana'
      end
      subject {sushi.check_application_parameters}
      it {should be_nil}
    end
    context 'when required_params is satisfied' do
      before do
        sushi.params['Name'] = 'name1'
        sushi.params['Bam']  = 'bam1'
        sushi.params['Species'] = 'A.thaliana'
      end
      let(:result) {
        {"cores"=>nil, "ram"=>nil, "scratch"=>nil, "node"=>"", "process_mode"=>"SAMPLE", "Name"=>"name1", "Bam"=>"bam1", "Species"=>"A.thaliana"}
      }
      subject {sushi.check_application_parameters}
      it {should eq result}
    end
  end
end
