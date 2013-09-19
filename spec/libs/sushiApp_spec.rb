#!/usr/bin/env ruby
# encoding: utf-8

require 'spec_helper'
require 'csv'

describe SushiApp do
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
    let(:sushi) {SushiApp.new}
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

  describe '#prepare_result_dir' do
    it 'making @result_dir'
  end
end
