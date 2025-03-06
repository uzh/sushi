#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric/sushiApp'

include SushiFabric
describe "Array monkey patch" do
  describe "[1,'a',2,'b'].to_h" do
    subject(:array_to_h){[1,"a",2,"b"].to_h}
    let(:hash){{1=>"a", 2=>"b"}}
    it {is_expected.to eq hash}
  end
end

describe "Hash monkey patch" do
  subject(:hash){{}}
  context "hash.set(1, 'a')" do
    describe "hash.get(1)" do
      before{hash.set(1, 'a')}
      specify{expect(hash.get(1)).to eq 'a'}
    end
  end
  context "hash[1, 'desc']='test'" do
    describe "hash[1, 'desc']" do
      before{hash[1, 'desc'] = 'test'}
      subject(:target){hash[1, 'desc']}
      it{is_expected.to eq 'test'}
    end
  end
  describe "#data_type" do
    context "hash[:key]=1" do
      before{hash[:key]=1}
      describe "hash.data_type(:key)" do
        subject{hash.data_type(:key)}
        it{is_expected.to eq Fixnum}
      end
    end
    context "hash[:key]=['a',1,2,3]" do
      before{hash[:key]=["a",1,2,3]}
      describe "hash.data_type(:key)" do
        subject{hash.data_type(:key)}
        it{is_expected.to eq String}
      end
    end
  end
  describe "#default_value" do
    context "hash.default_value(:key, :default_value)" do
      before do
        hash[:key] = 0
        hash.instance_variable_set(:@defaults, {})
        hash.default_value(:key, :default_value)
      end
      describe "hash.default_value(:key)" do
        subject{hash.default_value(:key)}
        it{is_expected.to eq :default_value}
      end
    end
  end
end

describe "String monkey patch" do
  describe "#tag?" do
    context "string = 'Read1 [File]'" do
      let(:string){"Read1 [File]"}
      describe "string.tag?('File')" do
        subject{string.tag?("File")}
        it{is_expected.to_not be_falsey}
      end
      describe "string.tag?('Hoge')" do
        subject{string.tag?("Hoge")}
        it{is_expected.to be_falsey}
      end
    end
  end
end

describe SushiApp do
  subject(:sushi_app) {SushiApp.new}
  context 'when new' do
    it {is_expected.to be_an_instance_of SushiApp} 
  end 
  describe "#job_header" do
    subject{sushi_app.job_header}
    let(:dataset) {{'Name' => 'Name'}}
    let(:out) {double('out')}
    before do
      allow(out).to receive_messages(:print => nil)
      sushi_app.instance_variable_set(:@out, out)
      sushi_app.instance_variable_set(:@scratch_result_dir, 'scratch_result_dir')
      sushi_app.instance_variable_set(:@dataset, dataset)
    end
    it {is_expected.to be_nil}
  end 
end
