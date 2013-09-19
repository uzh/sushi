require 'spec_helper'

describe DataSet do
  let(:data_set) {DataSet.new}
  let(:sample_hash) {
    {
        'Name'  => 'sample1',
        'Read1' => 'p1001/data/Cama_R1.fastq.gz'
    }
  }
  before do
    sample = Sample.new
    sample.key_value = sample_hash.to_s
    sample.save
    data_set.samples << sample 
  end
  describe "#headers" do
    let(:headers) {sample_hash.keys}
    subject {data_set.headers}
    it {should == headers}
  end
  describe "#saved?" do
    subject {data_set.saved?}
    context "before save" do
      it {should be_false}
    end
    context "after save" do
      before do
        data_set.md5 = data_set.md5hexdigest
        data_set.save
      end
      it {should be_true}
    end
  end
  describe "#export_tsv" do
    before do
      csv = double("csv")
    end
    pending
  end
end
