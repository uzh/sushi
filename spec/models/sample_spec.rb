require 'spec_helper'

describe Sample do
  let(:sample) {Sample.new}
  let(:sample_hash) {
    {
        'Name'  => 'sample1',
        'Read1' => 'p1001/data/Cama_R1.fastq.gz'
    }
  }
  describe "#to_hash" do
    before do
      sample.key_value = sample_hash.to_s
      sample.save
    end
    subject {sample.to_hash}
    it {should eq sample_hash}
  end
  describe "#saved?" do
    context "before save" do
      subject {sample.saved?}
      it {should be_false}
    end
    context "after save" do
      before do
        sample.save
      end
      subject {sample.saved?}
      it {should be_true}
    end
  end
end
