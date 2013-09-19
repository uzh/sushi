require 'spec_helper'

describe Sample do
  describe "#to_hash" do
    let(:sample) {Sample.new}
    let(:sample_hash) {
      {
          'Name'  => 'sample1',
          'Read1' => 'p1001/data/Cama_R1.fastq.gz'
      }
    }
    before do
      sample.key_value = sample_hash.to_s
      sample.save
    end
    subject {sample.to_hash}
    it {should eq sample_hash}
  end
end
