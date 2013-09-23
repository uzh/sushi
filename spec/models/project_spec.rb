require 'spec_helper'

describe Project do
  let(:project) {Project.new}
  let(:data_set){DataSet.new}
  before do
    project.data_sets << data_set
  end
  describe "#saved?" do
    context "when data_set.saved?" do
      before do
        data_set.stub(:saved?).and_return(true)
      end
      subject {project.saved?}
      it {should be_true}
    end
    context "when not data_set.saved?" do
      before do
        data_set.stub(:saved?).and_return(false)
      end
      subject {project.saved?}
      it {should be_false}
    end
  end
end
