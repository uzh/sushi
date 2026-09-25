# DataSet.pick_own_result_dir: which of a dataset's File/Link dirs is its OWN result dir.
# Pure logic, so it runs without booting Rails: `rspec spec/models/data_set_own_result_dir_spec.rb`
unless defined?(ActiveRecord::Base)
  module ActiveRecord; class Base; def self.method_missing(*); end; end; end
end
require_relative '../../app/models/data_set'

describe 'DataSet.pick_own_result_dir' do
  def pick(name, dirs)
    DataSet.pick_own_result_dir(name, dirs)
  end

  it 'picks the ScSeurat dir, not the CellBender dir its CountMatrix [Link] points into' do
    dirs = ['p28409/HD_blood_CellBender_2026-08-31--08-36-55',
            'p28409/HD_blood_ScSeurat_2026-09-01--14-08-37']
    expect(pick('HD_blood_ScSeurat', dirs)).to eq 'p28409/HD_blood_ScSeurat_2026-09-01--14-08-37'
  end

  it 'picks the dir matching the name even when it is not the newest' do
    dirs = ['p1001/rawMAD_input_2026-05-01--10-00-00',
            'p1001/rawMAD_2026-04-01--10-00-00',
            'p1001/other_2026-06-01--10-00-00']
    expect(pick('rawMAD', dirs)).to eq 'p1001/rawMAD_2026-04-01--10-00-00'
  end

  it 'falls back to the newest timestamp when no dir matches the name' do
    dirs = ['p1001/CellRangerMulti_2026-08-31--08-36-55',
            'p1001/ScSeurat_4242_2026-09-01--14-08-37']
    expect(pick('ScSeurat_4243', dirs)).to eq 'p1001/ScSeurat_4242_2026-09-01--14-08-37'
  end

  it 'falls back to the first dir when none carries a timestamp' do
    expect(pick('x', ['p1001/data', 'p1001/other'])).to eq 'p1001/data'
  end

  it 'returns nil for no dirs' do
    expect(pick('x', [])).to be_nil
  end
end
