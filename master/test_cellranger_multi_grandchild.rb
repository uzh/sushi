#!/usr/bin/env ruby
# Test script for CellRangerMultiApp grandchild_datasets() method (v3)
# Tests the corrected approach: iterating over @result_dataset instead of stale @dataset

require 'csv'

# Mock SushiFabric module
module SushiFabric
  GSTORE_DIR = '/srv/gstore/projects'
end

# Mock the CellRangerMultiApp with relevant attributes
# Simulates the real framework behavior where @dataset, @result_dataset, and @dataset_hash
# are all available after sample_mode completes
class TestCellRangerMulti
  attr_accessor :result_dataset, :dataset_hash, :dataset

  def initialize(params, dataset, result_dir, dataset_hash: nil, result_dataset: nil)
    @params = params
    @dataset = dataset
    @result_dir = result_dir
    @inherit_tags = ["Factor", "B-Fabric"]
    # @dataset_hash: original input rows WITH tagged keys (as in input_dataset.tsv)
    @dataset_hash = dataset_hash || []
    # @result_dataset: accumulated next_dataset results for selected samples
    @result_dataset = result_dataset || []
  end

  def extract_columns(tags)
    # Mock: return inherited columns from @dataset (stale in sample_mode!)
    result = {}
    if @dataset.is_a?(Hash)
      @dataset.each do |key, value|
        if key.include?('[Factor]') || key.include?('[B-Fabric]')
          result[key] = value
        end
      end
    end
    result
  end

  def grandchild_datasets
    # Only create grandchild datasets when multiplexing is used
    return [] unless @params['TenXLibrary'].to_s.include?('Multiplexing')

    # Detect CellRanger version for output path structure
    cr_match = @params['CellRangerVersion'].to_s.match(/(\d+)\./)
    is_v9_or_below = cr_match ? cr_match[1].to_i <= 9 : false

    # In sample_mode, @dataset points to the LAST input row after the loop,
    # not necessarily the processed sample. Use @result_dataset (populated
    # during the loop for selected samples only) to get correct sample names,
    # then look up input metadata from @dataset_hash.
    processed_names = if @result_dataset.is_a?(Array) && !@result_dataset.empty?
                        @result_dataset.map { |r| r['Name'] }.compact
                      else
                        [@dataset['Name']].compact
                      end

    grandchild_data = []
    processed_names.each do |sample_name|
      # Find the original input row in @dataset_hash (with tagged keys)
      input_row = if @dataset_hash.is_a?(Array)
                    @dataset_hash.find { |row|
                      stripped_name = row.find { |k, _| k.gsub(/\[.+\]/,'').strip == 'Name' }
                      stripped_name && stripped_name[1] == sample_name
                    }
                  end
      # Strip tags to get clean keys
      sample_data = if input_row
                      Hash[*input_row.map{|key,value| [key.gsub(/\[.+\]/,'').strip, value]}.flatten]
                    else
                      @dataset
                    end

      # Locate Sample2Barcode file in gStore metadata directory
      order_id = sample_data['Order Id']
      raw_data_dir = sample_data['RawDataDir']
      project_id = raw_data_dir.to_s.split('/').first

      next if project_id.to_s.empty? || order_id.to_s.empty?

      metadata_dir = File.join(SushiFabric::GSTORE_DIR, project_id, "o#{order_id}_metaData")
      next unless Dir.exist?(metadata_dir)

      # Find matching Sample2Barcode file for this sample
      sample2barcode_files = Dir.glob(File.join(metadata_dir, "*_Sample2Barcode.csv"))
      matching_file = sample2barcode_files.find { |f|
        prefix = File.basename(f).sub(/_Sample2Barcode\.csv$/i, '')
        sample_name.start_with?(prefix) || prefix == sample_name
      }
      next unless matching_file

      # Parse Sample2Barcode CSV to get demultiplexed sample names
      sample2barcode = CSV.read(matching_file, headers: true)
      sub_sample_names = sample2barcode['sample_id'].compact.reject(&:empty?)
      next if sub_sample_names.empty?

      # Construct grandchild datasets with predicted output paths
      sub_sample_names.each do |sub_sample|
        per_sample_dir = File.join(@result_dir, sample_name, 'per_sample_outs', sub_sample)

        dataset = {
          'Name' => sub_sample,
          'Species' => sample_data['Species'],
          'refBuild' => @params['refBuild'],
          'refFeatureFile' => @params['refFeatureFile'],
          'featureLevel' => @params['featureLevel'],
          'transcriptTypes' => @params['transcriptTypes'],
          'SCDataOrigin' => '10X',
          'Report [Link]' => File.join(per_sample_dir, 'web_summary.html')
        }

        if is_v9_or_below
          dataset['CountMatrix [Link]'] = File.join(per_sample_dir, 'count', 'sample_filtered_feature_bc_matrix')
          dataset['UnfilteredCountMatrix [Link]'] = File.join(per_sample_dir, 'count', 'sample_raw_feature_bc_matrix')
        else
          dataset['CountMatrix [Link]'] = File.join(per_sample_dir, 'sample_filtered_feature_bc_matrix')
          dataset['UnfilteredCountMatrix [Link]'] = File.join(per_sample_dir, 'sample_raw_feature_bc_matrix')
        end

        # Extract inherited columns directly from input row (not from stale @dataset)
        if input_row
          @inherit_tags.each do |tag|
            input_row.each do |key, value|
              dataset[key] = value if key.include?("[#{tag}]")
            end
          end
        else
          dataset.merge!(extract_columns(@inherit_tags))
        end

        grandchild_data << dataset
      end
    end

    grandchild_data
  end
end

# Test setup: simulate the real SUSHI framework scenario
# Input dataset has 3 samples (with tagged keys, as in input_dataset.tsv)
dataset_hash = [
  {
    'Name' => 'SKNBE2_undifferentiated_3',
    'Condition [Factor]' => 'NA',
    'RawDataDir [File]' => 'p38593/o39496_NovaSeq_260108_X416/SKNBE2_undifferentiated_3.tar',
    'Species' => 'Homo sapiens',
    'Sample Id [B-Fabric]' => '1081947',
    'Order Id [B-Fabric]' => '39496'
  },
  {
    'Name' => 'SKNSH_differentiated_2',
    'Condition [Factor]' => 'NA',
    'RawDataDir [File]' => 'p38593/o39496_NovaSeq_260108_X416/SKNSH_differentiated_2.tar',
    'Species' => 'Homo sapiens',
    'Sample Id [B-Fabric]' => '1081952',
    'Order Id [B-Fabric]' => '39496'
  },
  {
    'Name' => 'SKNSH_undifferentiated_1',
    'Condition [Factor]' => 'NA',
    'RawDataDir [File]' => 'p38593/o39496_NovaSeq_260108_X416/SKNSH_undifferentiated_1.tar',
    'Species' => 'Homo sapiens',
    'Sample Id [B-Fabric]' => '1081972',
    'Order Id [B-Fabric]' => '39496'
  }
]

params = {
  'TenXLibrary' => 'GEX,Multiplexing',
  'CellRangerVersion' => 'Aligner/CellRanger/10.0.0',
  'refBuild' => 'Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03',
  'refFeatureFile' => 'genes.gtf',
  'featureLevel' => 'gene',
  'transcriptTypes' => 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
}
result_dir = 'p38593/o39496_CellRangerMulti_2026-02-09--16-35-26'

# Test cases
puts "=" * 80
puts "Testing CellRangerMultiApp grandchild_datasets() v3"
puts "Fix: iterate over @result_dataset, not stale @dataset"
puts "=" * 80

# Test Case 1: CRITICAL BUG TEST - stale @dataset scenario
# In sample_mode, @dataset ends up as the LAST input row (SKNSH_undifferentiated_1)
# but only SKNBE2_undifferentiated_3 was selected for processing.
# The old code would use SKNSH_undifferentiated_1's Sample2Barcode file → WRONG paths.
puts "\n[Test 1] Stale @dataset - only SKNBE2_undifferentiated_3 processed"
puts "-" * 80

# @dataset is stale (last row = SKNSH_undifferentiated_1)
stale_dataset = {
  'Name' => 'SKNSH_undifferentiated_1',
  'Species' => 'Homo sapiens',
  'RawDataDir' => 'p38593/o39496_NovaSeq_260108_X416/SKNSH_undifferentiated_1.tar',
  'Order Id' => '39496',
  'Condition [Factor]' => 'NA',
  'Sample Id [B-Fabric]' => '1081972'
}

# @result_dataset has only the processed sample
result_dataset = [
  { 'Name' => 'SKNBE2_undifferentiated_3', 'Species' => 'Homo sapiens' }
]

app = TestCellRangerMulti.new(params, stale_dataset, result_dir,
                              dataset_hash: dataset_hash,
                              result_dataset: result_dataset)
children = app.grandchild_datasets

puts "Number of grandchild datasets: #{children.length}"
test1_count = children.length == 4

# Verify correct sample names (_3 suffix, NOT _1)
expected_names = ['SKNBE2_undifferentiated_3', 'SKNBE2_differentiated_3', 'SKNSH_undifferentiated_3', 'SKNSH_differentiated_3']
actual_names = children.map { |c| c['Name'] }.sort
test1_names = expected_names.sort == actual_names
if test1_names
  puts "PASS: Sample names have correct _3 suffix (not stale _1)"
else
  puts "FAIL: Sample names wrong"
  puts "  Expected: #{expected_names.sort}"
  puts "  Actual: #{actual_names}"
end

# Verify paths use correct parent directory
test1_paths = children.all? { |c|
  c['Report [Link]'].include?('SKNBE2_undifferentiated_3/per_sample_outs') &&
  !c['Report [Link]'].include?('SKNSH_undifferentiated_1/per_sample_outs')
}
if test1_paths
  puts "PASS: Paths use correct parent (SKNBE2_undifferentiated_3)"
else
  puts "FAIL: Paths use wrong parent"
  children.each { |c| puts "  #{c['Report [Link]']}" }
end

# Verify inherited columns come from correct sample (1081947, not stale 1081972)
test1_inherit = children.all? { |c|
  c['Sample Id [B-Fabric]'] == '1081947'
}
if test1_inherit
  puts "PASS: Inherited Sample Id is 1081947 (correct sample)"
else
  puts "FAIL: Inherited Sample Id wrong"
  children.each { |c| puts "  #{c['Name']}: Sample Id = #{c['Sample Id [B-Fabric]']}" }
end

puts "\nGrandchild details:"
children.each_with_index do |child, i|
  puts "\n  Child #{i+1}: #{child['Name']}"
  puts "    Report: #{child['Report [Link]']}"
  puts "    CountMatrix: #{child['CountMatrix [Link]']}"
  puts "    Sample Id: #{child['Sample Id [B-Fabric]']}"
  puts "    Order Id: #{child['Order Id [B-Fabric]']}"
end

# Test Case 2: Without Multiplexing
puts "\n" + "=" * 80
puts "[Test 2] WITHOUT Multiplexing - should return empty array"
puts "-" * 80
params_no_mult = params.merge('TenXLibrary' => 'GEX,VDJ-T')
app2 = TestCellRangerMulti.new(params_no_mult, stale_dataset, result_dir,
                               dataset_hash: dataset_hash,
                               result_dataset: result_dataset)
children2 = app2.grandchild_datasets
test2 = children2.length == 0
puts test2 ? "PASS: No grandchildren (correct)" : "FAIL: Expected 0, got #{children2.length}"

# Test Case 3: Missing metadata directory
puts "\n" + "=" * 80
puts "[Test 3] Missing metadata directory - should return empty array"
puts "-" * 80
bad_dataset_hash = dataset_hash.map { |row|
  row.merge('Order Id [B-Fabric]' => '99999')
}
bad_result_dataset = [{ 'Name' => 'SKNBE2_undifferentiated_3', 'Species' => 'Homo sapiens' }]
app3 = TestCellRangerMulti.new(params, stale_dataset, result_dir,
                               dataset_hash: bad_dataset_hash,
                               result_dataset: bad_result_dataset)
children3 = app3.grandchild_datasets
test3 = children3.length == 0
puts test3 ? "PASS: Graceful fallback (correct)" : "FAIL: Expected 0, got #{children3.length}"

# Test Case 4: CellRanger v9 paths
puts "\n" + "=" * 80
puts "[Test 4] CellRanger v9 - should use /count/ subdirectory"
puts "-" * 80
params_v9 = params.merge('CellRangerVersion' => 'Aligner/CellRanger/9.0.0')
app4 = TestCellRangerMulti.new(params_v9, stale_dataset, result_dir,
                               dataset_hash: dataset_hash,
                               result_dataset: result_dataset)
children4 = app4.grandchild_datasets
test4 = children4.length > 0 && children4.first['CountMatrix [Link]'].include?('/count/')
if test4
  puts "PASS: v9 paths include /count/ subdirectory"
else
  puts "FAIL: v9 paths should include /count/"
  puts "  CountMatrix: #{children4.first['CountMatrix [Link]']}" if children4.length > 0
end

# Test Case 5: No ResultDir [File] tag
puts "\n" + "=" * 80
puts "[Test 5] No ResultDir [File] tag - should use [Link] only"
puts "-" * 80
has_file_tag = children.any? { |c| c.keys.any? { |k| k.include?('[File]') } }
test5 = !has_file_tag
puts test5 ? "PASS: No [File] tags" : "FAIL: Found [File] tags: #{children.first.keys.select{|k| k.include?('[File]')}}"

# Test Case 6: Column completeness
puts "\n" + "=" * 80
puts "[Test 6] Column Verification - all required fields present"
puts "-" * 80
required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'featureLevel',
                   'transcriptTypes', 'SCDataOrigin', 'Report [Link]',
                   'CountMatrix [Link]', 'UnfilteredCountMatrix [Link]']
all_valid = true
children.each do |child|
  required_columns.each do |col|
    unless child.key?(col) && !child[col].to_s.empty?
      puts "FAIL: Missing or empty '#{col}' in '#{child['Name']}'"
      all_valid = false
    end
  end
end
test6 = all_valid
puts "PASS: All required columns present and non-empty" if test6

# Test Case 7: Multiple samples processed (all 3 in @result_dataset)
puts "\n" + "=" * 80
puts "[Test 7] All 3 samples processed - should produce 12 grandchildren"
puts "-" * 80
result_dataset_all = [
  { 'Name' => 'SKNBE2_undifferentiated_3', 'Species' => 'Homo sapiens' },
  { 'Name' => 'SKNSH_differentiated_2', 'Species' => 'Homo sapiens' },
  { 'Name' => 'SKNSH_undifferentiated_1', 'Species' => 'Homo sapiens' }
]
app7 = TestCellRangerMulti.new(params, stale_dataset, result_dir,
                               dataset_hash: dataset_hash,
                               result_dataset: result_dataset_all)
children7 = app7.grandchild_datasets
puts "Number of grandchild datasets: #{children7.length}"
test7 = children7.length == 12

# Verify each parent's grandchildren have correct paths and names
if test7
  # Check SKNBE2_undifferentiated_3 grandchildren have _3 names
  gc3 = children7.select { |c| c['Report [Link]'].include?('SKNBE2_undifferentiated_3/per_sample_outs') }
  gc2 = children7.select { |c| c['Report [Link]'].include?('SKNSH_differentiated_2/per_sample_outs') }
  gc1 = children7.select { |c| c['Report [Link]'].include?('SKNSH_undifferentiated_1/per_sample_outs') }

  puts "PASS: 12 grandchildren total (#{gc3.length} from _3, #{gc2.length} from _2, #{gc1.length} from _1)"

  # Verify Sample Ids are correct for each group
  gc3_ids = gc3.map { |c| c['Sample Id [B-Fabric]'] }.uniq
  gc2_ids = gc2.map { |c| c['Sample Id [B-Fabric]'] }.uniq
  gc1_ids = gc1.map { |c| c['Sample Id [B-Fabric]'] }.uniq
  puts "  _3 group Sample Ids: #{gc3_ids} (expected: [1081947])"
  puts "  _2 group Sample Ids: #{gc2_ids} (expected: [1081952])"
  puts "  _1 group Sample Ids: #{gc1_ids} (expected: [1081972])"
else
  puts "FAIL: Expected 12, got #{children7.length}"
  children7.each { |c| puts "  #{c['Name']}: #{c['Report [Link]']}" }
end

# Test Case 8: Fallback with empty @result_dataset (uses @dataset)
puts "\n" + "=" * 80
puts "[Test 8] Fallback - empty @result_dataset uses @dataset"
puts "-" * 80
correct_dataset = {
  'Name' => 'SKNBE2_undifferentiated_3',
  'Species' => 'Homo sapiens',
  'RawDataDir' => 'p38593/o39496_NovaSeq_260108_X416/SKNBE2_undifferentiated_3.tar',
  'Order Id' => '39496'
}
app8 = TestCellRangerMulti.new(params, correct_dataset, result_dir,
                               dataset_hash: [],
                               result_dataset: [])
children8 = app8.grandchild_datasets
test8 = children8.length == 4
puts test8 ? "PASS: Fallback works (4 grandchildren)" : "FAIL: Expected 4, got #{children8.length}"

# Summary
puts "\n" + "=" * 80
puts "TEST SUMMARY"
puts "=" * 80
tests = [
  ["Test 1a (stale @dataset - count)", test1_count],
  ["Test 1b (stale @dataset - names)", test1_names],
  ["Test 1c (stale @dataset - paths)", test1_paths],
  ["Test 1d (stale @dataset - inherit)", test1_inherit],
  ["Test 2 (no multiplexing)", test2],
  ["Test 3 (missing metadata)", test3],
  ["Test 4 (v9 paths)", test4],
  ["Test 5 (no [File] tags)", test5],
  ["Test 6 (column completeness)", test6],
  ["Test 7 (all 3 samples → 12 grandchildren)", test7],
  ["Test 8 (fallback with empty result_dataset)", test8]
]
tests.each do |name, passed|
  puts "#{passed ? 'PASS' : 'FAIL'}: #{name}"
end
all_passed = tests.all? { |_, passed| passed }
puts "=" * 80
puts all_passed ? "ALL TESTS PASSED" : "SOME TESTS FAILED"
puts "=" * 80
