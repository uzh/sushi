# Rails-free: ruby test/unit/methods_app_command_test.rb (from master/)
# The Methods job's R snippet must fall back to the base EzApp writer when the
# app's class is not in ezRun (EzAppNfCoreGeneric lives in a sourced file), and
# run_PyApp apps must get a Methods job written by EzApp.
require 'minitest/autorun'

module SushiFabric
  class SushiApp
    def initialize; @params = {}; end
  end
end
lib = File.expand_path('../../lib', __dir__)
load File.join(lib, 'MethodsApp.rb')
load File.join(lib, 'global_variables.rb')

class FakePyApp
  include GlobalVariables
  def initialize
    @params = { 'process_mode' => 'DATASET' }
    @gstore_dir = '/srv/gstore/projects'
    @result_dir = 'p1/x'
    @input_dataset_tsv_path = 'input_dataset.tsv'
  end
  def next_dataset; { 'Name' => 'x' }; end
  attr_reader :ezrun_class_name
end

class MethodsAppCommandTest < Minitest::Test
  def methods_command(cls)
    MethodsApp.new(ezrun_class_name: cls, analysis_name: 'A', next_dataset_id: 1,
                   gstore_result_dir: 'g', scratch_result_dir: 's', job_script_dir: 'j',
                   gstore_script_dir: '/gs/scripts', sushi_server: nil,
                   sample_count: 3, example_script: 'a.sh').commands
  end

  def assert_writer(cmd, cls)
    assert_includes cmd, %Q{cls <- if (exists("#{cls}")) get("#{cls}") else EzApp\n}
    assert_includes cmd, "cls\\$new()\\$write_methods(\n"
    assert_includes cmd, "  gstore_script_dir = '/gs/scripts',\n"
    assert_includes cmd, "  sample_count      = 3\n"
    refute_match(/^#{cls}\\\$new/, cmd)
  end

  def test_r_app
    assert_writer(methods_command('EzAppSTAR'), 'EzAppSTAR')
  end

  def test_class_not_in_ezrun
    assert_writer(methods_command('EzAppNfCoreGeneric'), 'EzAppNfCoreGeneric')
  end

  def test_py_app_gets_base_writer
    app = FakePyApp.new
    app.run_PyApp('Bin2cell', conda_env: 'gi_bin2cell')
    assert_equal 'EzApp', app.ezrun_class_name
    assert_writer(methods_command(app.ezrun_class_name), 'EzApp')
  end
end
