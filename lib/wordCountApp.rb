#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130419-101514'

require './sushiApp'

class WordCountApp < SushiApp
  def initialize
    super
    @name = 'Word_Count'
    @required_columns = ['Read1','Species']
    @params['lines_only'] = SushiApp::Param.new('boolean', true)
    @analysis_category = 'Stats'
    @output_files = ['Stats']
  end
  def output
    {'Sample'=>@dataset['Sample'], 'Stats'=>File.join(@result_dir, @dataset['Sample']+'.stats')}
  end
  def commands
    "gunzip -c $WORKSPACE_DIR/#{@dataset['Read1']} |wc > #{File.basename(output['Stats'])}"
  end
end

