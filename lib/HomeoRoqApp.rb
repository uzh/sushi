#!/usr/bin/env ruby
# encoding: utf-8
Version = '20151126-160659'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class HomeoRoqApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'HomeRoq'
    @description =<<-EOS
HomeoRoq detects the significant genes that homeolog ratio changes in a target condition.<br />
Note<br />
<ol>
<li>The first character of column name in Count_QC-tpm.txt should not be a number. R data.frame puts X prefix if it is a number character.</li>
</ol>

http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Polyploid'
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Report']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '10'
    @params['name'] = 'HomeoRoq'
    @params['iteration'] = '10'
    @params['tpm_txt'] = ''
    @params['input_dataset_tsv'] = ''
    @params['control_orig'] = ''
    @params['control_other'] = ''
    @params['target_orig'] = ''
    @params['target_other'] = ''
  end
  def set_default_parameters
    count_qc_dir = @dataset_hash.first['Report [File]']
    tpm_txt = File.join(count_qc_dir, "Count_QC-tpm.txt")
    tpm_txt = File.join(@gstore_dir, tpm_txt)
    input_dataset_tsv = File.join(count_qc_dir, "../input_dataset.tsv")
    input_dataset_tsv = File.join(@gstore_dir, input_dataset_tsv)
    @params['input_dataset_tsv'] = input_dataset_tsv
    @params['tpm_txt'] = tpm_txt
    groups = []
    CSV.foreach(input_dataset_tsv, :headers=>true, :col_sep=>"\t") do |row|
      groups <<  row["grouping [Factor]"]
    end
    groups.uniq!
    @params['control_orig'] = groups
    @params['control_other'] = groups
    @params['target_orig'] = groups
    @params['target_other'] = groups
  end
  def preprocess
  end
  def next_dataset
    {'Name'=>@params['name'], 
     'Results [File]'=>File.join(@result_dir, "homeoroq_results")
    }
  end
  def commands
    # make index.csv
    input_dataset_tsv = @params['input_dataset_tsv']
    co = @params['control_orig']
    ct = @params['control_other']
    to = @params['target_orig']
    tt = @params['target_other']

    command = "mkdir homeoroq_results\n"
    command << "ruby /usr/local/ngseq/src/HomeoRoq-1.0/bin/make_index_csv.rb #{input_dataset_tsv} #{co} #{ct} #{to} #{tt}\n"
    command << "cp index.csv homeoroq_results/\n"
    # calcpval_one.R
    @params['iteration'].to_i.times do |i|
      command << "/usr/local/ngseq/stow/R-3.2.2/bin/R --vanilla --slave --args homeoroq_results/tpms_pval_run#{i}.txt #{@params['tpm_txt']} index.csv #{@params['cores']} < /usr/local/ngseq/src/HomeoRoq-1.0/lib/calcpval_one.R\n"
    end
    command << "/usr/local/ngseq/stow/R-3.2.2/bin/Rscript /usr/local/ngseq/src/HomeoRoq-1.0/lib/calcpval_mean.R homeoroq_results\n"
    control = @params['control_orig'].gsub(/_orig/, '')
    target = @params['target_orig'].gsub(/_orig/, '')
    command << "Rscript /usr/local/ngseq/src/HomeoRoq-1.0/lib/plot_homeoroq.R homeoroq_results/pval_mean.xls #{control} #{target}\n"
    command
  end
end


