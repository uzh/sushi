#!/usr/bin/env ruby
# encoding: utf-8
Version = '20161105-122326'

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class CisTransApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'CisTrans'
    @description =<<-EOS
CisTrans detects cis- and trans-regurated genes by comparing parental diploids and polyploid expression data.<br />
Note<br />
<ol>
<li> Gene-expression divergence (A) between parental species was computed by log2(P1/P2). </li>
<li> Cis effects (B) were estimated by log2(F1P1/F1P2) </li>
<li> Trans effects were derived by subtracting cis effects from the expression divergence between species (A–B)</li>
<li> EdgeR is used to test for differences between F1P1 and F1P2 homeolog expression in F1 (cis effects) </li>
<li> HomeoRow is used to test for differences between parental expression ratio and homeolog expression ratio (trans effects). </li>
</ol>
<br />
Link<br />
<ul>
<li>EdgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html</li>
<li>HomeoRoq: http://seselab.org/homeoroq/</li>
</ul>

http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Polyploid'
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Results']
    @required_params = []
    # optional params
    @params['cores'] = '1'
    @params['ram'] = '10'
    @params['scratch'] = '10'
    @params['name'] = 'CisTrans'
    @params['pval_mean_xls'] = ''

    @params['DataSet'] = []
    @params['DataSet', 'description'] = "EdgeR result between homeologs"
    @required_params = ['DataSet']
  end
  def set_default_parameters
    homeoroq_dir = @dataset_hash.first['Results [File]']
    pval_mean_xls = File.join(homeoroq_dir, "pval_mean.xls")
    pval_mean_xls = File.join(@gstore_dir, pval_mean_xls)
    @params['pval_mean_xls'] = pval_mean_xls

    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['DataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.name, d.id]}.flatten].to_a.reverse
    end
  end
  def preprocess
    dataset_sushi_id = @params['DataSet']
    dataset = DataSet.find_by_id(dataset_sushi_id)
    sample = dataset.samples.first.to_hash
    edger_dir = sample['Report [File]']
    edger_dir = File.join(@gstore_dir, edger_dir)
    @edger_txt = Dir["#{edger_dir}/result*.txt"].to_a.first
  end
  def next_dataset
    {'Name'=>@params['name'], 
     'Results [File]'=>File.join(@result_dir, "cis_trans_results")
    }
  end
  def commands
    pval_mean_xls = @params['pval_mean_xls']

    command = "mkdir cis_trans_results\n"
    command << "ruby /usr/local/ngseq/src/CisTrans-1.0/bin/extract_only_cis_only_trans.rb cis trans #{pval_mean_xls} #{@edger_txt}\n"
    command << "Rscript /usr/local/ngseq/src/CisTrans-1.0/lib/plot_cistrans.R cis_trans_results/cis_trans.csv\n"
    command
  end
end


