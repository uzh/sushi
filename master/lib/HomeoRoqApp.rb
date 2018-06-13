#!/usr/bin/env ruby
# encoding: utf-8
Version = '20180613-135018'

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
<li>Sample name should not have space or -, because it will be '.' in R.</li>
<li>In case of failure, just try again without cds and obh filtering.</li>
</ol>

http://seselab.org/homeoroq/
    EOS

    @analysis_category = 'Polyploid'
    @params['process_mode'] = 'DATASET'
    @required_columns = ['Name','Report', 'Dummy']
    @required_params = []
    # optional params
    @params['cores'] = '8'
    @params['ram'] = '16'
    @params['scratch'] = '10'
    @params['name'] = 'HomeoRoq'
    @params['iteration'] = '10'
    @params['cds_filtering'] = true
    @params['obh_filtering'] = true
    @params['p1_ref'] = ref_selector
    @params['p2_ref'] = ref_selector
    @params['tpm_txt'] = ''
    @params['input_dataset_tsv'] = ''
    @params['control_orig'] = ''
    @params['control_other'] = ''
    @params['target_orig'] = ''
    @params['target_other'] = ''
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
    @modules = ["Tools/RBH/2018.2.8"]
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
      groups <<  row["Condition [Factor]"]
    end
    groups.uniq!
    @params['control_orig'] = groups
    @params['control_orig', 'description'] = "Diploid orig"
    @params['control_other'] = groups
    @params['control_other', 'description'] = "Diploid other"
    @params['target_orig'] = groups
    @params['target_orig', 'description'] = "Polyploid orig"
    @params['target_other'] = groups
    @params['target_other', 'description'] = "Polyploid other"
    @modules = ["Dev/R"]
  end
  def preprocess
  end
  def next_dataset
    {'Name'=>@params['name'], 
     'Results [File]'=>File.join(@result_dir, "homeoroq_results")
    }
  end
  def commands
    genomes_root = "/srv/GT/reference"
    # make index.csv
    input_dataset_tsv = @params['input_dataset_tsv']
    co = @params['control_orig']
    ct = @params['control_other']
    to = @params['target_orig']
    tt = @params['target_other']

    # filtering
    tpm_txt = @params['tpm_txt']
    command = "mkdir homeoroq_results\n"
    command << "echo 'HomeoRoq version:'\n"
    command << "ls -l /usr/local/ngseq/bin/make_index_csv.rb\n"
    if @params['cds_filtering'] 
      # output[['refBuild']] = 'Cardamine_hirsuta/KEN/DENOVO/Annotation/Version-2014-06-29'
      # ruby scripts/filter_tpm_by_cds_obh.rb -t data/Count_QC-tpm.txt -1 references/chir_genes.gtf -2 references/cama_genes.gtf -a references/chir_genome.fa -b references/cama_genome.fa
      p1_ref = File.join(genomes_root, @params['p1_ref'])
      p1_gtf = File.join(p1_ref, "Genes/genes.gtf")
      if @params['obh_filtering']
        p2_ref = File.join(genomes_root, @params['p2_ref'])
        p2_gtf = File.join(p2_ref, "Genes/genes.gtf")
        p1_genome_dir = File.expand_path('../../Sequence/WholeGenomeFasta', p1_ref)
        p1_genome = File.join(p1_genome_dir, "genome.fa")
        p2_genome_dir = File.expand_path('../../Sequence/WholeGenomeFasta', p2_ref)
        p2_genome = File.join(p2_genome_dir, "genome.fa")
        command << "ruby /usr/local/ngseq/lib/filter_tpm_by_cds_obh.rb -t #{@params['tpm_txt']} -1 #{p1_gtf} -2 #{p2_gtf} -a #{p1_genome} -b #{p2_genome}\n"
        command << "cp out.tsv homeoroq_results/filtered_tpm.txt\n"
        tpm_txt = "out.tsv"
      else
        command << "ruby /usr/local/ngseq/lib/filter_tpm_by_cds_obh.rb -t #{@params['tpm_txt']} -1 #{p1_gtf} --no-obh-filtering\n"
        command << "cp out.tsv homeoroq_results/filtered_tpm.txt\n"
        tpm_txt = "out.tsv"
      end
    else
      if @params['obh_filtering']
        p1_ref = File.join(genomes_root, @params['p1_ref'])
        p1_gtf = File.join(p1_ref, "Genes/genes.gtf")
        p2_ref = File.join(genomes_root, @params['p2_ref'])
        p2_gtf = File.join(p2_ref, "Genes/genes.gtf")
        p1_genome_dir = File.expand_path('../../Sequence/WholeGenomeFasta', p1_ref)
        p1_genome = File.join(p1_genome_dir, "genome.fa")
        p2_genome_dir = File.expand_path('../../Sequence/WholeGenomeFasta', p2_ref)
        p2_genome = File.join(p2_genome_dir, "genome.fa")

        command << "ruby /usr/local/ngseq/lib/filter_tpm_by_cds_obh.rb -t #{@params['tpm_txt']} -1 #{p1_gtf} -2 #{p2_gtf} -a #{p1_genome} -b #{p2_genome} --no-cds-filtering\n"
        command << "cp out.tsv homeoroq_results/filtered_tpm.txt\n"
        tpm_txt = "out.tsv"
      end
    end
    command << "ruby /usr/local/ngseq/bin/make_index_csv.rb #{input_dataset_tsv} #{co} #{ct} #{to} #{tt}\n"
    command << "cp index.csv homeoroq_results/\n"
    # calcpval_one.R
    @params['iteration'].to_i.times do |i|
      command << "/usr/local/ngseq/stow/R-3.2.2/bin/R --vanilla --slave --args homeoroq_results/tpms_pval_run#{i}.txt #{tpm_txt} index.csv #{@params['cores']} < /usr/local/ngseq/lib/calcpval_one.R\n"
    end
    command << "/usr/local/ngseq/stow/R-3.2.2/bin/Rscript /usr/local/ngseq/lib/calcpval_mean.R homeoroq_results\n"
    control = @params['control_orig'].gsub(/_orig/, '')
    target = @params['target_orig'].gsub(/_orig/, '')
    command << "Rscript /usr/local/ngseq/lib/plot_homeoroq.R homeoroq_results/pval_mean.xls #{control} #{target}\n"
    command
  end
end


