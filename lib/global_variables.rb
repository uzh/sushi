#!/usr/bin/env ruby
# encoding: utf-8

module GlobalVariables
  SUSHI = 'Supercalifragilisticexpialidocious!!'
  SAMTOOLS = '/usr/local/ngseq/stow/samtools-0.1.19/bin/samtools'
  BCFTOOLS = '/usr/local/ngseq/stow/samtools-0.1.19/bin/bcftools'
  PICARD_DIR = '/usr/local/ngseq/stow/picard-tools-1.108/bin'
  GATK_DIR = '/usr/local/ngseq/src/GenomeAnalysisTK-3.2-2'
  SNPEFF_DIR='/usr/local/ngseq/src/snpEff_v3.4/'
  BWA='/usr/local/ngseq/stow/bwa-0.7.10/bin/bwa'
  SAMSTAT='/usr/local/ngseq/stow/samstat_1.09/bin/samstat'
  QUALIMAP='/usr/local/ngseq/opt/qualimap_v1.0/qualimap'
  TRANSDECODER='/usr/local/ngseq/src/trinityrnaseq_r20131110/trinity-plugins/transdecoder/TransDecoder'
  HMMSCAN='/usr/local/ngseq/src/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmscan'
  UNIPROT_DIR='/srv/GT/databases/uniprot_blast_db/'
  BLAST_DIR='/usr/local/ngseq/src/ncbi-blast-2.2.26+/bin/'
  AUGUSTUS_DIR='/usr/local/ngseq/src/augustus.2.7/'
  HUMAN_SNP_DATABASES='/srv/GT/reference/Homo_sapiens/UCSC/hg19_karyotypic/SNPS/'
  def builder_selector(base_dir, shown_pattern=nil, value_pattern=nil)
    selector = {}
    Dir[base_dir].sort.select{|dir| File.directory?(dir)}.each do |dir|
      key = if shown_pattern
              dir.gsub(shown_pattern.keys.first, shown_pattern)
            else
              dir
            end
      value = if value_pattern
                dir.gsub(value_pattern.keys.first, value_pattern)
              else
                File.basename(dir)
              end
      selector[key] = value
    end
    selector
  end
  def ref_selector
    selector = {'select'=>''}

    base_pattern = "/srv/GT/reference/*/*/*"
    shown_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    value_replace_pattern = {/\/srv\/GT\/reference\//=>''}
    builds = builder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

    base_pattern = "/srv/GT/reference/*/*/*/Annotation/Version*"
    versions = builder_selector(base_pattern, shown_replace_pattern, value_replace_pattern)

    builds.keys.each do |build_key|
      if versions.keys.find{|version_key| version_key=~/#{build_key}/}
        builds.delete(build_key)
      end
    end
    builds.merge!(versions)
    builds = builds.sort.to_h

    selector.merge!(builds)
    selector
  end
  def factor_dataset
    factors = get_columns_with_tag 'Factor'
    dataset = {}
    factors.first.keys.each do |colname|
      dataset[colname+" [Factor]"] = @dataset[colname]
    end
    dataset
  rescue
    {}
  end
  def run(class_name)
    opt = OptionParser.new do |o|
      o.banner = "Usage: ruby #{__FILE__} [options]"
      o.on(:dataset_id, '-i dataset_id', '--dataset_id', Integer, 'DataSet ID in Sushi DB')
      o.on(:dataset_tsv, '-d dataset_tsv', '--dataset', String, 'DataSet file (.tsv) (This option is prior to dataset_id option)')
      o.on(:parameterset_tsv, '-m parameterset_tsv', '--parameterset', String, 'Parameterset file (.tsv)')
      o.on(:run_mode, '-r', '--run', 'Real run mode. without this option, it runs with test run mode which checks only DataSet and Parameters and no submittion')
      o.on(:project, '1001', '-p', '--project', String, 'Project Number (default: 1001)')
      o.on(:user, 'sushi_lover', '-u', '--user', String, 'Submit user (default: sushi_lover)')
      o.parse!(ARGV)
    end
    opt.project = 'p' + opt.project
    unless opt.dataset_id or opt.dataset_tsv
      puts
      warn "\e[31mERROR\e[0m: Either dataset_id or dataset_tsv is required"
      print opt.help
      exit
    end

    usecase = class_name.new

    usecase.project = opt.project
    usecase.user = opt.user
    usecase.dataset_sushi_id = opt.dataset_id
    usecase.dataset_tsv_file = opt.dataset_tsv
    usecase.parameterset_tsv_file = opt.parameterset_tsv
    if opt.run_mode
      usecase.run
    else
      usecase.test_run
    end
  end

end
