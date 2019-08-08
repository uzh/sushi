#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKv4GenerateTranscriptsFromVCFApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKv4GenerateTranscriptsFromVCF'
    @analysis_category = 'Variants'
    @description =<<-EOS
filtering out SNPs by the VCF coming from reference accession<br/>
    EOS
    @required_columns = ['Name', 'Filtered VCF', 'Species', 'refBuild']
    @required_params = []
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @modules = ["Dev/Ruby/2.4.3", "Tools/Cufflinks/2.2.1"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
  end
  def preprocess
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'New Genome Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genome.fa.gz"),
     'Transcripts Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.cds.fa.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@dataset['refBuild'],
     'Script [File]'=>File.join(@result_dir, "replace_snps_by_vcf.#{@dataset['Name']}.rb"),
     'Script log [File]'=>File.join(@result_dir, "replace_snps_by_vcf.#{@dataset['Name']}.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    command =<<-EOS
#!/bin/bash
# Version = '20190808-085705'

cat > replace_snps_by_vcf.#{@dataset['Name']}.rb <<-EOF
#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190808-085705'

unless genome_fa=ARGV[0] and filtered_vcf_gz=ARGV[1]
  puts <<-eos
  usage:
   \#{File.basename(__FILE__)} [genome.fa] [filtered.vcf.gz] > new.fa
  note:
   * Site becomes N if
    * hetero site (GT other than 1|1 or 0|0)
    * GQ < 20
    * DP<2 or DP>250
   * Site replaced from reference if
    * hard filtering == PASS
  eos
  exit
end

require 'zlib'

replaces = {}
warn "# loading \#{filtered_vcf_gz}: \#{Time.now}"
Zlib::GzipReader.open(filtered_vcf_gz) do |gz|
  while line=gz.gets
    unless line =~ /^#/
      # sa0001  12404264  . C T 37.60 PASSAC=1;AF=0.500;AN=2;BaseQRankSum=0.965;DP=14;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=58.75;MQRankSum=-3.328e+00;QD=2.69;ReadPosRankSum=-8.680e-01;SOR=1.022 GT:AD:DP:GQ:PGT:PID:PL:PS 0|1:13,1:14:45:0|1:12404260_C_T:45,0,793:12404260
      sid, pos, dot, ref, alt, qual, filter, info,  format, values = line.chomp.split
      ac_, af_, an_, bq_, dp_, *others = info.split(";")
      af = af_.split("=").last
      gt, ad, dp, gq, *others = values.split(":")
      #if alt.split(",").length > 1
      if filter == "PASS"
        if alt.split(",").length > 1 or
           af.split(",").first.to_f < 1.0 or
           gq.to_i < 20 or
           dp.to_i < 2 or
           dp.to_i > 250
          # N
          # puts [sid, pos, alt, af, gq, gt, dp].join("\\t")
          replaces[sid] ||= {}
          replaces[sid][pos.to_i] = "N"
        else
          # replace
          # puts [sid, pos, alt, af, gt].join("\\t")
          replaces[sid] ||= {}
          replaces[sid][pos.to_i] = alt
        end
      end
    end
  end
end

#p replaces["sa0071"]
#exit
# p replaces
require 'bio'
Bio::FlatFile.open(genome_fa).each do |e|
  sid = e.definition
  snps = 0
  ns   = 0
  new_seq = e.seq.clone
  warn "# processing \#{sid}: \#{Time.now}"
  replaces[sid]&.each do |pos, alt|
    if new_seq[pos-1] == alt
      raise "WWW"
    end
    new_seq[pos-1] = alt
    snps += 1
    if alt == "N"
      ns += 1
    end
  end
  warn "# \#{sid}: replaced SNPs \#{snps} (unclear SNPs \#{ns})"
  puts ">\#{sid}"
  puts new_seq.scan(/.{100}|.+\\Z/).join("\\n")
end
EOF
    EOS
    script_rb = "replace_snps_by_vcf.#{@dataset['Name']}.rb"
    script_log = "replace_snps_by_vcf.#{@dataset['Name']}.log"
    genome_fa = File.join(GENOME_REF_DIR, @dataset['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    filtered_vcf = File.join(@gstore_dir, @dataset['Filtered VCF'])
    refbuild_path = File.join(GENOME_REF_DIR, @dataset['refBuild'])
    genes_gtf = File.join(refbuild_path, "Genes/genes.gtf")
    new_genome_fa  = "#{@dataset['Name']}.genome.fa"
    transcripts_fa = "#{@dataset['Name']}.cds.fa"
    command << "ruby #{script_rb} #{genome_fa} #{filtered_vcf} > #{new_genome_fa} 2> #{script_log}\n"
    command << "gffread -x #{transcripts_fa} -g #{new_genome_fa} #{genes_gtf}\n" 
    command << "gzip #{new_genome_fa}\n"
    command << "gzip #{transcripts_fa}\n"
    command
  end
end
