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
    @params['scratch'] = '30'
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
     'New Genes Annotation [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genes.gtf.gz"),
     'Transcripts Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.cds.fa.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@dataset['refBuild'],
     'Script [File]'=>File.join(@result_dir, "replace_snps_indels_by_vcf.#{@dataset['Name']}.rb"),
     'Script log [File]'=>File.join(@result_dir, "replace_snps_indels_by_vcf.#{@dataset['Name']}.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    command =<<-EOS
#!/bin/bash
# Version = '20190813-152717'

cat > replace_snps_indels_by_vcf.#{@dataset['Name']}.rb <<-EOF
#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190813-152717'

unless genome_fa=ARGV[0] and genes_gtf=ARGV[1] and filtered_vcf_gz=ARGV[2]
  puts <<-eos
  usage:
   \#{File.basename(__FILE__)} [genome.fa] [genes.gtf] [filtered.vcf.gz] (-f genome.new.fa -t genes.new.gtf)
  options:
   without -f or -t, the default file names will be *.new.fa, *.new.gtf in the current directly
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

out_genome_fa = if idx = ARGV.index("-f")
                  ARGV[idx+1]
                else
                  File.basename(genome_fa).gsub(".fa", ".new.fa")
                end
out_genes_gtf = if idx = ARGV.index("-t")
                  ARGV[idx+1]
                else
                  File.basename(genes_gtf).gsub(".gtf", ".new.gtf")
                end

require 'zlib'

warn "# \#{Time.now}: loading \#{genes_gtf} "
genes_gtf_lines = []
File.readlines(genes_gtf).each do |line|
  genes_gtf_lines << line
end

replaces = {}
replaces_indel = {}
warn "# \#{Time.now}: loading \#{filtered_vcf_gz}"
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
          # puts [sid, pos, ref, alt, af, gq, gt, dp].join("\\t")
          replaces[sid] ||= {}
          new_alt = "N"*ref.length
          replaces[sid][pos.to_i] = [ref, new_alt]
        elsif ref.split(",").length > 1
          puts "WWW"
          exit
        else
          # replace
          # puts [sid, pos, alt, af, gt].join("\\t")
          replaces[sid] ||= {}
          replaces[sid][pos.to_i] = [ref, alt]
          unless ref.length == alt.length
            replaces_indel[sid] ||= {}
            replaces_indel[sid][pos.to_i] = [ref, alt]
          end
        end
      end
    end
  end
end

warn "# \#{Time.now}: processing \#{genome_fa}"
require 'bio'
open(out_genome_fa, "w") do |out|
  Bio::FlatFile.open(genome_fa).each do |e|
    sid = e.definition
    snps = 0
    indels = 0
    ns   = 0
    ni   = 0
    new_seq = e.seq.clone
    warn "# \#{Time.now}: processing \#{sid}"
    replaces[sid]&.to_a.reverse.each do |pos, ref_alt|
      ref, alt = ref_alt
      unless ref == new_seq[pos-1, ref.length]
        warn "# WARN: \#{sid}:\#{pos}"
        warn "#   ref org: \#{ref}" 
        warn "#   ref new: \#{new_seq[pos-1, ref.length]}"
        warn "#   alt new: \#{alt}"
      end
      new_seq[pos-1, ref.length] = alt
      if ref.length == 1 and alt.length == 1
        snps += 1
      else
        indels += 1
      end
      if alt =~ /N/
        if ref.length == 1 and alt.length == 1
          ns += 1
        else
          ni += 1
        end
      end
    end
    warn "# \#{sid}: replaced SNPs:   \#{snps}\\t(unclear SNPs: \#{ns})"
    warn "# \#{sid}: replaced InDels: \#{indels}\\t(unclear InDels: \#{ni})"
    out.puts ">\#{sid}"
    out.puts new_seq.scan(/.{100}|.+\\Z/).join("\\n")
  end
end

warn "# \#{Time.now}: processing \#{genes_gtf}"
replaces_i = 0
open(out_genes_gtf, "w") do |out|
  genes_gtf_lines.each do |line|
    unless line =~ /^#/
      # sa0001  maker exon  84805 85236 . + . transcript_id "ga00001"; gene_id "ga00001"; gene_name "ga00001"; gene_biotype "protein_coding"
      # sa0001  maker CDS 84805 85236 . + 0 transcript_id "ga00001"; gene_id "ga00001"; gene_name "ga00001"; gene_biotype "protein_coding"
      sid, maker, type, st_, ed_, dot1, strand, dot2, *others = line.chomp.split
      st = st_.to_i
      ed = ed_.to_i
      until_st = 0
      between_st_ed = 0
      replaces_indel[sid]&.each do |pos, ref_alt|
        ref, alt = ref_alt
        if pos < st
          until_st += (alt.length - ref.length)
        elsif pos <= ed
          between_st_ed += (alt.length - ref.length)
        else
          break
        end
      end
      new_st = st + until_st
      new_ed = ed + until_st + between_st_ed
      out.puts [sid, maker, type, new_st, new_ed, dot1, strand, dot2, others.join(" ")].join("\\t")
    else
      out.print line
    end
  end
end
EOF
    EOS
    script_rb = "replace_snps_indels_by_vcf.#{@dataset['Name']}.rb"
    script_log = "replace_snps_indels_by_vcf.#{@dataset['Name']}.log"
    genome_fa = File.join(GENOME_REF_DIR, @dataset['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    filtered_vcf = File.join(@gstore_dir, @dataset['Filtered VCF'])
    refbuild_path = File.join(GENOME_REF_DIR, @dataset['refBuild'])
    genes_gtf = File.join(refbuild_path, "Genes/genes.gtf")
    new_genome_fa  = "#{@dataset['Name']}.genome.fa"
    new_genes_gtf = "#{@dataset['Name']}.genes.gtf"
    transcripts_fa = "#{@dataset['Name']}.cds.fa"
    command << "ruby #{script_rb} #{genome_fa} #{genes_gtf} #{filtered_vcf} -f #{new_genome_fa} -t #{new_genes_gtf}  2> #{script_log}\n"
    command << "gffread -x #{transcripts_fa} -g #{new_genome_fa} #{new_genes_gtf}\n" 
    command << "gzip #{new_genome_fa}\n"
    command << "gzip #{new_genes_gtf}\n"
    command << "gzip #{transcripts_fa}\n"
    command
  end
end
