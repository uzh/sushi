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
    @required_columns = ['Name', 'BAM', 'Filtered VCF', 'Species', 'refBuild']
    @required_params = []
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '30'
    @modules = ["Dev/Ruby/2.4.3", "Tools/Cufflinks/2.2.1", "Tools/samtools/1.9"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
  end
  def preprocess
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'New Genome Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genome.fa.gz"),
     'New Genome Masked Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genome.masked.fa.gz"),
     'New Genes Annotation [File]'=>File.join(@result_dir, "#{@dataset['Name']}.genes.gtf.gz"),
     'Transcripts Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.cds.fa.gz"),
     'Transcripts Masked Fasta [File]'=>File.join(@result_dir, "#{@dataset['Name']}.cds.masked.fa.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@dataset['refBuild'],
     'Script1 [File]'=>File.join(@result_dir, "replace_N_with_low_high_coverage.#{@dataset['Name']}.rb"),
     'Script2 [File]'=>File.join(@result_dir, "replace_snps_indels_by_vcf.#{@dataset['Name']}.rb"),
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    command =<<-EOS
#!/bin/bash
# Version = '20190822-060725'

cat > replace_N_with_low_high_coverage.#{@dataset['Name']}.rb <<-EOF1

#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190822-060725'

unless bam_or_depth=ARGV[0] and genome_fa=ARGV[1]
  puts <<-eos
  usage:
   \#{File.basename(__FILE__)} [target.bam|samtools.depth] [genome.fa] (-o samtools.depth) > new_genome.fa

  note:
   .bam: the positions in coverage < 2 or coverage < 250 in samtools depth will be replaced to N in genome.fa
   .depth: all the positions will be replaced to N in genome.fa

  option:
   -o: keep samtools.depth
  eos
  exit
end

#command = "samtools depth -aa -r sa0001 \#{bam}"
sid2pos = {}
if bam_or_depth =~ /\\.bam$/

  bam = bam_or_depth
  bai = bam + ".bai"
  unless File.exist?(bai)
    command = "samtools index \#{bam}"
    warn "# \`which samtools\`"
    warn "# \#{Time.now}: \#{command}"
    system command
  end
  out = if idx = ARGV.index("-o")
          out_file = ARGV[idx+1]
          open(out_file, "w")
        end
  command = "samtools depth -aa \#{bam}"
  warn "# \#{Time.now}: \#{command}"
  IO.popen(command) do |io|
    while line=io.gets
      sid, pos, cov_ = line.chomp.split
      cov = cov_.to_i
      if cov < 2 or 250 < cov
        sid2pos[sid] ||= []
        sid2pos[sid] << pos.to_i
        out.print line if out
      end
    end
  end

elsif bam_or_depth =~ /\\.depth$/
  samtools_depth = bam_or_depth
  warn "# \#{Time.now}: loading: \#{samtools_depth}"
  File.readlines(samtools_depth).each do |line|
    sid, pos, cov = line.chomp.split
    sid2pos[sid] ||= []
    sid2pos[sid] << pos.to_i
  end

end

require 'bio'
warn "# \#{Time.now}: loading: \#{genome_fa}"
sid2seq = {}
Bio::FlatFile.open(genome_fa).each do |e|
  sid2seq[e.definition] = e.seq
end


total_bases = sid2pos.select{|sid, pos| sid2seq.keys.include?(sid)}.values.flatten.length
warn "# \#{Time.now}: start masking: Total \#{total_bases} bases replacement"

replaced_bases = 0
sid2pos.each do |sid, poss|
  if seq = sid2seq[sid]
    poss.each do |pos|
      seq[pos-1] = "N"
    end
    replaced_bases += poss.length
    warn "# \#{Time.now}: \#{sid} done, \#{poss.length} bases replaced with N (\#{"%.1f" % (replaced_bases/total_bases*100.0)}%)"
  end
end

sid2seq.each do |sid, seq|
  puts ">\#{sid}"
  puts seq.scan(/.{100}|.+\\Z/).join("\\n")
end
warn "# \#{Time.now}: done"

EOF1

cat > replace_snps_indels_by_vcf.#{@dataset['Name']}.rb <<-EOF2

#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190820-160716'

unless genome_fa=ARGV[0] and genome_masked_fa=ARGV[1] and genes_gtf=ARGV[2] and filtered_vcf_gz=ARGV[3]
  puts <<-eos
  usage:
   \#{File.basename(__FILE__)} [genome.fa] [genome.masked.fa] [genes.gtf] [filtered.vcf.gz] (-f genome.new.fa -m genome.masked.new.fa -t genes.new.gtf)
  options:
   without -f, -m or -t, the default file names will be *.new.fa, *.new.gtf in the current directly
  note:
   * masked_genome.fa: N replaced by low/high coverage
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

out_genome_masked_fa = if idx = ARGV.index("-m")
                  ARGV[idx+1]
                else
                  File.basename(genome_masked_fa).gsub(".fa", ".new.fa")
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

# puts "# check: replaces['sa0071'] = \#{replaces["sa0071"]}"
# exit
# p replaces.values.length
# p replaces_indel.values.length
# exit

require 'bio'

warn "# \#{Time.now}: loading \#{genome_masked_fa}"
sid2seq_masked = {}
Bio::FlatFile.open(genome_masked_fa).each do |e|
  sid = e.definition
  sid2seq_masked[sid] = e.seq
end

warn "# \#{Time.now}: processing \#{genome_fa}"
open(out_genome_fa, "w") do |out|
open(out_genome_masked_fa, "w") do |out_masked|
  Bio::FlatFile.open(genome_fa).each do |e|
    sid = e.definition
    snps = 0
    indels = 0
    ns   = 0
    ni   = 0
    new_seq = e.seq.clone
    new_seq_masked = sid2seq_masked[sid]
    warn "# \#{Time.now}: processing \#{sid}"
    replaces[sid]&.to_a&.reverse&.each do |pos, ref_alt|
      ref, alt = ref_alt
      unless ref == new_seq[pos-1, ref.length]
        warn "# WARN: \#{sid}:\#{pos}"
        warn "#   ref org: \#{ref}" 
        warn "#   ref new: \#{new_seq[pos-1, ref.length]}"
        warn "#   alt new: \#{alt}"
      end
      new_seq[pos-1, ref.length] = alt
      new_seq_masked[pos-1, ref.length] = alt
      if ref.length == 1 and alt.length == 1
        snps += 1
      else
        indels += 1
        warn "# InDel: \#{[sid, pos, ref_alt].join(",")}"
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
    out.puts new_seq.scan(/.{100}|.+\Z/).join("\\n")
    out_masked.puts ">\#{sid}"
    out_masked.puts new_seq_masked.scan(/.{100}|.+\Z/).join("\\n")
  end
end
end
warn "# \#{Time.now}: generated \#{out_genome_fa}"
warn "# \#{Time.now}: generated \#{out_genome_masked_fa}"

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
warn "# \#{Time.now}: generated \#{out_genes_gtf}"

EOF2
    EOS
    script1_rb = "replace_N_with_low_high_coverage.#{@dataset['Name']}.rb"
    genome_fa = File.join(GENOME_REF_DIR, @dataset['refBuild'].split('/')[0,3].join("/")+"/Sequence/WholeGenomeFasta/genome.fa")
    gstore_bam = File.join(@gstore_dir, @dataset['BAM'])
    bam = File.basename(@dataset['BAM'])
    org_genome_masked_fa = "#{@dataset['Name']}.genome.masked.org.fa"
    command << "cp #{gstore_bam} ./#{bam}\n"
    # ruby replace_N_with_low_high_coverage.rb data/Ecor_GE12_DENOVO_v2.0_A_subgenome_ref.bam references/sa0050.fa > sa0050.new2.fa
    command << "ruby #{script1_rb} #{bam} #{genome_fa} > #{org_genome_masked_fa}\n"

    script2_rb = "replace_snps_indels_by_vcf.#{@dataset['Name']}.rb"
    filtered_vcf = File.join(@gstore_dir, @dataset['Filtered VCF'])
    refbuild_path = File.join(GENOME_REF_DIR, @dataset['refBuild'])
    genes_gtf = File.join(refbuild_path, "Genes/genes.gtf")

    new_genome_fa  = "#{@dataset['Name']}.genome.fa"
    new_genome_masked_fa = "#{@dataset['Name']}.genome.masked.fa"
    new_genes_gtf = "#{@dataset['Name']}.genes.gtf"
    transcripts_fa = "#{@dataset['Name']}.cds.fa"
    transcripts_masked_fa = "#{@dataset['Name']}.cds.masked.fa"
    # ruby replace_snps_indels_by_vcf.rb references/sa0049.fa genome.sa0049.masked.fa references/genes.gtf A_Ecor_GE12.snps_indels.sa0049.vcf.gz -f genome.sa0049.new.fa -m genome.sa0049.masked.new.fa -t genes.new.gtf
    command << "ruby #{script2_rb} #{genome_fa} #{org_genome_masked_fa} #{genes_gtf} #{filtered_vcf} -f #{new_genome_fa} -m #{new_genome_masked_fa} -t #{new_genes_gtf}\n"
    command << "gffread -x #{transcripts_fa} -g #{new_genome_fa} #{new_genes_gtf}\n" 
    command << "gffread -x #{transcripts_masked_fa} -g #{new_genome_masked_fa} #{new_genes_gtf}\n" 
    command << "gzip #{new_genome_fa}\n"
    command << "gzip #{new_genome_masked_fa}\n"
    command << "gzip #{new_genes_gtf}\n"
    command << "gzip #{transcripts_fa}\n"
    command << "gzip #{transcripts_masked_fa}\n"
    command
  end
end
