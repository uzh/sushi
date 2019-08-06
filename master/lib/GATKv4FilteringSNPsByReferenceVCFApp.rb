#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables

class GATKv4FilteringSNPsByReferenceVCFApp <  SushiFabric::SushiApp
  def initialize
    super
    @name = 'GATKv4FilteringSNPsByReferenceVCF'
    @analysis_category = 'Variants'
    @description =<<-EOS
genotype,merge and annotate gvcf-Files<br/>
    EOS
    @required_columns = ['Name', 'Raw VCF', 'Filtered VCF', 'Species', 'refBuild']
    @required_params = ['name']
    @params['cores'] = '1'
    @params['ram'] = '50'
    @params['scratch'] = '100'
    @params['name'] = 'GATKv4Filtering'
    @params['DataSet'] = []
    @modules = ["Dev/Ruby/2.4.3"]
    @inherit_tags = ["Factor", "B-Fabric", "Characteristic"]
  end
  def set_default_parameters
    if data_set = DataSet.find_by_id(@dataset_sushi_id)
      @params['DataSet'] = Hash[*data_set.project.data_sets.map{|d| [d.id.to_s + ":" + d.name, d.id]}.sort_by{|name, id| id}.reverse.flatten]
    end
  end
  def preprocess
  end
  def sample_path(data_set)
    paths = []
    data_set.samples.each do |sample|
      sample.to_hash.each do |header, file|
        if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
          paths << File.dirname(file)
        end
      end
    end
    paths.uniq!
    if path = paths.first
      path.split('/')[0,2].join('/')
    end
  end
  def next_dataset
    {'Name'=>@dataset['Name'],
     'Filtered VCF [File]'=>File.join(@result_dir, "#{@dataset['Name']}.filtered_by_reference_vcf.vcf.gz"),
     'Species'=>@dataset['Species'],
     'refBuild'=>@params['refBuild'],
     'Script [File]'=>File.join(@result_dir, 'filter_out_snps_by_reference_vcf.rb'),
     'Script log [File]'=>File.join(@result_dir, "filter_out_snps_by_reference_vcf.#{@dataset['Name']}.log")
    }
  end
  def commands
    if filter_dataset = DataSet.find_by_id(params['DataSet'])
      filter_vcf_gz = File.join(sample_path(filter_dataset), "GATKv4_Genotyping.filtered.vcf.gz")
      filter_vcf_gz = File.join(@gstore_dir, filter_vcf_gz)
    end
    command =<<-EOS
cat > filter_out_snps_by_reference_vcf.rb <<-EOF
#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190806-113030'

unless target_snps_vcf_gz=ARGV[0] and false_snps_vcf_gz=ARGV[1]
  puts <<-eos
  usage:
   \#{File.basename(__FILE__)} [target.vcf.gz] [false_positive.vcf.gz] (-o [new_filtered.vcf.gz]) (2> log.txt)
  note:
   * without -o option, flat vcf file will be printed in stdout
  eos
  exit
end

out_file = if idx = ARGV.index("-o")
             ARGV[idx+1]
           end

require 'zlib'

false_snps = {}
Zlib::GzipReader.open(false_snps_vcf_gz) do |io|
  io.each do |line|
    # CHROM  POSID REF ALT QUAL  FILTER  INFO  FORMAT  ref_bam_A/Ecor_PR2021.ref
    # ["sa0002", "23965969", ".", "T", "C", "175.77", "PASS", ".", "GT:AD:DP:GQ:PGT:PID:PL", "0/1:2,5:7:69:0|1:23965969_T_C:204,0,69"]
    unless line =~ /^#/
      sid, pos, dot, ref, alt, qual, filter, *others = line.split
      key = [sid, pos].join("_")
      false_snps[key] = alt.split(",")
    end
  end
end

total_snps = 0
false_positives = 0

output =->(out) do
  Zlib::GzipReader.open(target_snps_vcf_gz) do |io|
    io.each do |line|
      unless line =~ /^#/
        sid, pos, dot, ref, alt, qual, filter, *others = line.split
        if filter == "PASS"
          total_snps += 1
          key = [sid, pos].join("_")
          unless (alt.split(",") & false_snps[key].to_a).empty?
            false_positives += 1
          else
            out.print line
          end
        end
      else
        out.print line
      end
    end
  end
end

if out_file and out_file =~ /.gz$/
  Zlib::GzipWriter.open(out_file, &output)
else
  output.(\$stdout)
end

true_positives = total_snps-false_positives
warn
warn "# target vcf:          \#{target_snps_vcf_gz}"
warn "# filter vcf:          \#{false_snps_vcf_gz}"
warn "# total_snps (target): \#{total_snps}"
warn "# false_snps (filter): \#{false_snps.keys.length}"
warn "# false_positives:     \#{false_positives}"
warn "# true_positives:      \#{true_positives}"

precision = true_positives.to_f/(true_positives+false_positives) # PPV (positive predictive value)
fdr = false_positives.to_f/(true_positives+false_positives)      # FDR (false discovery rate)
warn "# precision (TP/total):\#{"%.2f" % (true_positives.to_f/total_snps)}"
warn "# precision (TP/total):\#{"%.2f" % (precision)} (validation)"
warn "# FDR (FP/total):      \#{"%.2f" % (fdr)}"  
warn "# Ref: https://en.wikipedia.org/wiki/Confusion_matrix"
EOF
    EOS
    command << "ruby filter_out_snps_by_reference_vcf.rb #{File.join(@gstore_dir, @dataset['Filtered VCF'])} #{filter_vcf_gz} -o #{@dataset['Name']}.filtered_by_reference_vcf.vcf.gz 2> filter_out_snps_by_reference_vcf.#{@dataset['Name']}.log\n"
    command
  end
end
