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
filtering out SNPs by the VCF coming from reference accession<br/>
    EOS
    @required_columns = ['Name', 'Raw VCF', 'Filtered VCF', 'Species', 'refBuild']
    @required_params = []
    @params['cores'] = '1'
    @params['cores', \"context\"] = \"slurm\"
    @params['ram'] = '50'
    @params['ram', \"context\"] = \"slurm\"
    @params['scratch'] = '100'
    @params['scratch', "context"] = "slurm"
    @params['DataSet'] = []
    @params['DataSet', "context"] = "GATKv4FilteringSNPsByReferenceVCF"
    @modules = ["Dev/Ruby/3.1.3"]
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
     'refBuild'=>@dataset['refBuild'],
     'Script [File]'=>File.join(@result_dir, "filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.rb"),
     'Script log [File]'=>File.join(@result_dir, "filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.log")
    }.merge(extract_columns(@inherit_tags))
  end
  def commands
    if filter_dataset = DataSet.find_by_id(params['DataSet']) and
       filter_vcf_gz_dir = File.join(@gstore_dir, sample_path(filter_dataset)) and
       filter_vcf_gz = Dir[File.join(filter_vcf_gz_dir, "*.filtered.vcf.gz")].to_a.first

      command =<<-EOS
cat > filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.rb <<-EOF
#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20190812-110125'

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

false_snps_indels = {}
false_snps = 0
false_indels = 0
Zlib::GzipReader.open(false_snps_vcf_gz) do |io|
  io.each do |line|
    # CHROM  POSID REF ALT QUAL  FILTER  INFO  FORMAT  ref_bam_A/Ecor_PR2021.ref
    # ["sa0002", "23965969", ".", "T", "C", "175.77", "PASS", ".", "GT:AD:DP:GQ:PGT:PID:PL", "0/1:2,5:7:69:0|1:23965969_T_C:204,0,69"]
    unless line =~ /^#/
      sid, pos, dot, ref, alt, qual, filter, *others = line.split
      key = [sid, pos].join("_")
      false_snps_indels[key] = alt.split(",")
      if alt.split(",").map{|a| a.length}.max > 1 or ref.split(",").map{|r| r.length}.max > 1
        false_indels += 1
      else
        false_snps += 1
      end
    end
  end
end

total_snps = 0
false_positive_snps = 0
total_indels = 0
false_positive_indels = 0

output =->(out) do
  Zlib::GzipReader.open(target_snps_vcf_gz) do |io|
    io.each do |line|
      unless line =~ /^#/
        sid, pos, dot, ref, alt, qual, filter, *others = line.split
        flag_indel = false
        if filter == "PASS"
          if alt.split(",").map{|a| a.length}.max > 1 or ref.split(",").map{|r| r.length}.max > 1
            total_indels += 1
            flag_indel = true
          else
            total_snps += 1
          end
          key = [sid, pos].join("_")
          #unless (alt.split(",") & false_snps[key].to_a).empty?
          if (alt.split(",") - false_snps_indels[key].to_a).empty?
            if flag_indel
              false_positive_indels += 1
            else
              false_positive_snps += 1
            end
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

total_snps_indels = total_snps + total_indels
warn
warn "# target vcf:                         \#{target_snps_vcf_gz}"
warn "# filter vcf:                         \#{false_snps_vcf_gz}"
warn "# total_snps_indels (target):         \#{total_snps_indels}"
warn "#   total_snps (target):              \#{total_snps}"
warn "#   total_indels (target):            \#{total_indels}"
warn "# false_snps_indels (filter):         \#{false_snps_indels.keys.length}"
warn "#   false_snps (filter):              \#{false_snps}"
warn "#   false_indels (filter):            \#{false_indels}"
warn

false_positives = false_positive_snps + false_positive_indels
true_positives = total_snps_indels-false_positives
warn "# false_positives (filtered):         \#{false_positives}"
warn "#   false_positive_snps (filtered):   \#{false_positive_snps}"
warn "#   false_positive_indels (filtered): \#{false_positive_indels}"
warn "# true_positives:                     \#{true_positives}"
warn "#   = total_snps_indels - false_positives"
warn "#   = \#{total_snps_indels} - \#{false_positives} = \#{total_snps_indels-false_positives}"
warn

precision = true_positives.to_f/(true_positives+false_positives) # PPV (positive predictive value)
fdr = false_positives.to_f/(true_positives+false_positives)      # FDR (false discovery rate)
warn "# precision (TP/total):             \#{"%.2f" % (true_positives.to_f/total_snps_indels)}"
warn "# precision (TP/total):             \#{"%.2f" % (precision)} (validation)"
warn "# FDR (FP/total):                   \#{"%.2f" % (fdr)}"  
warn "# Ref: https://en.wikipedia.org/wiki/Confusion_matrix"
EOF
      EOS
      command << "ruby filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.rb #{File.join(@gstore_dir, @dataset['Filtered VCF'])} #{filter_vcf_gz} -o #{@dataset['Name']}.filtered_by_reference_vcf.vcf.gz 2> filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.log\n"
      command
    else
      command =  "echo 'Error'\n"
      command << "touch filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.rb\n"
      command << "touch #{@dataset['Name']}.filtered_by_reference_vcf.vcf.gz\n"
      command << "touch filter_out_snps_indels_by_reference_vcf.#{@dataset['Name']}.log\n"
    end
  end
end
