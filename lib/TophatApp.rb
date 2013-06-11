#!/usr/bin/env ruby
# encoding: utf-8
Version = '20130610-163048'

require 'sushiApp'

class TophatApp < SushiApp
  def initialize
    super
    @name = 'Tophat'
    @analysis_category = 'Map'
    @required_columns = ['Sample','Read1','Species']
    @required_params = ['build','paired','cores']
    # optional params
    @params['is_stranded'] = ''
    @params['paired'] = false
    @params['build'] = ''
    @output_files = ['BAM','BAI']
  end
  def preprocess
    if @params['paired']
      @required_columns << 'Read2'
    end
  end
  def next_dataset
    {'Sample'=>@dataset['Sample'], 
     'BAM'=>File.join(@result_dir, "#{@dataset['Sample']}.bam"), 
     'BAI'=>File.join(@result_dir, "#{@dataset['Sample']}.bam.bai"),
     'Build'=>@params['build']
    }
  end
  def build_dir
    @build_dir ||= Dir["/srv/GT/reference/*/*/#{@params['build']}"].to_a[0]
  end
  def bowtie2_index
    @bowtie2_index ||= if build_dir
                         File.join(build_dir, 'Sequence/BOWTIE2Index/genome')
                       end
  end
  def transcripts_index
    @transcripts_index ||= if build_dir
                             File.join(build_dir, 'Annotation/Genes/genes_BOWTIE2Index/transcripts')
                           end
  end
  def library_type
    if @params['is_stranded'].empty?
      "--library-type fr-unstranded"
    else
     if  @params['is_stranded'] == 'sense'
      "--library-type fr-secondstrand"
     else
       "--library-type fr-firststrand"
     end
    end
  end
  def commands
    if bowtie2_index and transcripts_index
      command = "/usr/local/ngseq/bin/tophat -o . --num-threads #{@params['cores']} #{library_type} --transcriptome-index #{transcripts_index} #{bowtie2_index} $WORKSPACE_DIR/#{@dataset['Read1']}"
      if @params['paired']
        command << ",$WORKSPACE_DIR/#{@dataset['Read2']}\n"
      else
        command << "\n"
      end
      command << "mv accepted_hits.bam #{@dataset['Sample']}.bam\n"
      command << "samtools index #{@dataset['Sample']}.bam\n"
    end
    command
  end
end

if __FILE__ == $0
  usecase = TophatApp.new

  usecase.project = "p1001"

  # set user parameter
  # for GUI sushi
  #usecase.params['process_mode'].value = 'SAMPLE'
  #usecase.params['build'] = 'TAIR10'
  #usecase.params['paired'] = true
  #usecase.params['cores'] = 2
  #usecase.params['node'] = 'fgcz-c-048'

  # also possible to load a parameterset csv file
  # mainly for CUI sushi
  usecase.parameterset_tsv_file = 'tophat_parameterset.tsv'

  # set input dataset
  # mainly for CUI sushi
  usecase.dataset_tsv_file = 'tophat_dataset.tsv'

  # also possible to load a input dataset from Sushi DB
  #usecase.dataset_sushi_id = 1

  # run (submit to workflow_manager)
  usecase.run
  #usecase.test_run

end

