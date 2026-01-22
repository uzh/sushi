class DataSet < ActiveRecord::Base
#  attr_accessible :name, :md5, :comment
  has_many :samples
  has_many :jobs, :foreign_key => :next_dataset_id
  belongs_to :project
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id
  serialize :runnable_apps, Hash
  belongs_to :user
  serialize :order_ids, Array
  serialize :job_parameters, Hash

  def headers
    self.samples.map{|sample| sample.to_hash.keys}.flatten.uniq
  end
  def factor_first_headers
    headers.sort_by do |col| 
      if col == 'Name' 
        0
      elsif col.scan(/\[(.*)\]/).flatten.join =~ /Factor/
        1
      else
        2
      end
    end
  end
  def saved?
    if DataSet.find_by_md5(md5hexdigest)
      true
    else
      false
    end
  end
  def md5hexdigest
    key_value = self.samples.map{|sample| sample.key_value}.join + self.parent_id.to_s + self.project_id.to_s
    Digest::MD5.hexdigest(key_value)
  end
  def tsv_string
    string = CSV.generate(:col_sep=>"\t") do |out|
      out << headers
      self.samples.each do |sample|
        out << headers.map{|header| 
          val = sample.to_hash[header]
          val.to_s.empty? ? nil:val}
      end
    end
    string
  end
  
  # Executes external command and returns [stdout, stderr, exit_status]
  def run_external_command(command)
    r_out, w_out = IO.pipe
    r_err, w_err = IO.pipe
    pid = Process.spawn(command, out: w_out, err: w_err)
    w_out.close
    w_err.close

    stdout_buf = +""
    stderr_buf = +""
    t_out = Thread.new { r_out.each_line { |l| stdout_buf << l } }
    t_err = Thread.new { r_err.each_line { |l| stderr_buf << l } }

    _, status = Process.wait2(pid)
    t_out.join
    t_err.join

    [stdout_buf, stderr_buf, status.exitstatus]
  end
  private :run_external_command
  
  # Returns true if the dataset has at least one column tagged with "[File]"
  def has_file_tag_column?
    headers.any? { |header| header && header.tag?('File') }
  end
  private :has_file_tag_column?

  # Detect duplicate column headers when ignoring tag suffixes like "[File]" or "[Link]".
  # Returns an array of base header names that appear more than once after tag removal.
  def duplicate_headers_ignoring_tags
    normalized_headers = headers.map { |header| header.to_s.gsub(/\s*\[.+\]/, '').strip }
    normalized_headers.tally.select { |_, count| count > 1 }.keys
  end
  private :duplicate_headers_ignoring_tags
  def samples_length
    unless self.num_samples
      self.num_samples=self.samples.length
      self.save
    end
    self.num_samples
  end
  def check_order_ids
    order_ids_ = {}
    self.samples.each do |sample_|
      sample = sample_.to_hash
      if order_ids__ = sample["Order Id [B-Fabric]"]
        order_ids__.strip!
        order_ids__.split(",").each do |order_id|
          order_id.strip!
          order_ids_[order_id] = true
        end
      end
    end
    unless order_ids_.empty?
      self.order_ids.concat(order_ids_.keys)
      self.order_ids.uniq!
      self.order_ids.sort!

      # OrderID save
      if self.order_ids.uniq.length == 1 and order_id = self.order_ids.first.to_i
        self.order_id = order_id
      end

      self.save
    end
  end
  def register_bfabric(op = 'new', bfabric_application_number: nil, register_child_dataset_too: nil, update_completed_samples: false)
    register_command = "register_sushi_dataset_into_bfabric"
    check_command = "check_dataset_bfabric"
    parent_dataset = self.data_set
    if parent_dataset.nil? or parent_dataset.bfabric_id
      if SushiFabric::Application.config.fgcz? and system("which #{register_command} > /dev/null 2>&1") and system("which #{check_command} > /dev/null 2>&1")
        time = Time.new.strftime("%Y%m%d-%H%M%S")
        dataset_tsv = File.join(SushiFabric::Application.config.scratch_dir, "dataset.#{self.id}_#{time}.tsv")
        # Abort early if there are duplicate headers ignoring tag suffixes (e.g., "[File]", "[Link]")
        if (duplicate_columns = duplicate_headers_ignoring_tags).any?
          warn "# Not run DataSet#register_bfabric for DataSet(id=#{self.id}, name=\"#{self.name}\") due to duplicate column names (ignoring tags): #{duplicate_columns.join(', ')}"
          return false
        end
        # Abort if there is no [File]-tagged column in the dataset
        unless has_file_tag_column?
          warn "# Not run DataSet#register_bfabric for DataSet(id=#{self.id}, name=\"#{self.name}\") due to no [File] tag column in dataset"
          return false
        end
        option_check = if ((op == 'new' or op == 'only_one') and !self.bfabric_id) or op == 'renewal'
                         true
                       elsif op == 'update' and bfabric_id = self.bfabric_id
                         com = "#{check_command} #{bfabric_id}"
                         warn "$ #{com}"
                         if out = `#{com}`
                           warn "# returned: #{out.chomp.downcase}"
                           eval(out.chomp.downcase)
                         end
                       end

        self.check_order_ids

        warn "parent_dataset.nil?= #{parent_dataset.nil?}"
        warn "self.order_ids= #{self.order_ids}"
        command = if self.order_ids.uniq.length == 1 and order_id = self.order_ids.first.to_i
                    if parent_dataset.nil? # root dataset
                      if order_id > 8000
                        [register_command, "o#{self.order_ids.first}", dataset_tsv, self.name, self.id, "--skip-file-check"].join(" ")
                      else
                        [register_command, "p#{self.project.number}", dataset_tsv, self.name, self.id, "--skip-file-check"].join(" ")
                      end
                    elsif parent_dataset and bfabric_id = parent_dataset.bfabric_id and register_child_dataset_too # child dataset
                      if order_id > 8000
                        [register_command, "o#{self.order_ids.first}", dataset_tsv, self.name, self.id, bfabric_id, "--sushi-app #{self.sushi_app_name} --skip-file-check"].join(" ")
                      else
                        [register_command, "p#{self.project.number}", dataset_tsv, self.name, self.id, bfabric_id, "--sushi-app #{self.sushi_app_name} --skip-file-check"].join(" ")
                      end
                    end
                  elsif self.order_ids.uniq.length > 1 # multi order dataset
                    if parent_dataset.nil? # root dataset
                        [register_command, "p#{self.project.number}", dataset_tsv, self.name, self.id, "--skip-file-check"].join(" ")
                    elsif parent_dataset and bfabric_id = parent_dataset.bfabric_id and register_child_dataset_too # child dataset
                        [register_command, "p#{self.project.number}", dataset_tsv, self.name, self.id, bfabric_id, "--sushi-app #{self.sushi_app_name} --skip-file-check"].join(" ")
                    end
                  end

        if command and bfabric_application_number
          command << " -a #{bfabric_application_number}"
        end
        if command and option_check
          open(dataset_tsv, "w") do |out|
            out.print self.tsv_string
          end
          warn "# created: #{dataset_tsv}"
          if File.exist?(dataset_tsv)
          bfabric_ids, err, exit_status = run_external_command(command)
          warn "$ #{command}"
          warn "# mode: #{op}"
          warn "# bfabric_ids: #{bfabric_ids}"
          # Treat as success if an IDs line ("<int>,<int>") exists, even when exit status is non-zero.
          ids_line = bfabric_ids.lines.find { |l| l.strip.match?(/^\d+\s*,\s*\d+$/) }
          if exit_status != 0 && !ids_line
            warn "# register_sushi_dataset_into_bfabric exited with status #{exit_status}; stderr: #{err}"
            File.unlink dataset_tsv
            warn "# removed: #{dataset_tsv}"
            return false
          end
          parsed_line = (ids_line ? ids_line.strip : bfabric_ids.chomp)
          if parsed_line.split(/\n/).uniq.length < 2
            workunit_id, dataset_id = parsed_line.split(',').map(&:strip)
              if workunit_id.to_i > 0
                self.bfabric_id = dataset_id.to_i
                self.workunit_id = workunit_id.to_i
                warn "# DataSetID  (BFabric): #{self.bfabric_id}"
                warn "# WorkunitID (BFabric): #{self.workunit_id}"
                self.save
              end
            else
              warn "# Not executed properly:"
              warn "# BFabricID: #{bfabric_id}"
            end
            File.unlink dataset_tsv
            warn "# removed: #{dataset_tsv}"
          end
        else
          warn "# Not run DataSet#register_bfabric due to OrderID is missing in DataSet"
        end

        # Update completed_samples if requested
        if update_completed_samples
          sample_available = 0
          if self.completed_samples.to_i != self.samples_length
            self.samples.each do |sample|
              if sample_file = sample.to_hash.select{|header, file| header and header.tag?('File')}.first
                file_list = sample_file.last.split(",") ## sample_file is an array holding the header and the file
                all_files_exist = file_list.all? { |f| File.exist?(File.join(SushiFabric::GSTORE_DIR, f)) }
                if all_files_exist
                  sample_available+=1
                end
              else # in case of no [File] tag sample
                sample_available+=1
              end
            end
          else
            sample_available = self.samples_length
          end
          self.completed_samples = sample_available
          self.save
        end

        unless op == 'only_one'
          if child_data_sets = self.data_sets
            child_data_sets.each do |child_data_set|
              child_data_set.register_bfabric(op, update_completed_samples: update_completed_samples)
            end
          end
        end
      end
    else
      warn "# Not run DataSet#register_bfabric because its parental dataset is not registered in bfabric"
    end
  end
  def update_resource_size
    #command = "update_resource_size -w #{self.workunit_id} --test"
    command = "update_resource_size -w #{self.workunit_id}"
    #p command
    `#{command}`
  end
  def paths
    dirs = []
    samples.each do |sample|
      sample.to_hash.each do |header, file|
        if (header.tag?('File') or header.tag?('Link')) and file !~ /http/
          dirs << File.dirname(file)
        end
      end
    end
    dirs = dirs.map{|path| path.split('/')[0,2].join('/')}.uniq
  end
  def save_as_tsv(out_file=nil)
    file_name = out_file||"#{self.name}_dataset.tsv"
    File.write(file_name, tsv_string)
  end
  def sample_paths
    paths = []
    self.samples.each do |sample|
      begin
        sample.to_hash.each do |header, file|
          if (header.tag?('File') or header.tag?('Link')) and file !~ /http/ and !file.nil?
            paths << File.dirname(file)
          end
        end
      rescue => e
        warn "Error in sample #{sample}: #{e.message}"
      end
    end
    paths.uniq!
    paths.map{|path| path.split('/')[0,2].join('/')}
  end

  def self.save_dataset_to_database(data_set_arr:, headers:, rows:, user: nil, child: nil, sushi_app_name: nil)
    data_set_hash = Hash[*data_set_arr]
    
    # Find or create project
    project = Project.find_by_number(data_set_hash['ProjectNumber'].to_i)
    unless project
      project = Project.new
      project.number = data_set_hash['ProjectNumber'].to_i
      project.save
    end

    # Create new dataset
    data_set = DataSet.new
    data_set.user = user if user
    data_set.name = data_set_hash['DataSetName']
    data_set.project = project
    data_set.child = child if child

    # Set parent dataset if exists
    if parent_id = data_set_hash['ParentID']
      parent_data_set = DataSet.find_by_id(parent_id.to_i)
      data_set.data_set = parent_data_set if parent_data_set
    end

    # Set comment if exists
    if comment = data_set_hash['Comment'] and !comment.to_s.empty?
      data_set.comment = comment
    end

    # Set sushi_app_name if exists
    data_set.sushi_app_name = sushi_app_name if sushi_app_name

    # Create samples
    sample_hash = {}
    rows.each do |row|
      headers.each_with_index do |header, i|
        sample_hash[header] = row[i]
      end
      sample = Sample.new
      sample.key_value = sample_hash.to_s
      sample.save
      data_set.samples << sample
    end

    data_set.save
    data_set.id
  end
  def file_paths
    samples.flat_map do |sample|
      begin
        sample.to_hash.select{|header, file| header and header.tag?('File')}.values
      rescue => e
        Rails.logger.warn "Failed to parse key_value for Sample ID=#{sample.id}: #{e}"
        []
      end
    end.uniq
  end
  
  # Instance method: Merge this dataset with another dataset
  # Creates a new child dataset with concatenated Read files
  # @param other_dataset [DataSet] The dataset to merge with
  # @param options [Hash] Merge options
  # @return [DataSet] The newly created merged dataset
  def merge_with(other_dataset, options: {})
    require 'set'
    
    raise "Other dataset is required" unless other_dataset
    raise "Other dataset must be a DataSet object" unless other_dataset.is_a?(DataSet)
    
    # Parse options
    matching_column = options[:matching_column] || 'Name'
    excluded_columns = options[:excluded_columns] || ['Sample Id [B-Fabric]']
    merged_dataset_name = options[:merged_dataset_name] || "MergedReads_#{self.id}_#{other_dataset.id}"
    user = options[:user]
    parent_id = options[:parent_id] || self.id
    
    # Build hash of samples from other dataset
    dataset_hash2 = {}
    other_dataset.samples.each do |sample|
      sample_hash = sample.to_hash
      key = sample_hash[matching_column]
      dataset_hash2[key] = sample_hash if key
    end
    
    # Process this dataset's samples and merge with other dataset
    merged_samples = []
    processed_samples = Set.new
    all_headers = Set.new
    
    self.samples.each do |sample|
      sample_hash = sample.to_hash.clone
      key = sample_hash[matching_column]
      processed_samples.add(key)
      
      # Collect all headers
      sample_hash.keys.each { |h| all_headers.add(h) }
      
      # Merge with dataset2 if match exists
      if key && dataset_hash2[key]
        sample2 = dataset_hash2[key]
        
        # Merge all columns from dataset2
        sample2.keys.each do |col|
          # Skip matching column (already same)
          next if col == matching_column
          
          # Handle Read1/Read2 specially (concatenate with comma)
          if col =~ /^Read[12]\s*\[File\]/
            if sample_hash[col] && sample2[col]
              sample_hash[col] = "#{sample_hash[col]},#{sample2[col]}"
            elsif sample2[col]
              sample_hash[col] = sample2[col]
            end
          # Handle Read Count (sum)
          elsif col == 'Read Count'
            if sample_hash[col] && sample2[col]
              count1 = sample_hash[col].to_i
              count2 = sample2[col].to_i
              sample_hash[col] = (count1 + count2).to_s
            elsif sample2[col]
              sample_hash[col] = sample2[col]
            end
          # For other columns, take value from dataset2 if dataset1 doesn't have it
          else
            sample_hash[col] = sample2[col] unless sample_hash[col]
          end
        end
        
        # Collect headers from sample2
        sample2.keys.each { |h| all_headers.add(h) }
      end
      
      # Remove excluded columns
      excluded_columns.each { |col| sample_hash.delete(col) }
      
      merged_samples << sample_hash
    end
    
    # Add samples that exist only in dataset2
    dataset_hash2.each do |key, sample_data|
      unless processed_samples.include?(key)
        sample_hash = sample_data.clone
        
        # Collect headers
        sample_hash.keys.each { |h| all_headers.add(h) }
        
        # Remove excluded columns
        excluded_columns.each { |col| sample_hash.delete(col) }
        
        merged_samples << sample_hash
      end
    end
    
    # Sort by matching column
    merged_samples.sort_by! { |row| row[matching_column] || '' }
    
    # Remove excluded columns from headers
    excluded_columns.each { |col| all_headers.delete(col) }
    headers = all_headers.to_a
    
    # Build rows array for save_dataset_to_database
    rows = []
    merged_samples.each do |sample_hash|
      row = headers.map { |header| sample_hash[header] }
      rows << row
    end
    
    # Prepare dataset array
    data_set_arr = [
      'DataSetName', merged_dataset_name,
      'ProjectNumber', self.project.number,
      'ParentID', parent_id,
      'Comment', "Merged datasets: #{self.name} + #{other_dataset.name}"
    ]
    
    # Save to database
    merged_dataset_id = DataSet.save_dataset_to_database(
      data_set_arr: data_set_arr,
      headers: headers,
      rows: rows,
      user: user,
      child: true,
      sushi_app_name: 'MergeReadDatasets'
    )
    
    # Return the new dataset
    DataSet.find_by_id(merged_dataset_id)
  end
  
  # Class method: Merge two datasets by ID (for backward compatibility)
  # @param dataset1_id [Integer] ID of the first dataset
  # @param dataset2_id [Integer] ID of the second dataset
  # @param options [Hash] Merge options
  # @return [DataSet] The newly created merged dataset
  def self.merge_datasets(dataset1_id:, dataset2_id:, options: {})
    dataset1 = DataSet.find_by_id(dataset1_id.to_i)
    dataset2 = DataSet.find_by_id(dataset2_id.to_i)
    
    raise "Dataset1 (ID: #{dataset1_id}) not found" unless dataset1
    raise "Dataset2 (ID: #{dataset2_id}) not found" unless dataset2
    
    # Call instance method
    dataset1.merge_with(dataset2, options: options)
  end

  # Generate Materials & Methods document using LLM
  # @param use_llm [Boolean] Whether to use LLM for generation (default: true)
  # @return [String] M&M content in Markdown
  def generate_materials_and_methods(use_llm: true)
    require Rails.root.join('lib', 'llm_client')
    
    # Collect all analysis data
    analysis_data = collect_analysis_data
    
    if use_llm
      begin
        # Use LLM for generating natural language M&M
        LlmClient.generate_mm(analysis_data)
      rescue => e
        Rails.logger.error("LLM M&M generation failed: #{e.message}")
        # Fallback to template-based generation
        generate_template_mm(analysis_data)
      end
    else
      generate_template_mm(analysis_data)
    end
  end

  # Collect all analysis data from DataSet and gStore
  # @return [Hash] Analysis data for M&M generation
  def collect_analysis_data
    data = {
      name: self.name,
      id: self.id,
      project: self.project ? "p#{self.project.number}" : nil,
      date: self.created_at&.strftime('%Y-%m-%d'),
      sushi_app: self.sushi_app_name,
      sample_count: self.samples.length,
      bfabric_id: self.bfabric_id,
      workunit_id: self.workunit_id,
      order_id: self.order_id
    }
    
    # Filter job parameters (exclude internal ones)
    if self.job_parameters.present?
      internal_params = ['cores', 'ram', 'scratch', 'node', 'process_mode', 'mail', 'partition']
      data[:parameters] = self.job_parameters.reject { |k, v| internal_params.include?(k) || v.to_s.empty? }
      data[:ref_build] = self.job_parameters['refBuild']
      data[:ref_feature] = self.job_parameters['refFeatureFile']
    end
    
    # Get gStore path
    if paths = self.sample_paths and !paths.empty?
      gstore_path = File.join(SushiFabric::GSTORE_DIR, paths.first)
      data[:gstore_path] = gstore_path
      
      # Read job script if available
      data[:job_script] = read_job_script(gstore_path)
      
      # Read input samples from input_dataset.tsv
      data[:input_samples] = read_input_dataset(gstore_path)
      
      # Read output files from dataset.tsv
      data[:output_files] = read_output_dataset(gstore_path)
      
      # Read parameters.tsv for additional info
      data[:parameters_file] = read_parameters_file(gstore_path)
    end
    
    data
  end

  # Read job script content from gStore
  # @param gstore_path [String] Path to dataset directory in gStore
  # @return [String, nil] Job script content or nil
  def read_job_script(gstore_path)
    scripts_dir = File.join(gstore_path, 'scripts')
    return nil unless File.directory?(scripts_dir)
    
    # Find shell scripts
    script_files = Dir.glob(File.join(scripts_dir, '*.sh'))
    return nil if script_files.empty?
    
    # Read the first (or main) script
    # Prefer scripts that are not just wrappers
    main_script = script_files.find { |f| File.basename(f) !~ /^wrapper/ } || script_files.first
    
    begin
      content = File.read(main_script)
      # Limit size to avoid overwhelming the LLM
      content.length > 10000 ? content[0, 10000] + "\n... (truncated)" : content
    rescue => e
      Rails.logger.warn("Failed to read job script: #{e.message}")
      nil
    end
  end
  private :read_job_script

  # Read input sample names from input_dataset.tsv
  # @param gstore_path [String] Path to dataset directory in gStore
  # @return [Array<String>] List of sample names
  def read_input_dataset(gstore_path)
    input_file = File.join(gstore_path, 'input_dataset.tsv')
    return [] unless File.exist?(input_file)
    
    begin
      samples = []
      lines = File.readlines(input_file)
      return [] if lines.empty?
      
      headers = lines.first.chomp.split("\t")
      name_idx = headers.index('Name') || 0
      
      lines[1..-1].each do |line|
        cols = line.chomp.split("\t")
        samples << cols[name_idx] if cols[name_idx]
      end
      
      # Limit to first 20 samples
      samples.length > 20 ? samples[0, 20] + ["... (#{samples.length - 20} more)"] : samples
    rescue => e
      Rails.logger.warn("Failed to read input_dataset.tsv: #{e.message}")
      []
    end
  end
  private :read_input_dataset

  # Read output file information from dataset.tsv
  # @param gstore_path [String] Path to dataset directory in gStore
  # @return [Array<String>] List of output file descriptions
  def read_output_dataset(gstore_path)
    output_file = File.join(gstore_path, 'dataset.tsv')
    return [] unless File.exist?(output_file)
    
    begin
      lines = File.readlines(output_file)
      return [] if lines.empty?
      
      headers = lines.first.chomp.split("\t")
      # Find columns with [File] tag
      file_columns = headers.select { |h| h.include?('[File]') }
      
      file_columns.map { |col| col.gsub(/\s*\[File\]/, '') }
    rescue => e
      Rails.logger.warn("Failed to read dataset.tsv: #{e.message}")
      []
    end
  end
  private :read_output_dataset

  # Read parameters.tsv for additional analysis info
  # @param gstore_path [String] Path to dataset directory in gStore
  # @return [Hash, nil] Parameters hash or nil
  def read_parameters_file(gstore_path)
    params_file = File.join(gstore_path, 'parameters.tsv')
    return nil unless File.exist?(params_file)
    
    begin
      params = {}
      File.readlines(params_file).each do |line|
        next if line.strip.empty?
        key, value = line.chomp.split("\t", 2)
        params[key] = value if key && value
      end
      params
    rescue => e
      Rails.logger.warn("Failed to read parameters.tsv: #{e.message}")
      nil
    end
  end
  private :read_parameters_file

  # Generate template-based M&M (fallback when LLM unavailable)
  # @param data [Hash] Analysis data
  # @return [String] M&M content in Markdown
  def generate_template_mm(data)
    mm_content = []
    mm_content << "# Materials and Methods"
    mm_content << ""
    mm_content << "## Analysis Information"
    mm_content << ""
    mm_content << "- **Dataset Name**: #{data[:name]}"
    mm_content << "- **Dataset ID**: #{data[:id]}"
    mm_content << "- **Project**: #{data[:project]}" if data[:project]
    mm_content << "- **Analysis Date**: #{data[:date]}" if data[:date]
    mm_content << "- **SUSHI Application**: #{data[:sushi_app]}" if data[:sushi_app]
    mm_content << "- **Number of Samples**: #{data[:sample_count]}"
    mm_content << ""

    if data[:bfabric_id]
      mm_content << "## B-Fabric Registration"
      mm_content << ""
      mm_content << "- **B-Fabric Dataset ID**: #{data[:bfabric_id]}"
      mm_content << "- **B-Fabric Workunit ID**: #{data[:workunit_id]}" if data[:workunit_id]
      mm_content << "- **Order ID**: #{data[:order_id]}" if data[:order_id]
      mm_content << ""
    end

    if data[:parameters] && !data[:parameters].empty?
      mm_content << "## Analysis Parameters"
      mm_content << ""
      data[:parameters].each do |key, value|
        mm_content << "- **#{key}**: #{value}"
      end
      mm_content << ""
    end

    if data[:ref_build] || data[:ref_feature]
      mm_content << "## Reference Information"
      mm_content << ""
      mm_content << "- **Reference Build**: #{data[:ref_build]}" if data[:ref_build]
      mm_content << "- **Reference Feature File**: #{data[:ref_feature]}" if data[:ref_feature]
      mm_content << ""
    end

    mm_content << "## Software and Tools"
    mm_content << ""
    mm_content << "- **Analysis Platform**: SUSHI (https://fgcz-sushi.uzh.ch)"
    mm_content << "- **Execution Backend**: ezRun (R library)"
    mm_content << ""

    if data[:gstore_path]
      mm_content << "## Data Location"
      mm_content << ""
      mm_content << "- **gStore Path**: #{data[:gstore_path]}"
      mm_content << ""
    end

    mm_content << "---"
    mm_content << ""
    mm_content << "*This document was automatically generated by SUSHI on #{Time.now.strftime('%Y-%m-%d %H:%M:%S')}*"
    mm_content << "*Note: LLM generation was unavailable. This is a template-based output.*"

    mm_content.join("\n")
  end
  private :generate_template_mm

  # Push M&M file to GitLab repository
  # @param file_path [String] Path to the M&M file in gstore
  # @param commit_message [String] Git commit message
  # @return [Hash] Result with success status and message
  def push_to_gitlab(file_path, commit_message = nil)
    project_dir = File.join(SushiFabric::GSTORE_DIR, "p#{self.project.number}")
    git_dir = File.join(project_dir, '.git')

    unless File.exist?(git_dir)
      return { success: false, message: "Git repository not found at #{project_dir}" }
    end

    commit_message ||= "Add M&M for #{self.name} (DataSet ID: #{self.id})"

    begin
      Dir.chdir(project_dir) do
        # Get relative path for git add
        relative_path = file_path.sub("#{project_dir}/", '')
        
        # Git operations
        system("git add #{relative_path}")
        system("git commit -m '#{commit_message}'")
        push_result = system("git push origin master 2>&1")
        
        if push_result
          return { success: true, message: "Successfully pushed to GitLab" }
        else
          return { success: false, message: "Git push failed" }
        end
      end
    rescue => e
      return { success: false, message: "Git operation failed: #{e.message}" }
    end
  end
end
