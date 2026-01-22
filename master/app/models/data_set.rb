class DataSet < ActiveRecord::Base
#  attr_accessible :name, :md5, :comment
  has_many :samples
  has_many :jobs, :foreign_key => :next_dataset_id
  belongs_to :project
  has_many :data_sets, :foreign_key => :parent_id
  belongs_to :data_set, :foreign_key => :parent_id
  serialize :runnable_apps, Hash
  serialize :nfcore_apps, Array
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
end
