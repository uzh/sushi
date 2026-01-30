require 'open-uri'
require 'json'
require 'fileutils'

module NfCoreInfoFetcher
  NFCORE_PIPELINES_URL = "https://nf-co.re/pipelines.json"
  GITHUB_API_BASE = "https://api.github.com/repos/nf-core"
  SCHEMA_INPUT_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/master/assets/schema_input.json"
  NEXTFLOW_SCHEMA_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/%s/nextflow_schema.json"
  
  # Parameters that are auto-managed by SUSHI (excluded from user params)
  AUTO_MANAGED_PARAMS = %w[outdir input email publish_dir_mode]
  
  # Input type detection patterns based on schema_input.json
  INPUT_TYPE_PATTERNS = {
    'id_list' => {
      # Pipelines that take accession IDs
      'required' => ['id'],
      'description_patterns' => [/accession/i, /SRA|ENA|GEO|DDBJ/i]
    },
    'samplesheet' => {
      # Standard samplesheet with sample + fastq columns
      'required' => ['sample'],
      'properties' => ['fastq_1', 'fastq_2', 'fasta']
    },
    'fasta_list' => {
      # Pipelines that take FASTA files
      'required' => ['fasta'],
      'description_patterns' => [/fasta/i, /genome/i]
    }
  }
  
  # Cache directory relative to this file: ../../tmp/nfcore_cache
  CACHE_DIR = File.expand_path('../../tmp/nfcore_cache', __FILE__)
  CACHE_TTL = 86400  # 24 hours

  COLUMN_MAPPING = {
    'sample' => 'Name',
    'fastq_1' => 'Read1',
    'fastq_2' => 'Read2',
    'fasta' => 'FASTA',
    'strandedness' => 'Strandedness'
  }
  
  # Category mapping from nf-core topics
  CATEGORY_MAPPING = {
    'rna-seq' => 'Transcriptomics',
    'dna-seq' => 'Genomics',
    'chip-seq' => 'Epigenetics',
    'metagenomics' => 'Metagenomics',
    'quality-control' => 'QC',
    'sc-rna-seq' => 'SingleCell',
    'single-cell' => 'SingleCell'
  }
  
  def self.fetch_all(pipeline_name)
    {
      description: fetch_description(pipeline_name),
      category: fetch_category(pipeline_name),
      latest_version: fetch_latest_version(pipeline_name),
      required_columns: infer_required_columns(pipeline_name)
    }
  end

  # Fetch basic info from pipelines.json (fast, single API call)
  def self.fetch_basic_info(pipeline_name)
    pipelines = fetch_pipelines_json
    return nil unless pipelines && pipelines['remote_workflows']
    
    pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
    return nil unless pipeline
    
    # Map topic to category
    category = nil
    if pipeline['topics']
      pipeline['topics'].each do |topic|
        if CATEGORY_MAPPING[topic]
          category = CATEGORY_MAPPING[topic]
          break
        end
      end
      category ||= pipeline['topics'].first
    end
    
    {
      description: pipeline['description'] ? "#{pipeline['description']}\n<a href='#{pipeline['url']}'>nf-core/#{pipeline_name}</a>" : "nf-core/#{pipeline_name} pipeline",
      category: category,
      latest_version: pipeline['releases']&.first&.dig('tag_name') || 'master',
      required_columns: ['Name']  # Minimal - will be overridden by YAML config if needed
    }
  end

  # Get all pipelines with basic info (efficient bulk load)
  # Note: We use minimal required_columns ['Name'] for maximum compatibility.
  # nf-core pipelines accept --input samplesheet, and the actual required columns
  # depend on each pipeline. Users should ensure their dataset has the appropriate
  # columns for the selected pipeline.
  def self.fetch_all_pipelines_with_info
    puts "NfCoreInfoFetcher: fetch_all_pipelines_with_info called"
    $stdout.flush
    pipelines = fetch_pipelines_json
    puts "NfCoreInfoFetcher: pipelines.json loaded, remote_workflows: #{pipelines&.dig('remote_workflows')&.size || 'nil'}"
    $stdout.flush
    return {} unless pipelines && pipelines['remote_workflows']
    
    result = {}
    pipelines['remote_workflows'].each do |pipeline|
      name = pipeline['name']
      
      # Map topic to category
      category = nil
      if pipeline['topics']
        pipeline['topics'].each do |topic|
          if CATEGORY_MAPPING[topic]
            category = CATEGORY_MAPPING[topic]
            break
          end
        end
        category ||= pipeline['topics'].first
      end
      
      result[name] = {
        'nf_core_name' => name,
        'description' => pipeline['description'] ? "#{pipeline['description']}\n<a href='#{pipeline['url']}'>nf-core/#{name}</a>" : "nf-core/#{name} pipeline",
        'category' => category,
        'latest_version' => pipeline['releases']&.first&.dig('tag_name') || 'master',
        'required_columns' => ['Name']  # Minimal - only require Name column
      }
    end
    
    result
  end

  def self.fetch_all_pipelines
    pipelines = fetch_pipelines_json
    return [] unless pipelines && pipelines['remote_workflows']
    
    pipelines['remote_workflows'].map { |p| p['name'] }
  end
  
  def self.fetch_pipelines_json
    cache_path = File.join(CACHE_DIR, 'pipelines.json')
    fetch_with_cache(NFCORE_PIPELINES_URL, cache_path)
  end
  
  def self.fetch_description(pipeline_name)
    pipelines = fetch_pipelines_json
    return nil unless pipelines
    
    pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
    if pipeline
      <<~EOS
        #{pipeline['description']}
        <a href='#{pipeline['url']}'>nf-core/#{pipeline_name}</a>
      EOS
    else
      "nf-core/#{pipeline_name} pipeline"
    end
  rescue => e
    warn "Failed to fetch description for #{pipeline_name}: #{e.message}"
    "nf-core/#{pipeline_name} pipeline"
  end
  
  def self.fetch_category(pipeline_name)
    pipelines = fetch_pipelines_json
    return nil unless pipelines
    
    pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
    if pipeline && pipeline['topics']
      pipeline['topics'].each do |topic|
        return topic if CATEGORY_MAPPING.keys.include?(topic)
      end
      pipeline['topics'].first
    else
      nil
    end
  rescue => e
    warn "Failed to fetch category for #{pipeline_name}: #{e.message}"
    nil
  end
  
  def self.fetch_latest_version(pipeline_name)
    # First try to get from pipelines.json (faster)
    pipelines = fetch_pipelines_json
    if pipelines && pipelines['remote_workflows']
      pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
      if pipeline && pipeline['releases'] && !pipeline['releases'].empty?
        return pipeline['releases'].first['tag_name']
      end
    end
    
    # Fallback to GitHub API
    url = "#{GITHUB_API_BASE}/#{pipeline_name}/releases/latest"
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_release.json")
    
    data = fetch_with_cache(url, cache_path)
    if data
      data['tag_name']
    else
      "master"
    end
  rescue => e
    warn "Failed to fetch latest version for #{pipeline_name}: #{e.message}"
    "master"
  end
  
  # Note: This method is kept for potential future use but currently not called
  # during bulk loading. We use minimal required_columns for maximum compatibility.
  def self.infer_required_columns(pipeline_name)
    url = SCHEMA_INPUT_URL_TEMPLATE % pipeline_name
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_schema_input.json")
    
    schema = fetch_with_cache(url, cache_path)
    
    required_cols = ['Name'] # Minimal default fallback
    
    if schema && schema['required']
      mapped_cols = schema['required'].map { |col| COLUMN_MAPPING[col] || col }
      mapped_cols << 'Name' unless mapped_cols.include?('Name')
      
      if schema['properties']
        ['fastq_2', 'fasta', 'strandedness'].each do |col|
          if schema['properties'][col]
             mapped_cols << COLUMN_MAPPING[col] if COLUMN_MAPPING[col]
          end
        end
      end
      
      required_cols = mapped_cols.uniq.compact
    end
    
    required_cols
  rescue => e
    warn "Failed to fetch schema for #{pipeline_name}: #{e.message}"
    ['Name']
  end
  
  # Detect input type from schema_input.json
  # Returns: 'samplesheet', 'id_list', 'fasta_list', or 'unknown'
  def self.detect_input_type(pipeline_name)
    url = SCHEMA_INPUT_URL_TEMPLATE % pipeline_name
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_schema_input.json")
    
    schema = fetch_with_cache(url, cache_path)
    return 'unknown' unless schema
    
    required = schema['required'] || []
    properties = schema['properties'] || {}
    
    # Check for ID list type (fetchngs-like)
    if required.include?('id') || properties.key?('id')
      id_prop = properties['id'] || {}
      desc = id_prop['description'] || ''
      if desc =~ /accession|SRA|ENA|GEO|DDBJ/i
        return 'id_list'
      end
    end
    
    # Check for standard samplesheet type
    if required.include?('sample') || properties.key?('sample')
      if properties.key?('fastq_1') || properties.key?('fastq_2')
        return 'samplesheet'
      end
    end
    
    # Check for FASTA list type
    if required.include?('fasta') && !properties.key?('fastq_1')
      return 'fasta_list'
    end
    
    # Default to samplesheet for most pipelines
    'samplesheet'
  rescue => e
    warn "NfCoreInfoFetcher: Failed to detect input type for #{pipeline_name}: #{e.message}"
    'samplesheet'  # Default fallback
  end
  
  # Get full input configuration for a pipeline
  # Returns hash with input_type, samplesheet_columns, custom_params
  def self.get_input_config(pipeline_name)
    input_type = detect_input_type(pipeline_name)
    
    config = {
      'input_type' => input_type,
      'samplesheet_columns' => [],
      'custom_params' => []
    }
    
    case input_type
    when 'id_list'
      config['custom_params'] = [
        {
          'name' => 'accession_ids',
          'type' => 'text_area',
          'required' => true,
          'description' => 'Enter accession IDs (SRA/ENA/GEO/DDBJ), one per line'
        }
      ]
    when 'samplesheet'
      config['samplesheet_columns'] = infer_samplesheet_columns(pipeline_name)
    when 'fasta_list'
      config['custom_params'] = [
        {
          'name' => 'fasta_files',
          'type' => 'text_area',
          'required' => true,
          'description' => 'Enter FASTA file paths, one per line'
        }
      ]
    end
    
    config
  rescue => e
    warn "NfCoreInfoFetcher: Failed to get input config for #{pipeline_name}: #{e.message}"
    { 'input_type' => 'samplesheet', 'samplesheet_columns' => [], 'custom_params' => [] }
  end
  
  # Infer samplesheet columns from schema_input.json
  def self.infer_samplesheet_columns(pipeline_name)
    url = SCHEMA_INPUT_URL_TEMPLATE % pipeline_name
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_schema_input.json")
    
    schema = fetch_with_cache(url, cache_path)
    return [] unless schema && schema['properties']
    
    columns = []
    schema['properties'].each do |col_name, col_def|
      col_info = {
        'name' => col_name,
        'type' => col_def['type'] || 'string',
        'required' => (schema['required'] || []).include?(col_name),
        'description' => col_def['description'] || ''
      }
      col_info['format'] = col_def['format'] if col_def['format']
      col_info['pattern'] = col_def['pattern'] if col_def['pattern']
      col_info['enum'] = col_def['enum'] if col_def['enum']
      
      columns << col_info
    end
    
    columns
  rescue => e
    warn "NfCoreInfoFetcher: Failed to infer samplesheet columns for #{pipeline_name}: #{e.message}"
    []
  end
  
  # Fetch required parameters from nextflow_schema.json
  # Returns array of parameter hashes: [{name: 'param', type: 'string', description: '...', default: '...', enum: [...]}]
  def self.fetch_required_params(pipeline_name, version = 'master')
    url = NEXTFLOW_SCHEMA_URL_TEMPLATE % [pipeline_name, version]
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_#{version.gsub('.', '_')}_nextflow_schema.json")
    
    schema = fetch_with_cache(url, cache_path)
    return [] unless schema
    
    # Handle both old ('definitions') and new ('$defs') schema formats
    definitions = schema['definitions'] || schema['$defs']
    return [] unless definitions
    
    required_params = []
    
    definitions.each do |section_name, section|
      next unless section['properties']
      
      # Get required params for this section
      section_required = section['required'] || []
      
      section['properties'].each do |param_name, param_def|
        # Skip auto-managed and hidden params
        next if AUTO_MANAGED_PARAMS.include?(param_name)
        next if param_def['hidden'] == true
        
        # Only include if required in this section
        next unless section_required.include?(param_name)
        
        param_info = {
          'name' => param_name,
          'type' => param_def['type'] || 'string',
          'description' => param_def['description'] || '',
          'section' => section['title'] || section_name
        }
        
        # Add default if present
        param_info['default'] = param_def['default'] if param_def.key?('default')
        
        # Add enum options if present
        param_info['enum'] = param_def['enum'] if param_def['enum']
        
        # Add format info
        param_info['format'] = param_def['format'] if param_def['format']
        
        required_params << param_info
      end
    end
    
    required_params
  rescue => e
    warn "NfCoreInfoFetcher: Failed to fetch required params for #{pipeline_name}: #{e.message}"
    []
  end
  
  # Fetch all user-configurable parameters (not just required)
  # Useful for showing optional params that users might want to set
  def self.fetch_all_params(pipeline_name, version = 'master')
    url = NEXTFLOW_SCHEMA_URL_TEMPLATE % [pipeline_name, version]
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_#{version.gsub('.', '_')}_nextflow_schema.json")
    
    schema = fetch_with_cache(url, cache_path)
    return { 'required' => [], 'optional' => [] } unless schema
    
    # Handle both old ('definitions') and new ('$defs') schema formats
    definitions = schema['definitions'] || schema['$defs']
    return { 'required' => [], 'optional' => [] } unless definitions
    
    required_params = []
    optional_params = []
    
    definitions.each do |section_name, section|
      next unless section['properties']
      
      # Skip institutional/generic options sections (usually hidden)
      next if section_name =~ /institutional|generic/i
      
      section_required = section['required'] || []
      
      section['properties'].each do |param_name, param_def|
        # Skip auto-managed and hidden params
        next if AUTO_MANAGED_PARAMS.include?(param_name)
        next if param_def['hidden'] == true
        
        param_info = {
          'name' => param_name,
          'type' => param_def['type'] || 'string',
          'description' => param_def['description'] || '',
          'section' => section['title'] || section_name
        }
        
        param_info['default'] = param_def['default'] if param_def.key?('default')
        param_info['enum'] = param_def['enum'] if param_def['enum']
        param_info['format'] = param_def['format'] if param_def['format']
        
        if section_required.include?(param_name)
          required_params << param_info
        else
          optional_params << param_info
        end
      end
    end
    
    { 'required' => required_params, 'optional' => optional_params }
  rescue => e
    warn "NfCoreInfoFetcher: Failed to fetch all params for #{pipeline_name}: #{e.message}"
    { 'required' => [], 'optional' => [] }
  end
  
  def self.fetch_with_cache(url, cache_path)
    FileUtils.mkdir_p(File.dirname(cache_path))
    
    if File.exist?(cache_path) && File.size(cache_path) > 0 && (Time.now - File.mtime(cache_path)) < CACHE_TTL
      puts "NfCoreInfoFetcher: Using cache for #{url}"
      $stdout.flush
      return JSON.parse(File.read(cache_path, encoding: 'UTF-8'))
    end
    
    begin
      puts "NfCoreInfoFetcher: Fetching #{url}..."
      $stdout.flush
      data = URI.open(url, read_timeout: 30, open_timeout: 10).read
      # Force UTF-8 encoding for proper handling
      data = data.force_encoding('UTF-8')
      File.write(cache_path, data, encoding: 'UTF-8')
      puts "NfCoreInfoFetcher: Successfully fetched and cached #{url}"
      $stdout.flush
      JSON.parse(data)
    rescue OpenURI::HTTPError => e
      puts "NfCoreInfoFetcher: HTTP Error fetching #{url}: #{e.message}"
      $stdout.flush
      if File.exist?(cache_path) && File.size(cache_path) > 0
        puts "NfCoreInfoFetcher: Using stale cache for #{url}"
        $stdout.flush
        JSON.parse(File.read(cache_path, encoding: 'UTF-8'))
      else
        puts "NfCoreInfoFetcher: No valid cache available for #{url}"
        $stdout.flush
        nil
      end
    rescue => e
      puts "NfCoreInfoFetcher: Error fetching #{url}: #{e.message}"
      puts e.backtrace.first(5).join("\n")
      $stdout.flush
      if File.exist?(cache_path) && File.size(cache_path) > 0
        puts "NfCoreInfoFetcher: Using stale cache for #{url}"
        $stdout.flush
        begin
          JSON.parse(File.read(cache_path, encoding: 'UTF-8'))
        rescue => parse_err
          puts "NfCoreInfoFetcher: Cache file is invalid: #{parse_err.message}"
          $stdout.flush
          nil
        end
      else
        puts "NfCoreInfoFetcher: No valid cache available for #{url}"
        $stdout.flush
        nil
      end
    end
  end
end
