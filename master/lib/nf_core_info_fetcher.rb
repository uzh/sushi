require 'open-uri'
require 'json'
require 'fileutils'
require_relative 'nf_core_conventions'

module NfCoreInfoFetcher
  NFCORE_PIPELINES_URL = "https://nf-co.re/pipelines.json"
  GITHUB_API_BASE = "https://api.github.com/repos/nf-core"
  SCHEMA_INPUT_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/master/assets/schema_input.json"
  NEXTFLOW_SCHEMA_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/%s/nextflow_schema.json"
  NEXTFLOW_CONFIG_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/%s/nextflow.config"
  
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

  # Legacy column mapping - use NfCoreConventions::COLUMN_MAPPING instead
  COLUMN_MAPPING = {
    'sample' => 'Name',
    'fastq_1' => 'Read1',
    'fastq_2' => 'Read2',
    'fasta' => 'FASTA',
    'strandedness' => 'Strandedness'
  }
  
  # Legacy category mapping - use NfCoreConventions::CATEGORY_MAPPING instead
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
        'releases' => pipeline['releases'] || [],
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
  
  # ============================================================
  # Nextflow Version Compatibility
  # ============================================================
  
  # Cached installed Nextflow version (detected once per process)
  @@installed_nextflow_version = nil
  
  # Detect the installed Nextflow version
  # @return [String, nil] version string (e.g., "24.10.3") or nil if not found
  def self.installed_nextflow_version
    return @@installed_nextflow_version if @@installed_nextflow_version
    
    begin
      output = `nextflow -v 2>/dev/null`.strip
      # Parse "nextflow version 24.10.3.5933" -> "24.10.3"
      if output =~ /version\s+([\d]+\.[\d]+\.[\d]+)/
        @@installed_nextflow_version = $1
        puts "NfCoreInfoFetcher: Detected installed Nextflow version: #{@@installed_nextflow_version}"
        $stdout.flush
      else
        puts "NfCoreInfoFetcher: Could not parse Nextflow version from: #{output}"
        $stdout.flush
      end
    rescue => e
      warn "NfCoreInfoFetcher: Failed to detect Nextflow version: #{e.message}"
    end
    
    @@installed_nextflow_version
  end
  
  # Fetch the minimum required Nextflow version from a pipeline's nextflow.config
  # Parses the line: nextflowVersion = '!>=24.10.5'
  # @param pipeline_name [String] nf-core pipeline name
  # @param version [String] pipeline version tag
  # @return [String, nil] minimum required version (e.g., "24.10.5") or nil
  def self.fetch_required_nextflow_version(pipeline_name, version)
    url = NEXTFLOW_CONFIG_URL_TEMPLATE % [pipeline_name, version]
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_#{version.gsub('.', '_')}_nextflow.config")
    
    # Use text caching (not JSON)
    FileUtils.mkdir_p(File.dirname(cache_path))
    
    config_text = nil
    if File.exist?(cache_path) && File.size(cache_path) > 0 && (Time.now - File.mtime(cache_path)) < CACHE_TTL
      config_text = File.read(cache_path, encoding: 'UTF-8')
    else
      begin
        config_text = URI.open(url, read_timeout: 15, open_timeout: 10).read
        config_text = config_text.force_encoding('UTF-8')
        File.write(cache_path, config_text, encoding: 'UTF-8')
      rescue => e
        # Try stale cache
        if File.exist?(cache_path) && File.size(cache_path) > 0
          config_text = File.read(cache_path, encoding: 'UTF-8')
        else
          return nil
        end
      end
    end
    
    return nil unless config_text
    
    # Parse: nextflowVersion = '!>=24.10.5' or '>=24.10.5'
    if config_text =~ /nextflowVersion\s*=\s*['"]!?>?=?([\d]+\.[\d]+\.[\d]+)['"]/
      required = $1
      return required
    end
    
    nil
  end
  
  # Compare two version strings (e.g., "24.10.3" vs "24.10.5")
  # @return [Integer] -1, 0, or 1 (like <=>)
  def self.compare_versions(v1, v2)
    parts1 = v1.split('.').map(&:to_i)
    parts2 = v2.split('.').map(&:to_i)
    
    # Pad to same length
    max_len = [parts1.length, parts2.length].max
    parts1 += [0] * (max_len - parts1.length)
    parts2 += [0] * (max_len - parts2.length)
    
    parts1 <=> parts2
  end
  
  # Check if the installed Nextflow version meets the pipeline's requirement
  # @param required_version [String] minimum required version
  # @return [Boolean]
  def self.nextflow_version_compatible?(required_version)
    installed = installed_nextflow_version
    return true unless installed  # If we can't detect, assume compatible
    return true unless required_version  # If no requirement, assume compatible
    
    compare_versions(installed, required_version) >= 0
  end
  
  # Find the latest compatible version of a pipeline for the installed Nextflow
  # Iterates through releases from newest to oldest, checking compatibility
  # @param pipeline_name [String] nf-core pipeline name
  # @param releases [Array<Hash>] release info from pipelines.json
  # @return [Hash] { 'version' => '1.0.1', 'required_nf_version' => '24.04.2', 'is_fallback' => true/false }
  def self.find_compatible_version(pipeline_name, releases = nil)
    installed = installed_nextflow_version
    
    # Get releases if not provided
    unless releases
      pipelines = fetch_pipelines_json
      if pipelines && pipelines['remote_workflows']
        pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
        releases = pipeline['releases'] if pipeline
      end
    end
    
    return { 'version' => 'master', 'required_nf_version' => nil, 'is_fallback' => false } unless releases
    
    # Filter out 'dev' releases
    valid_releases = releases.select { |r| r['tag_name'] && r['tag_name'] =~ /^\d/ }
    
    return { 'version' => 'master', 'required_nf_version' => nil, 'is_fallback' => false } if valid_releases.empty?
    
    latest_version = valid_releases.first['tag_name']
    
    # If no Nextflow detected, return latest
    unless installed
      return { 'version' => latest_version, 'required_nf_version' => nil, 'is_fallback' => false }
    end
    
    # Check each release from newest to oldest
    valid_releases.each do |release|
      version = release['tag_name']
      required_nf = fetch_required_nextflow_version(pipeline_name, version)
      
      if required_nf.nil? || nextflow_version_compatible?(required_nf)
        is_fallback = (version != latest_version)
        if is_fallback
          puts "NfCoreInfoFetcher: Pipeline #{pipeline_name} latest (#{latest_version}) requires Nextflow >= #{fetch_required_nextflow_version(pipeline_name, latest_version)}, " \
               "but installed is #{installed}. Falling back to compatible version: #{version} (requires >= #{required_nf || 'unknown'})"
        else
          puts "NfCoreInfoFetcher: Pipeline #{pipeline_name} v#{version} is compatible with Nextflow #{installed} (requires >= #{required_nf || 'unknown'})"
        end
        $stdout.flush
        
        return {
          'version' => version,
          'required_nf_version' => required_nf,
          'is_fallback' => is_fallback
        }
      end
    end
    
    # No compatible version found - return latest anyway with warning
    puts "NfCoreInfoFetcher: WARNING: No compatible version found for #{pipeline_name} with Nextflow #{installed}. Using latest: #{latest_version}"
    $stdout.flush
    { 'version' => latest_version, 'required_nf_version' => nil, 'is_fallback' => false, 'incompatible' => true }
  end
  
  # ============================================================
  # Auto-detection Methods (Convention over Configuration)
  # ============================================================
  
  # Detect if pipeline needs reference selector based on nextflow_schema.json
  # Looks for fasta/gtf/genome parameters in reference-related sections
  # @param pipeline_name [String] nf-core pipeline name
  # @param version [String] pipeline version (default: 'master')
  # @return [Boolean]
  def self.needs_ref_selector?(pipeline_name, version = 'master')
    schema = fetch_nextflow_schema(pipeline_name, version)
    return false unless schema
    
    # Handle both old ('definitions') and new ('$defs') schema formats
    definitions = schema['definitions'] || schema['$defs']
    return false unless definitions
    
    definitions.each do |section_name, section|
      # Check if section name indicates reference genome options
      next unless NfCoreConventions.is_ref_section?(section_name) ||
                  section_name =~ /reference|genome/i
      
      properties = section['properties'] || {}
      # Check if any reference selector parameters exist
      if NfCoreConventions::REF_SELECTOR_PARAMS.any? { |p| properties.key?(p) }
        puts "NfCoreInfoFetcher: Pipeline #{pipeline_name} needs ref_selector (found in section: #{section_name})"
        $stdout.flush
        return true
      end
    end
    
    puts "NfCoreInfoFetcher: Pipeline #{pipeline_name} does not need ref_selector"
    $stdout.flush
    false
  rescue => e
    warn "NfCoreInfoFetcher: Failed to check ref_selector for #{pipeline_name}: #{e.message}"
    false
  end
  
  # Auto-generate samplesheet mapping from schema_input.json using conventions
  # @param pipeline_name [String] nf-core pipeline name
  # @return [Hash] mapping from nf-core column to SUSHI column
  def self.auto_samplesheet_mapping(pipeline_name)
    schema = fetch_schema_input(pipeline_name)
    return {} unless schema
    
    # schema_input.json can have properties at top level or under items.properties
    properties = schema['properties'] || (schema['items'] && schema['items']['properties'])
    return {} unless properties
    
    mapping = {}
    properties.each do |nf_col, _|
      sushi_col = NfCoreConventions.sushi_column_for(nf_col)
      mapping[nf_col] = sushi_col if sushi_col
    end
    
    puts "NfCoreInfoFetcher: Auto-generated samplesheet mapping for #{pipeline_name}: #{mapping.inspect}"
    $stdout.flush
    mapping
  rescue => e
    warn "NfCoreInfoFetcher: Failed to auto-generate samplesheet mapping for #{pipeline_name}: #{e.message}"
    {}
  end
  
  # Get column defaults for nf-core columns (e.g., strandedness -> 'auto')
  # @param pipeline_name [String] nf-core pipeline name
  # @return [Hash] column name => default value
  def self.get_column_defaults(pipeline_name)
    schema = fetch_schema_input(pipeline_name)
    return {} unless schema
    
    # schema_input.json can have properties at top level or under items.properties
    properties = schema['properties'] || (schema['items'] && schema['items']['properties'])
    return {} unless properties
    
    defaults = {}
    properties.each do |nf_col, _|
      default_val = NfCoreConventions.default_for(nf_col)
      defaults[nf_col] = default_val if default_val
    end
    
    puts "NfCoreInfoFetcher: Column defaults for #{pipeline_name}: #{defaults.inspect}"
    $stdout.flush
    defaults
  rescue => e
    warn "NfCoreInfoFetcher: Failed to get column defaults for #{pipeline_name}: #{e.message}"
    {}
  end
  
  # Fetch schema_input.json for a pipeline
  # @param pipeline_name [String] nf-core pipeline name
  # @return [Hash, nil] parsed JSON or nil on failure
  def self.fetch_schema_input(pipeline_name)
    url = SCHEMA_INPUT_URL_TEMPLATE % pipeline_name
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_schema_input.json")
    fetch_with_cache(url, cache_path)
  end
  
  # Fetch nextflow_schema.json for a pipeline
  # @param pipeline_name [String] nf-core pipeline name
  # @param version [String] pipeline version
  # @return [Hash, nil] parsed JSON or nil on failure
  def self.fetch_nextflow_schema(pipeline_name, version = 'master')
    url = NEXTFLOW_SCHEMA_URL_TEMPLATE % [pipeline_name, version]
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_#{version.gsub('.', '_')}_nextflow_schema.json")
    fetch_with_cache(url, cache_path)
  end
  
  # Get resource defaults based on category from conventions
  # @param category [String] SUSHI analysis category
  # @return [Hash] resource defaults (cores, ram, scratch)
  def self.resource_defaults_for_category(category)
    NfCoreConventions.resources_for(category)
  end
  
  # Map nf-core topic to SUSHI category using conventions
  # @param topic [String] nf-core topic
  # @return [String, nil] SUSHI category
  def self.category_from_topic(topic)
    # First try conventions, then legacy mapping
    NfCoreConventions.category_for(topic) || CATEGORY_MAPPING[topic]
  end
  
  # Infer SUSHI category from pipeline topics
  # @param pipeline_name [String] nf-core pipeline name
  # @return [String, nil] SUSHI category
  def self.infer_category(pipeline_name)
    pipelines = fetch_pipelines_json
    return nil unless pipelines && pipelines['remote_workflows']
    
    pipeline = pipelines['remote_workflows'].find { |p| p['name'] == pipeline_name }
    return nil unless pipeline && pipeline['topics']
    
    # Find first matching topic
    pipeline['topics'].each do |topic|
      category = category_from_topic(topic)
      return category if category
    end
    
    # Return first topic if no mapping found
    pipeline['topics'].first
  rescue => e
    warn "NfCoreInfoFetcher: Failed to infer category for #{pipeline_name}: #{e.message}"
    nil
  end
end
