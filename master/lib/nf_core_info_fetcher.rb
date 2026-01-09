require 'open-uri'
require 'json'
require 'fileutils'

module NfCoreInfoFetcher
  NFCORE_PIPELINES_URL = "https://nf-co.re/pipelines.json"
  GITHUB_API_BASE = "https://api.github.com/repos/nf-core"
  SCHEMA_URL_TEMPLATE = "https://raw.githubusercontent.com/nf-core/%s/master/nextflow_schema.json"
  
  # Cache directory relative to this file: ../../tmp/nfcore_cache
  CACHE_DIR = File.expand_path('../../tmp/nfcore_cache', __FILE__)
  CACHE_TTL = 86400  # 24 hours
  
  def self.fetch_all(pipeline_name)
    {
      description: fetch_description(pipeline_name),
      category: fetch_category(pipeline_name),
      latest_version: fetch_latest_version(pipeline_name),
      required_columns: infer_required_columns(pipeline_name)
    }
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
    # Map topics to categories
    # This is a best-effort mapping since nf-core topics are free-form
    if pipeline && pipeline['topics']
      pipeline['topics'].each do |topic|
        return topic if ['rna-seq', 'dna-seq', 'chip-seq', 'metagenomics', 'sc-rna-seq', 'single-cell'].include?(topic)
      end
      # Return the first topic if no standard one found
      pipeline['topics'].first
    else
      nil
    end
  rescue => e
    warn "Failed to fetch category for #{pipeline_name}: #{e.message}"
    nil
  end
  
  def self.fetch_latest_version(pipeline_name)
    # Use GitHub API to get latest release
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
  
  def self.infer_required_columns(pipeline_name)
    url = SCHEMA_URL_TEMPLATE % pipeline_name
    cache_path = File.join(CACHE_DIR, "#{pipeline_name}_schema.json")
    
    schema = fetch_with_cache(url, cache_path)
    
    # Default columns
    required_cols = ['Name', 'Read1']
    
    if schema && schema['definitions'] && schema['definitions']['input_output_options']
      input_prop = schema['definitions']['input_output_options']['properties']['input']
      # If input is a samplesheet (common in nf-core), we assume Name/Read1/Read2 structure
      # But strictly speaking we can't fully infer SUSHI columns just from "input: string"
      # So we stick to a sensible default.
      
      # However, if we could parse the samplesheet structure... but that's not in the schema.
      # The schema just says "input: path to csv".
      
      # Future improvement: parse assets/schema_input.json or similar if available?
    end
    
    required_cols
  rescue => e
    warn "Failed to fetch schema for #{pipeline_name}: #{e.message}"
    ['Name', 'Read1']
  end
  
  def self.fetch_with_cache(url, cache_path)
    FileUtils.mkdir_p(File.dirname(cache_path))
    
    if File.exist?(cache_path) && (Time.now - File.mtime(cache_path)) < CACHE_TTL
      return JSON.parse(File.read(cache_path))
    end
    
    begin
      puts "Fetching #{url}..."
      data = URI.open(url).read
      File.write(cache_path, data)
      JSON.parse(data)
    rescue OpenURI::HTTPError => e
      warn "HTTP Error fetching #{url}: #{e.message}"
      if File.exist?(cache_path)
        warn "Using stale cache for #{url}"
        JSON.parse(File.read(cache_path))
      else
        nil
      end
    rescue => e
      warn "Error fetching #{url}: #{e.message}"
      nil
    end
  end
end
