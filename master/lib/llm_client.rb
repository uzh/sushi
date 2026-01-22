# LLM Client for Materials & Methods Generation
# Uses OpenAI-compatible API (local vLLM server)

require 'net/http'
require 'json'
require 'uri'

class LlmClient
  DEFAULT_ENDPOINT = "http://fgcz-c-056:8000/v1/chat/completions"
  DEFAULT_MODEL = "LLM"
  DEFAULT_MAX_TOKENS = 4096
  DEFAULT_TEMPERATURE = 0.3

  attr_accessor :endpoint, :model, :max_tokens, :temperature, :timeout

  def initialize(options = {})
    @endpoint = options[:endpoint] || DEFAULT_ENDPOINT
    @model = options[:model] || DEFAULT_MODEL
    @max_tokens = options[:max_tokens] || DEFAULT_MAX_TOKENS
    @temperature = options[:temperature] || DEFAULT_TEMPERATURE
    @timeout = options[:timeout] || 120  # seconds
  end

  # Generate Materials & Methods from analysis data
  def generate_mm(analysis_data)
    system_prompt = build_system_prompt
    user_prompt = build_user_prompt(analysis_data)
    
    response = chat_completion(system_prompt, user_prompt)
    
    if response[:success]
      response[:content]
    else
      # Fallback to template-based generation
      Rails.logger.error("LLM API failed: #{response[:error]}")
      generate_fallback_mm(analysis_data)
    end
  end

  # Raw chat completion API call
  def chat_completion(system_prompt, user_prompt)
    uri = URI.parse(@endpoint)
    
    request_body = {
      model: @model,
      messages: [
        { role: "system", content: system_prompt },
        { role: "user", content: user_prompt }
      ],
      max_tokens: @max_tokens,
      temperature: @temperature
    }

    begin
      http = Net::HTTP.new(uri.host, uri.port)
      http.read_timeout = @timeout
      http.open_timeout = 10
      
      request = Net::HTTP::Post.new(uri.path)
      request["Content-Type"] = "application/json"
      request.body = request_body.to_json
      
      response = http.request(request)
      
      if response.code.to_i == 200
        result = JSON.parse(response.body)
        content = result.dig("choices", 0, "message", "content")
        {
          success: true,
          content: content,
          usage: result["usage"]
        }
      else
        {
          success: false,
          error: "HTTP #{response.code}: #{response.body}"
        }
      end
    rescue => e
      {
        success: false,
        error: "#{e.class}: #{e.message}"
      }
    end
  end

  private

  def build_system_prompt
    <<~PROMPT
      You are a scientific writing assistant specializing in bioinformatics and genomics research.
      Your task is to write a Materials and Methods section for a scientific publication based on the provided analysis information.
      
      Guidelines:
      - Write in formal scientific prose suitable for peer-reviewed journals
      - Use past tense and passive voice where appropriate
      - Include specific software versions, parameters, and reference genome information when available
      - Organize the content into logical subsections (e.g., Data Processing, Quality Control, Alignment, Quantification)
      - Be precise and reproducible - another researcher should be able to replicate the analysis
      - Do not include placeholder text or instructions
      - Output in Markdown format with appropriate headers (##, ###)
      - Keep the text concise but comprehensive (typically 300-600 words)
      - If certain information is not provided, do not fabricate it - only include what is given
      
      Focus on:
      1. Sample information and experimental design
      2. Data processing pipeline and tools used
      3. Reference genome and annotation details
      4. Key analysis parameters
      5. Quality control measures
    PROMPT
  end

  def build_user_prompt(data)
    sections = []
    
    sections << "Please write a Materials and Methods section based on the following analysis information:\n"
    
    # Basic analysis info
    sections << "## Analysis Overview"
    sections << "- Analysis Name: #{data[:name]}" if data[:name]
    sections << "- Analysis Type/Application: #{data[:sushi_app]}" if data[:sushi_app]
    sections << "- Analysis Date: #{data[:date]}" if data[:date]
    sections << "- Number of Samples: #{data[:sample_count]}" if data[:sample_count]
    sections << "- Project: #{data[:project]}" if data[:project]
    sections << ""
    
    # Parameters
    if data[:parameters] && !data[:parameters].empty?
      sections << "## Analysis Parameters"
      data[:parameters].each do |key, value|
        next if value.to_s.empty?
        sections << "- #{key}: #{value}"
      end
      sections << ""
    end
    
    # Reference information
    if data[:ref_build] || data[:ref_feature]
      sections << "## Reference Information"
      sections << "- Reference Build: #{data[:ref_build]}" if data[:ref_build]
      sections << "- Annotation/Feature File: #{data[:ref_feature]}" if data[:ref_feature]
      sections << ""
    end
    
    # Job script content (shows actual commands)
    if data[:job_script] && !data[:job_script].empty?
      sections << "## Job Script (actual commands executed)"
      sections << "```bash"
      # Truncate if too long
      script = data[:job_script]
      if script.length > 8000
        script = script[0, 8000] + "\n... (truncated)"
      end
      sections << script
      sections << "```"
      sections << ""
    end
    
    # Input dataset info
    if data[:input_samples] && !data[:input_samples].empty?
      sections << "## Input Samples"
      sections << "Sample names from input dataset:"
      data[:input_samples].each do |sample|
        sections << "- #{sample}"
      end
      sections << ""
    end
    
    # Output info
    if data[:output_files] && !data[:output_files].empty?
      sections << "## Output Files Generated"
      data[:output_files].each do |file|
        sections << "- #{file}"
      end
      sections << ""
    end
    
    sections.join("\n")
  end

  def generate_fallback_mm(data)
    # Simple template-based fallback when LLM is unavailable
    mm = []
    mm << "# Materials and Methods"
    mm << ""
    mm << "## Data Analysis"
    mm << ""
    
    if data[:sushi_app]
      mm << "Data analysis was performed using #{data[:sushi_app]} via the SUSHI bioinformatics platform (Functional Genomics Center Zurich)."
      mm << ""
    end
    
    if data[:ref_build]
      mm << "Sequences were aligned to the #{data[:ref_build]} reference genome."
      mm << ""
    end
    
    if data[:parameters] && !data[:parameters].empty?
      mm << "### Analysis Parameters"
      mm << ""
      data[:parameters].each do |key, value|
        next if value.to_s.empty?
        next if ['cores', 'ram', 'scratch', 'node', 'process_mode', 'mail'].include?(key)
        mm << "- #{key}: #{value}"
      end
      mm << ""
    end
    
    mm << "---"
    mm << ""
    mm << "*Note: LLM generation was unavailable. This is a template-based output.*"
    mm << "*Generated by SUSHI on #{Time.now.strftime('%Y-%m-%d %H:%M:%S')}*"
    
    mm.join("\n")
  end

  # Class method for quick access
  def self.generate_mm(analysis_data, options = {})
    client = new(options)
    client.generate_mm(analysis_data)
  end
end
