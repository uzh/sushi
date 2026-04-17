# NfCoreConventions - Convention over Configuration for nf-core SUSHI Apps
#
# This module defines the conventions used to automatically generate SUSHI Apps
# from nf-core pipelines without requiring manual YAML configuration.
#
# Based on the Rails philosophy: "Convention over Configuration"
# - Most nf-core pipelines follow similar patterns
# - We define sensible defaults based on these patterns
# - YAML configuration is only needed for exceptions
#
module NfCoreConventions
  # ============================================================
  # Samplesheet Column Mapping (nf-core -> SUSHI)
  # ============================================================
  # Maps nf-core samplesheet column names to SUSHI dataset column names.
  # For columns with special handling (like strandedness), use a Hash
  # with :sushi_col and optional :default value.
  #
  COLUMN_MAPPING = {
    # Standard columns (used by most pipelines)
    'sample' => 'Name',
    'fastq_1' => 'Read1 [File]',
    'fastq_2' => 'Read2 [File]',
    'fasta' => 'FASTA [File]',
    
    # RNA-seq specific
    'strandedness' => { sushi_col: 'StrandMode', default: 'auto' },
    
    # ChIP-seq / ATAC-seq specific
    'antibody' => 'Antibody',
    'control' => 'Control',
    'control_type' => 'ControlType',
    
    # Variant calling specific (sarek, etc.)
    'patient' => 'Patient',
    'sex' => 'Sex',
    'status' => 'Status',
    'lane' => 'Lane',
    
    # General metadata
    'condition' => 'Condition',
    'replicate' => 'Replicate',
    'single_end' => 'PairedEnd',  # Note: inverse logic
    
    # Long-read sequencing
    'bam' => 'BAM [File]',
    'bai' => 'BAI [File]',
    
    # Proteomics
    'experiment' => 'Experiment',
    'fraction' => 'Fraction'
  }.freeze

  # ============================================================
  # Reference Selector Trigger Parameters
  # ============================================================
  # If any of these parameters exist in the "reference_genome_options"
  # section of nextflow_schema.json, the pipeline needs a reference selector.
  #
  REF_SELECTOR_PARAMS = %w[fasta gtf genome gff transcript_fasta].freeze
  
  # Section name patterns that indicate reference genome options
  REF_SELECTOR_SECTION_PATTERNS = [
    /reference.*genome/i,
    /genome.*options/i,
    /reference.*options/i
  ].freeze

  # ============================================================
  # Category-based Resource Defaults
  # ============================================================
  # Different pipeline categories have different typical resource requirements.
  # These defaults can be overridden in nf_core_pipelines.yml for specific pipelines.
  #
  RESOURCE_DEFAULTS = {
    'Transcriptomics' => { 'cores' => 8, 'ram' => 64, 'scratch' => 200 },
    'Genomics' => { 'cores' => 8, 'ram' => 64, 'scratch' => 200 },
    'SingleCell' => { 'cores' => 16, 'ram' => 128, 'scratch' => 500 },
    'Epigenetics' => { 'cores' => 8, 'ram' => 30, 'scratch' => 100 },
    'Metagenomics' => { 'cores' => 8, 'ram' => 64, 'scratch' => 200 },
    'Proteomics' => { 'cores' => 8, 'ram' => 64, 'scratch' => 200 },
    'QC' => { 'cores' => 4, 'ram' => 16, 'scratch' => 100 },
    'DataFetch' => { 'cores' => 4, 'ram' => 16, 'scratch' => 500 },
    'default' => { 'cores' => 8, 'ram' => 30, 'scratch' => 100 }
  }.freeze

  # ============================================================
  # Category Mapping (nf-core topics -> SUSHI categories)
  # ============================================================
  # Maps nf-core pipeline topics (from pipelines.json) to SUSHI analysis categories.
  #
  CATEGORY_MAPPING = {
    # RNA-seq related
    'rna-seq' => 'Transcriptomics',
    'rna' => 'Transcriptomics',
    'transcriptomics' => 'Transcriptomics',
    'gene-expression' => 'Transcriptomics',
    
    # DNA-seq related
    'dna-seq' => 'Genomics',
    'dna' => 'Genomics',
    'genomics' => 'Genomics',
    'variant-calling' => 'Genomics',
    'whole-genome-sequencing' => 'Genomics',
    'exome-sequencing' => 'Genomics',
    
    # Epigenetics
    'chip-seq' => 'Epigenetics',
    'atac-seq' => 'Epigenetics',
    'methylation' => 'Epigenetics',
    'epigenetics' => 'Epigenetics',
    'histone-modification' => 'Epigenetics',
    
    # Single cell
    'sc-rna-seq' => 'SingleCell',
    'single-cell' => 'SingleCell',
    'spatial-transcriptomics' => 'SingleCell',
    
    # Metagenomics
    'metagenomics' => 'Metagenomics',
    'microbiome' => 'Metagenomics',
    'amplicon-sequencing' => 'Metagenomics',
    
    # Proteomics
    'proteomics' => 'Proteomics',
    'mass-spectrometry' => 'Proteomics',
    
    # QC
    'quality-control' => 'QC',
    'qc' => 'QC',
    
    # Data fetching
    'data-management' => 'DataFetch',
    'download' => 'DataFetch'
  }.freeze

  # ============================================================
  # Helper Methods
  # ============================================================
  
  # Get SUSHI column name for nf-core column
  # @param nf_col [String] nf-core column name
  # @return [String, nil] SUSHI column name or nil if not mapped
  def self.sushi_column_for(nf_col)
    mapping = COLUMN_MAPPING[nf_col]
    return nil unless mapping
    mapping.is_a?(Hash) ? mapping[:sushi_col] : mapping
  end
  
  # Get default value for nf-core column (if any)
  # @param nf_col [String] nf-core column name
  # @return [String, nil] default value or nil
  def self.default_for(nf_col)
    mapping = COLUMN_MAPPING[nf_col]
    return nil unless mapping.is_a?(Hash)
    mapping[:default]
  end
  
  # Get resource defaults for a category
  # @param category [String] SUSHI analysis category
  # @return [Hash] resource defaults (cores, ram, scratch)
  def self.resources_for(category)
    RESOURCE_DEFAULTS[category] || RESOURCE_DEFAULTS['default']
  end
  
  # Map nf-core topic to SUSHI category
  # @param topic [String] nf-core topic
  # @return [String, nil] SUSHI category or nil if not mapped
  def self.category_for(topic)
    CATEGORY_MAPPING[topic.to_s.downcase]
  end
  
  # Check if section name indicates reference genome options
  # @param section_name [String] section name from nextflow_schema.json
  # @return [Boolean]
  def self.is_ref_section?(section_name)
    REF_SELECTOR_SECTION_PATTERNS.any? { |pattern| section_name =~ pattern }
  end
end
