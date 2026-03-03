EzAppNfCoreGeneric <- setRefClass(
  "EzAppNfCoreGeneric",
  contains = "EzApp",
  methods = list(
    initialize = function() {
      "Initializes the application"
      callSuper()
    },
    run = function(input, output, param) {
      pipeline <- param[["nfcorePipeline"]]
      version <- param[["pipelineVersion"]]
      dataRoot <- param[["dataRoot"]]
      
      # Output to current directory (scratch)
      # Results will be copied to gstore by the job script footer
      outdir <- "_result"
      
      # Ensure output directory exists
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      # Build sample sheet in current directory (scratch)
      sampleSheetPath <- "samplesheet.csv"
      writeNfCoreSampleSheet(input, sampleSheetPath, param)
      
      message("Sample sheet written to: ", sampleSheetPath)
      message("Sample sheet contents:")
      print(read.csv(sampleSheetPath))
      
      # Get resource limits from param (set by SUSHI, auto-detected from pipeline schema)
      maxCpus <- ifelse(is.null(param[["cores"]]), 16, as.integer(param[["cores"]]))
      maxMemory <- ifelse(is.null(param[["ram"]]), 128, as.integer(param[["ram"]]))
      
      customConfig <- "custom_resources.config"
      configContent <- paste0(
        "// Executor: total available resources for this job\n",
        "executor {\n",
        "  name = 'local'\n",
        "  cpus = ", maxCpus, "\n",
        "  memory = '", maxMemory, ".GB'\n",
        "}\n",
        "\n",
        "// Cap all process labels to fit within executor limits\n",
        "process {\n",
        "  resourceLimits = [\n",
        "    cpus: ", maxCpus, ",\n",
        "    memory: '", maxMemory, ".GB'\n",
        "  ]\n",
        "}\n"
      )
      writeLines(configContent, customConfig)
      message("Custom config written: ", customConfig)
      message(configContent)
      
      # Load Nextflow module via lmod (FGCZ environment)
      # Source the lmod R init script to get the module() function
      source("/usr/share/lmod/lmod/init/R")

      # Get module version from param (default to 25.10)
      nfModuleVersion <- param[["nextflowModuleVersion"]]
      if (is.null(nfModuleVersion) || nfModuleVersion == "") {
        nfModuleVersion <- "25.10"
      }

      # Load the nextflow module
      # Note: lmod's module() uses match.call() which captures variable names literally,
      # so we must use do.call() to pass evaluated values
      moduleName <- paste0("nextflow/", nfModuleVersion)
      message("Loading module: ", moduleName)
      do.call(module, list("load", moduleName))

      # Construct nextflow command
      # Output to local directory, will be copied to gstore by footer
      cmd <- paste(
        "nextflow run",
        paste0("nf-core/", pipeline),
        "-r", version,
        "--input", sampleSheetPath,
        "--outdir", outdir,
        "-profile singularity",
        "-c", customConfig
      )
      
      # Resolve reference genome from refBuild parameter
      refBuild <- param[["refBuild"]]
      if (!is.null(refBuild) && refBuild != "" && refBuild != "select") {
        genomesRoots <- strsplit(GENOMES_ROOT, ":")[[1]]
        refParts <- strsplit(refBuild, "/")[[1]]
        genomeBase <- paste(refParts[1:min(3, length(refParts))], collapse = "/")
        
        refRoot <- NULL
        for (root in genomesRoots) {
          if (dir.exists(file.path(root, genomeBase))) {
            refRoot <- root
            break
          }
        }
        
        if (!is.null(refRoot)) {
          fastaPath <- file.path(refRoot, genomeBase, "Sequence", "WholeGenomeFasta", "genome.fa")
          if (file.exists(fastaPath)) {
            cmd <- paste(cmd, "--fasta", fastaPath)
            message("Using FASTA: ", fastaPath)
          }
          
          gtfPath <- NULL
          if (length(refParts) > 3) {
            gtfPath <- file.path(refRoot, refBuild, "Genes", "genes.gtf")
          }
          if (is.null(gtfPath) || !file.exists(gtfPath)) {
            gtfPath <- file.path(refRoot, genomeBase, "Annotation", "Genes", "genes.gtf")
          }
          if (file.exists(gtfPath)) {
            cmd <- paste(cmd, "--gtf", gtfPath)
            message("Using GTF: ", gtfPath)
          }
        } else {
          message("WARNING: Reference genome not found for refBuild: ", refBuild)
        }
      }
      
      message("Running command: ", cmd)
      
      ezSystem(cmd)
    },
    writeNfCoreSampleSheet = function(input, filePath, param) {
      # Handle input
      if (is.character(input) && file.exists(input)) {
        # Dataset mode: input is path to TSV
        ds <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
      } else {
        # Sample mode: input is list
        ds <- as.data.frame(input, stringsAsFactors=FALSE)
      }
      
      dataRoot <- param[["dataRoot"]]
      
      # Helper to resolve file paths
      resolve_path <- function(path) {
        if (is.na(path) || path == "") return(NA)
        if (startsWith(path, "/")) return(path)
        return(file.path(dataRoot, path))
      }
      
      # Map SUSHI columns to nf-core samplesheet columns
      # nf-core/demo expects: sample, fastq_1, fastq_2
      
      # Find Read1 column (may have suffix like [File])
      read1Col <- grep("^Read1", names(ds), value = TRUE)[1]
      read2Col <- grep("^Read2", names(ds), value = TRUE)[1]
      
      samples <- data.frame(
        sample = ds$Name,
        fastq_1 = sapply(ds[[read1Col]], resolve_path),
        stringsAsFactors = FALSE
      )
      
      if (!is.null(read2Col) && !is.na(read2Col)) {
        read2Values <- ds[[read2Col]]
        if (!all(is.na(read2Values)) && !all(read2Values == "")) {
          samples$fastq_2 <- sapply(read2Values, resolve_path)
        }
      }
      
      # Inject strandedness from parameter if set (nf-core/rnaseq requires this column)
      strandedness <- param[["strandedness"]]
      if (!is.null(strandedness) && strandedness != "") {
        samples$strandedness <- rep(strandedness, nrow(samples))
        message("Added strandedness column: ", strandedness)
      }
      
      write.csv(samples, file=filePath, row.names=FALSE, quote=FALSE)
    }
  )
)
