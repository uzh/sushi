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
      
      # Get resource limits from param (set by SUSHI)
      maxCpus <- ifelse(is.null(param[["cores"]]), 4, as.integer(param[["cores"]]))
      maxMemory <- ifelse(is.null(param[["ram"]]), 8, as.integer(param[["ram"]]))
      
      # Create custom config to override resource limits
      # executor memory should be >= process memory to avoid "exceeds available" errors
      executorMemory <- max(maxMemory, 64)  # At least 64 GB for nf-core pipelines
      
      customConfig <- "custom_resources.config"
      configContent <- paste0(
        "// Override executor to allow sufficient resources\n",
        "executor {\n",
        "  name = 'local'\n",
        "  cpus = ", max(maxCpus, 8), "\n",
        "  memory = '", executorMemory, ".GB'\n",
        "}\n",
        "\n",
        "// Override all process resource requirements\n",
        "process {\n",
        "  cpus = ", maxCpus, "\n",
        "  memory = '", maxMemory, ".GB'\n",
        "}\n"
      )
      writeLines(configContent, customConfig)
      message("Custom config written: ", customConfig)
      message(configContent)
      
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
      
      write.csv(samples, file=filePath, row.names=FALSE, quote=FALSE)
    }
  )
)
