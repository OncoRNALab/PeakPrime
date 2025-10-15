#!/usr/bin/env Rscript

# Fetch MANE Select transcript sequences using biomaRt
# This script reads a list of Ensembl gene IDs and fetches their MANE Select
# transcript sequences from Ensembl using the biomaRt package.

suppressPackageStartupMessages({
  library(biomaRt)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--gene-ids"), type = "character", default = NULL,
              help = "Input file with one Ensembl gene ID per line", metavar = "FILE"),
  make_option(c("--output-fasta"), type = "character", default = "mane_transcripts.fasta",
              help = "Output FASTA file with transcript sequences", metavar = "FILE"),
  make_option(c("--output-mapping"), type = "character", default = "mane_mapping.tsv",
              help = "Output TSV file mapping gene IDs to transcript IDs", metavar = "FILE"),
  make_option(c("--ensembl-version"), type = "integer", default = NULL,
              help = "Ensembl version (default: current version)", metavar = "INT"),
  make_option(c("--allow-fallback"), action = "store_true", default = TRUE,
              help = "Allow fallback to canonical transcript if no MANE Select [default]"),
  make_option(c("--no-fallback"), action = "store_false", dest = "allow_fallback",
              help = "Require MANE Select transcripts only (no fallback)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure allow_fallback is properly set
if (is.null(opt$allow_fallback)) {
  opt$allow_fallback <- TRUE
}

# Validate required arguments
if (is.null(opt$`gene-ids`)) {
  print_help(opt_parser)
  stop("--gene-ids is required", call. = FALSE)
}

# Read gene IDs
cat("Reading gene IDs from:", opt$`gene-ids`, "\n", file = stderr())
gene_ids <- readLines(opt$`gene-ids`)
gene_ids <- gene_ids[nchar(gene_ids) > 0]  # Remove empty lines
cat("Found", length(gene_ids), "gene IDs\n", file = stderr())

# Connect to Ensembl
cat("Connecting to Ensembl biomaRt...\n", file = stderr())
if (is.null(opt$`ensembl-version`)) {
  ensembl <- useEnsembl(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl")
  cat("Using current Ensembl version\n", file = stderr())
} else {
  ensembl <- useEnsembl(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl", 
                        version = opt$`ensembl-version`)
  cat("Using Ensembl version:", opt$`ensembl-version`, "\n", file = stderr())
}

# Define attributes to retrieve
attrs <- c("ensembl_gene_id",
           "hgnc_symbol",
           "ensembl_transcript_id",
           "transcript_mane_select",
           "transcript_is_canonical",
           "cdna")

cat("Fetching transcript data from Ensembl...\n", file = stderr())
cat("This may take a while for large gene lists...\n", file = stderr())

# Fetch all transcript data for the gene list
results <- getBM(
  attributes = attrs,
  filters    = "ensembl_gene_id",
  values     = gene_ids,
  mart       = ensembl
)

cat("Retrieved", nrow(results), "transcript records\n", file = stderr())

# Process results to select MANE or fallback transcripts
selected_transcripts <- data.frame()
mane_count <- 0
canonical_count <- 0
failed_genes <- character()

for (gene_id in gene_ids) {
  gene_data <- results[results$ensembl_gene_id == gene_id, ]
  
  if (nrow(gene_data) == 0) {
    cat("Warning: No data found for", gene_id, "\n", file = stderr())
    failed_genes <- c(failed_genes, gene_id)
    next
  }
  
  # Strategy 1: Try MANE Select
  mane_transcripts <- gene_data[gene_data$transcript_mane_select != "", ]
  
  if (nrow(mane_transcripts) > 0) {
    # Use MANE Select transcript
    selected <- mane_transcripts[1, ]
    selected$transcript_type <- "MANE_Select"
    selected$transcript_id_used <- selected$transcript_mane_select
    selected_transcripts <- rbind(selected_transcripts, selected)
    mane_count <- mane_count + 1
    cat("✓", gene_id, "->", selected$transcript_mane_select, "(MANE Select)\n", file = stderr())
    
  } else if (opt$allow_fallback) {
    # Strategy 2: Try canonical transcript
    canonical_transcripts <- gene_data[gene_data$transcript_is_canonical == 1, ]
    
    if (nrow(canonical_transcripts) > 0) {
      selected <- canonical_transcripts[1, ]
      selected$transcript_type <- "canonical"
      selected$transcript_id_used <- selected$ensembl_transcript_id
      selected_transcripts <- rbind(selected_transcripts, selected)
      canonical_count <- canonical_count + 1
      cat("Info:", gene_id, "->", selected$ensembl_transcript_id, 
          "(canonical, no MANE)\n", file = stderr())
      
    } else {
      # Strategy 3: Use first available transcript
      selected <- gene_data[1, ]
      selected$transcript_type <- "first_available"
      selected$transcript_id_used <- selected$ensembl_transcript_id
      selected_transcripts <- rbind(selected_transcripts, selected)
      cat("Info:", gene_id, "->", selected$ensembl_transcript_id, 
          "(first available, no MANE/canonical)\n", file = stderr())
    }
    
  } else {
    # No MANE and fallback not allowed
    cat("Warning: No MANE Select transcript for", gene_id, "(fallback disabled)\n", file = stderr())
    failed_genes <- c(failed_genes, gene_id)
  }
}

# Check if we have any successful results
if (nrow(selected_transcripts) == 0) {
  cat("\nERROR: No transcripts were fetched successfully\n", file = stderr())
  cat("Failed genes:", paste(failed_genes, collapse = ", "), "\n", file = stderr())
  quit(status = 1)
}

# Write FASTA output
cat("\nWriting FASTA file:", opt$`output-fasta`, "\n", file = stderr())
fasta_conn <- file(opt$`output-fasta`, "w")

for (i in 1:nrow(selected_transcripts)) {
  row <- selected_transcripts[i, ]
  
  # Skip if no sequence available
  if (is.na(row$cdna) || row$cdna == "") {
    cat("Warning: No cDNA sequence for", row$transcript_id_used, "\n", file = stderr())
    next
  }
  
  # Write FASTA header with | separator for primer3 compatibility
  # Format: >GENE_ID|transcript:TRANSCRIPT_ID gene_name:GENE_NAME type:TYPE
  header <- paste0(">", row$ensembl_gene_id, "|", 
                   "transcript:", row$transcript_id_used, " ",
                   "gene_name:", row$hgnc_symbol, " ",
                   "type:", row$transcript_type)
  writeLines(header, fasta_conn)
  
  # Write sequence
  writeLines(row$cdna, fasta_conn)
}

close(fasta_conn)

# Write mapping file
cat("Writing mapping file:", opt$`output-mapping`, "\n", file = stderr())
mapping <- data.frame(
  gene_id = selected_transcripts$ensembl_gene_id,
  transcript_id = selected_transcripts$transcript_id_used,
  gene_name = selected_transcripts$hgnc_symbol,
  transcript_type = selected_transcripts$transcript_type,
  stringsAsFactors = FALSE
)

write.table(mapping, 
            file = opt$`output-mapping`,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Print summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n", file = stderr())
cat("Summary:\n", file = stderr())
cat("  Total genes queried:", length(gene_ids), "\n", file = stderr())
cat("  Successful:", nrow(selected_transcripts), "\n", file = stderr())
cat("  Failed:", length(failed_genes), "\n", file = stderr())
cat("  MANE Select:", mane_count, "\n", file = stderr())
cat("  Canonical:", canonical_count, "\n", file = stderr())
cat("  Other:", nrow(selected_transcripts) - mane_count - canonical_count, "\n", file = stderr())
cat(paste(rep("=", 60), collapse = ""), "\n", file = stderr())

if (length(failed_genes) > 0) {
  cat("\nFailed genes:", paste(failed_genes, collapse = ", "), "\n", file = stderr())
}

cat("\nSUCCESS: Fetched", nrow(selected_transcripts), "transcripts\n", file = stderr())
