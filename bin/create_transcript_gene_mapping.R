#!/usr/bin/env Rscript
# Create a pre-built transcript-to-gene mapping file from transcriptome FASTA
# This avoids slow FASTA parsing during every pipeline run

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

opt_list <- list(
  make_option("--transcriptome_fasta", type="character", help="Transcriptome FASTA file"),
  make_option("--output", type="character", help="Output TSV file with transcript-to-gene mapping")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$transcriptome_fasta)) {
  stop("--transcriptome_fasta is required")
}

if (is.null(opt$output)) {
  stop("--output is required")
}

cat("Creating transcript-to-gene mapping from", opt$transcriptome_fasta, "\n")

# Read FASTA sequences
fasta_seqs <- readDNAStringSet(opt$transcriptome_fasta)

transcript_ids <- character(length(fasta_seqs))
gene_names <- character(length(fasta_seqs))

# Parse each header to extract gene name
for (i in seq_along(fasta_seqs)) {
  full_header <- names(fasta_seqs)[i]
  
  # Extract transcript ID (first part before space)
  transcript_id <- strsplit(full_header, " ")[[1]][1]
  transcript_ids[i] <- transcript_id
  
  # Extract gene name from header (assuming format like "ENST00000456328 gene=DDX11L2")
  gene_name <- "Unknown"
  if (grepl("gene=", full_header)) {
    # Use robust regex to extract gene name
    gene_match <- regmatches(full_header, regexpr("gene=([A-Za-z0-9_.-]+)", full_header))
    if (length(gene_match) > 0) {
      gene_name <- sub("gene=", "", gene_match)
    }
  }
  gene_names[i] <- gene_name
  
  # Progress indicator for large files
  if (i %% 10000 == 0) {
    cat(sprintf("Processed %d/%d sequences\n", i, length(fasta_seqs)))
  }
}

# Create mapping data frame
mapping <- data.frame(
  transcript_id = transcript_ids,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

# Write mapping file
write.table(mapping, opt$output, sep="\t", quote=FALSE, row.names=FALSE)

cat("Created mapping file with", nrow(mapping), "entries\n")
cat("Written to:", opt$output, "\n")

# Summary statistics
cat("\nSummary:\n")
cat("  Total transcripts:", nrow(mapping), "\n")
cat("  Transcripts with gene names:", sum(mapping$gene_name != "Unknown"), "\n")
cat("  Transcripts without gene names:", sum(mapping$gene_name == "Unknown"), "\n")
cat("  Unique genes:", length(unique(mapping$gene_name[mapping$gene_name != "Unknown"])), "\n")
