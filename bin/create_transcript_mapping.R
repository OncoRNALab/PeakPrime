#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

opt_list <- list(
  make_option("--transcriptome_fasta", type="character", help="Transcriptome FASTA file"),
  make_option("--out_mapping", type="character", default="transcript_to_gene_mapping.tsv", help="Output mapping file")
)

opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$transcriptome_fasta)) {
  stop("--transcriptome_fasta is required")
}

cat("Creating transcript-to-gene mapping from", opt$transcriptome_fasta, "\n")

# Read FASTA sequences
fasta_seqs <- readDNAStringSet(opt$transcriptome_fasta)

# Parse headers
transcript_ids <- character()
gene_names <- character()

for (i in seq_along(fasta_seqs)) {
  full_header <- names(fasta_seqs)[i]
  
  # Extract transcript ID (first part before space)
  transcript_id <- strsplit(full_header, " ")[[1]][1]
  
  # Extract gene name
  gene_name <- "Unknown"
  if (grepl("gene=", full_header)) {
    gene_match <- regmatches(full_header, regexpr("gene=([A-Za-z0-9_.-]+)", full_header))
    if (length(gene_match) > 0) {
      gene_name <- sub("gene=", "", gene_match)
    }
  }
  
  transcript_ids <- c(transcript_ids, transcript_id)
  gene_names <- c(gene_names, gene_name)
  
  if (i %% 10000 == 0) {
    cat("Processed", i, "transcripts\n")
  }
}

# Create mapping dataframe
mapping_df <- data.frame(
  transcript_id = transcript_ids,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

# Write mapping file
write.table(mapping_df, opt$out_mapping, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Created mapping file with", nrow(mapping_df), "entries:", opt$out_mapping, "\n")
