#!/usr/bin/env Rscript

# Extract cDNA-complementary primers from Primer3 output
# This script selects primers that will match the mature mRNA sequence
# (and thus be complementary to cDNA synthesized from that mRNA)
# For negative-strand genes: extract LEFT primers (match antisense mRNA directly)
# For positive-strand genes: extract RIGHT primers (match sense mRNA directly)

suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option("--primer3_output", type="character", help="Primer3 output file"),
  make_option("--peaks_tsv", type="character", help="peaks.tsv file with gene strand information"),
  make_option("--out_tsv", type="character", default="cdna_primers.tsv", help="Output TSV file for cDNA-complementary primers")
)

opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$primer3_output) || is.null(opt$peaks_tsv)) {
  stop("Both --primer3_output and --peaks_tsv are required")
}

# Read peaks data to get strand information
peaks <- read.delim(opt$peaks_tsv, stringsAsFactors = FALSE)
gene_strands <- setNames(peaks$strand, peaks$gene)

# Parse Primer3 output
parse_primer3_output <- function(file_path) {
  lines <- readLines(file_path)
  results <- list()
  current_gene <- NULL
  current_data <- list()
  
  for (line in lines) {
    if (line == "=") {
      # End of current gene block
      if (!is.null(current_gene) && length(current_data) > 0) {
        results[[current_gene]] <- current_data
      }
      current_gene <- NULL
      current_data <- list()
    } else if (startsWith(line, "SEQUENCE_ID=")) {
      current_gene <- sub("SEQUENCE_ID=", "", line)
      current_data <- list()
    } else if (grepl("^PRIMER_(LEFT|RIGHT)_[0-9]+_SEQUENCE=", line)) {
      # Parse primer sequence lines
      parts <- strsplit(line, "=")[[1]]
      key <- parts[1]
      sequence <- parts[2]
      
      # Extract primer type (LEFT/RIGHT) and index
      if (grepl("PRIMER_LEFT_", key)) {
        primer_type <- "LEFT"
        primer_idx <- as.numeric(sub("PRIMER_LEFT_([0-9]+)_SEQUENCE", "\\1", key))
      } else if (grepl("PRIMER_RIGHT_", key)) {
        primer_type <- "RIGHT"
        primer_idx <- as.numeric(sub("PRIMER_RIGHT_([0-9]+)_SEQUENCE", "\\1", key))
      }
      
      current_data[[paste0(primer_type, "_", primer_idx)]] <- list(
        type = primer_type,
        index = primer_idx,
        sequence = sequence
      )
    }
  }
  
  # Don't forget the last gene if file doesn't end with =
  if (!is.null(current_gene) && length(current_data) > 0) {
    results[[current_gene]] <- current_data
  }
  
  return(results)
}

# Parse primer data
primer_data <- parse_primer3_output(opt$primer3_output)

# Extract cDNA-complementary primers
cdna_primers <- data.frame(
  gene_id = character(),
  primer_index = integer(),
  primer_type = character(),
  primer_sequence = character(),
  gene_strand = character(),
  rationale = character(),
  stringsAsFactors = FALSE
)

for (gene_id in names(primer_data)) {
  if (!gene_id %in% names(gene_strands)) {
    cat("Warning: Gene", gene_id, "not found in peaks data, skipping\n")
    next
  }
  
  strand <- gene_strands[[gene_id]]
  primers <- primer_data[[gene_id]]
  
  # Determine which primer type to extract based on strand
  # Since transcriptome FASTA contains mature mRNA sequences in 5' to 3' orientation,
  # and our genomic template extraction is strand-aware (getSeq returns mRNA-matching sequences),
  # we always want primers that match the template directly = LEFT primers
  # This ensures primers match the mRNA sequence (to be complementary to cDNA)
  
  # Both positive and negative strand genes should use LEFT primers
  # because the extracted genomic template always matches the mRNA orientation
  target_type <- "LEFT"
  rationale <- paste0("LEFT primer matches mRNA sequence (", 
                     ifelse(strand == "+", "positive", "negative"), 
                     "-strand gene: genomic template matches mRNA)")
  
  # Note: This logic assumes strand-aware sequence extraction by getSeq()
  
  # Extract primers of target type
  for (primer_name in names(primers)) {
    primer_info <- primers[[primer_name]]
    if (primer_info$type == target_type) {
      cdna_primers <- rbind(cdna_primers, data.frame(
        gene_id = gene_id,
        primer_index = primer_info$index,
        primer_type = primer_info$type,
        primer_sequence = primer_info$sequence,
        gene_strand = strand,
        rationale = rationale,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sort by gene_id and primer_index
cdna_primers <- cdna_primers[order(cdna_primers$gene_id, cdna_primers$primer_index), ]

# Write output
write.table(cdna_primers, opt$out_tsv, sep="\t", quote=FALSE, row.names=FALSE)

cat("Extracted", nrow(cdna_primers), "cDNA-complementary primers for", length(unique(cdna_primers$gene_id)), "genes\n")
cat("Output written to:", opt$out_tsv, "\n")

# Summary by strand
neg_count <- sum(cdna_primers$gene_strand == "-")
pos_count <- sum(cdna_primers$gene_strand == "+")
cat("Summary:\n")
cat("  Negative-strand genes (LEFT primers):", neg_count, "\n")
cat("  Positive-strand genes (RIGHT primers):", pos_count, "\n")
