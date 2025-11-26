#!/usr/bin/env Rscript

# Extract cDNA-complementary primers from Primer3 output
# This script selects primers that will match the mature mRNA sequence
# (and thus be complementary to cDNA synthesized from that mRNA)
# 
# ALWAYS extracts LEFT primers for ALL genes (both + and - strand)
# This is correct because the genomic template sequence is already strand-corrected
# in the previous step (process_macs2_peaks.R uses getSeq with strand-aware extraction)
# so the template always matches the mRNA orientation (5' to 3')

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
peaks <- read.delim(opt$peaks_tsv, sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Validate peaks table
if (is.null(peaks) || nrow(peaks) == 0) {
  cat("No peaks found in", opt$peaks_tsv, "- writing empty output and exiting\n")
  write.table(data.frame(), opt$out_tsv, sep="\t", quote=FALSE, row.names=FALSE)
  quit(status=0)
}
cat("Peaks table columns:", paste(colnames(peaks), collapse=", "), "\n")
required_cols <- c("gene", "strand")
missing_cols <- setdiff(required_cols, colnames(peaks))
if (length(missing_cols) > 0) {
  stop("Peaks file is missing required columns: ", paste(missing_cols, collapse=", "))
}

# Create lookup for strand information
# In multi-peak mode, we need to handle multiple peaks per gene
# Key format: "ENSG00000123456|peak_N" or just "ENSG00000123456" for single-peak
if ("peak_rank" %in% colnames(peaks)) {
  # Multi-peak mode: create keys with peak identifiers
  peak_keys <- paste0(peaks$gene, "|peak_", peaks$peak_rank)
  gene_strands <- setNames(as.character(peaks$strand), peak_keys)
  # Also keep gene-only keys for backward compatibility (use first occurrence)
  gene_only <- setNames(as.character(peaks$strand), as.character(peaks$gene))
  gene_strands <- c(gene_strands, gene_only[!names(gene_only) %in% names(gene_strands)])
} else {
  # Single-peak mode: original behavior
  gene_strands <- setNames(as.character(peaks$strand), as.character(peaks$gene))
}

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

# If Primer3 produced no records, write empty output and exit gracefully
if (is.null(primer_data) || length(primer_data) == 0) {
  cat("No primer records parsed from", opt$primer3_output, "- writing empty output\n")
  write.table(data.frame(), opt$out_tsv, sep="\t", quote=FALSE, row.names=FALSE)
  quit(status=0)
}

# Extract cDNA-complementary primers
cdna_primers <- data.frame(
  gene_id = character(),
  peak_id = character(),
  primer_index = integer(),
  primer_type = character(),
  primer_sequence = character(),
  gene_strand = character(),
  rationale = character(),
  stringsAsFactors = FALSE
)

for (sequence_id in names(primer_data)) {
  # Parse sequence_id to extract gene_id and peak_id (if present)
  # Format can be:
  #   - "ENSG00000123456|peak_1" (multi-peak mode)
  #   - "ENSG00000123456" (single-peak mode)
  if (grepl("\\|peak_", sequence_id)) {
    # Multi-peak format
    parts <- strsplit(sequence_id, "\\|")[[1]]
    gene_id <- parts[1]
    peak_id <- parts[2]
  } else {
    # Single-peak format
    gene_id <- sequence_id
    peak_id <- NA
  }
  
  # Look up strand information
  # Try sequence_id first (includes peak_id), fall back to gene_id
  if (sequence_id %in% names(gene_strands)) {
    strand <- gene_strands[[sequence_id]]
  } else if (gene_id %in% names(gene_strands)) {
    strand <- gene_strands[[gene_id]]
  } else {
    cat("Warning: Neither sequence_id '", sequence_id, "' nor gene_id '", gene_id, 
        "' found in peaks data, skipping\n", sep="")
    next
  }
  
  primers <- primer_data[[sequence_id]]
  
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
        peak_id = peak_id,
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

# Sort by gene_id, peak_id, and primer_index
cdna_primers <- cdna_primers[order(cdna_primers$gene_id, cdna_primers$peak_id, cdna_primers$primer_index), ]

# Write output
write.table(cdna_primers, opt$out_tsv, sep="\t", quote=FALSE, row.names=FALSE)

cat("Extracted", nrow(cdna_primers), "cDNA-complementary primers for", length(unique(cdna_primers$gene_id)), "genes\n")

# Count peaks if peak_id column has data
n_peaks <- sum(!is.na(cdna_primers$peak_id))
if (n_peaks > 0) {
  cat("  Multi-peak mode:", length(unique(cdna_primers$peak_id[!is.na(cdna_primers$peak_id)])), "unique peaks\n")
}

cat("Output written to:", opt$out_tsv, "\n")

# Summary by strand
neg_count <- sum(cdna_primers$gene_strand == "-")
pos_count <- sum(cdna_primers$gene_strand == "+")
cat("Summary:\n")
cat("  Negative-strand genes (LEFT primers):", neg_count, "\n")
cat("  Positive-strand genes (LEFT primers):", pos_count, "\n")
