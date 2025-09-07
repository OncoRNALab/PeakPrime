#!/usr/bin/env Rscript

# Filter gene list to only include genes that have peaks in selected_peaks.tsv
# This prevents plotting errors for genes that were filtered out by peak quality thresholds

library(optparse)
library(data.table)

option_list <- list(
  make_option("--peaks_tsv", type="character", help="selected_peaks.tsv file"),
  make_option("--genes_file", type="character", help="Original gene list file"),
  make_option("--output", type="character", help="Output filtered gene list file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read files
peaks <- fread(opt$peaks_tsv)
original_genes <- readLines(opt$genes_file)

# Find genes with peaks
genes_with_peaks <- unique(peaks$gene)
filtered_genes <- intersect(original_genes, genes_with_peaks)
missing_genes <- setdiff(original_genes, genes_with_peaks)

# Report filtering
cat("Original genes:", length(original_genes), "\n")
cat("Genes with peaks:", length(filtered_genes), "\n")
cat("Genes removed (no peaks):", length(missing_genes), "\n")

if (length(missing_genes) > 0) {
  cat("Missing genes:\n")
  cat(paste(missing_genes, collapse="\n"), "\n")
}

# Write filtered list
writeLines(filtered_genes, opt$output)
cat("Filtered gene list written to:", opt$output, "\n")
