#!/usr/bin/env Rscript

# Analyze comprehensive QC summary from MACS2 peak processing
# Shows exactly why each gene failed at each filtering step

library(optparse)
library(data.table)

option_list <- list(
  make_option("--qc_file", type="character", help="QC summary file (peaks_qc_summary.tsv)"),
  make_option("--output_summary", type="character", default="filtering_summary.txt", help="Output summary file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read QC data
qc <- fread(opt$qc_file)

# Summary statistics
total_genes <- nrow(qc)
final_selected <- sum(qc$final_selection)

cat("=== MACS2 Peak Selection Summary ===\n")
cat("Total target genes:", total_genes, "\n")
cat("Final selected genes:", final_selected, "\n")
cat("Success rate:", round(100 * final_selected / total_genes, 1), "%\n\n")

# Filtering cascade analysis
cat("=== Filtering Cascade ===\n")
cat("1. Genes with ANY MACS2 peaks:", sum(qc$has_macs2_peaks), "/", total_genes, 
    "(", round(100 * sum(qc$has_macs2_peaks) / total_genes, 1), "%)\n")

cat("2. Genes with SIGNIFICANT peaks:", sum(qc$has_significant_peaks), "/", total_genes,
    "(", round(100 * sum(qc$has_significant_peaks) / total_genes, 1), "%)\n")

cat("3. Genes with OVERLAPPING peaks:", sum(qc$has_overlapping_peaks), "/", total_genes,
    "(", round(100 * sum(qc$has_overlapping_peaks) / total_genes, 1), "%)\n")

cat("4. Genes with BEST peak selected:", sum(qc$is_best_peak), "/", total_genes,
    "(", round(100 * sum(qc$is_best_peak) / total_genes, 1), "%)\n")

cat("5. Genes passing EXONIC filter:", sum(qc$passes_exonic_filter), "/", total_genes,
    "(", round(100 * sum(qc$passes_exonic_filter) / total_genes, 1), "%)\n")

cat("6. Final SELECTED genes:", final_selected, "/", total_genes,
    "(", round(100 * final_selected / total_genes, 1), "%)\n\n")

# Failure reasons summary
cat("=== Failure Reasons ===\n")
failure_counts <- table(qc$failure_reason)
failure_counts <- failure_counts[order(failure_counts, decreasing = TRUE)]

for(i in 1:length(failure_counts)) {
  reason <- names(failure_counts)[i]
  count <- failure_counts[i]
  pct <- round(100 * count / total_genes, 1)
  cat(sprintf("%-50s: %3d genes (%4.1f%%)\n", reason, count, pct))
}

cat("\n=== Peak Statistics (genes with peaks) ===\n")
genes_with_peaks <- qc[qc$has_macs2_peaks == TRUE]
if(nrow(genes_with_peaks) > 0) {
  cat("Raw peak counts per gene:\n")
  cat("  Mean:", round(mean(genes_with_peaks$peak_count_raw), 1), "\n")
  cat("  Median:", median(genes_with_peaks$peak_count_raw), "\n")
  cat("  Range:", min(genes_with_peaks$peak_count_raw), "-", max(genes_with_peaks$peak_count_raw), "\n\n")
  
  genes_with_best <- qc[!is.na(qc$best_peak_score)]
  if(nrow(genes_with_best) > 0) {
    cat("Best peak scores:\n")
    cat("  Mean:", round(mean(genes_with_best$best_peak_score), 1), "\n")
    cat("  Median:", round(median(genes_with_best$best_peak_score), 1), "\n")
    cat("  Range:", round(min(genes_with_best$best_peak_score), 1), "-", round(max(genes_with_best$best_peak_score), 1), "\n\n")
  }
  
  genes_with_exonic <- qc[!is.na(qc$exonic_fraction)]
  if(nrow(genes_with_exonic) > 0) {
    cat("Exonic fractions:\n")
    cat("  Mean:", round(mean(genes_with_exonic$exonic_fraction), 3), "\n")
    cat("  Median:", round(median(genes_with_exonic$exonic_fraction), 3), "\n")
    cat("  Range:", round(min(genes_with_exonic$exonic_fraction), 3), "-", round(max(genes_with_exonic$exonic_fraction), 3), "\n")
  }
}

# Write detailed report
cat("\n=== Writing detailed report to", opt$output_summary, "===\n")

# Create detailed output
detailed_output <- c(
  "=== MACS2 Peak Selection Detailed Report ===",
  paste("Generated:", Sys.time()),
  paste("Input file:", opt$qc_file),
  "",
  "=== Summary Statistics ===",
  paste("Total target genes:", total_genes),
  paste("Final selected genes:", final_selected),
  paste("Success rate:", round(100 * final_selected / total_genes, 1), "%"),
  "",
  "=== Genes by Failure Reason ===",
  "",
  sapply(names(failure_counts), function(reason) {
    genes <- qc$gene_id[qc$failure_reason == reason]
    paste0("## ", reason, " (", length(genes), " genes):")
  }),
  "",
  "=== Detailed Gene List ===",
  "",
  apply(qc, 1, function(row) {
    sprintf("%s\t%s\t%s\tpeaks:%d/%d/%d\tscore:%.1f\texonic:%.3f", 
            row["gene_id"], 
            row["failure_reason"],
            paste0(row["gene_chr"], ":", row["gene_start"], "-", row["gene_end"]),
            as.numeric(row["peak_count_raw"]),
            as.numeric(row["peak_count_significant"]), 
            as.numeric(row["peak_count_overlapping"]),
            if(is.na(row["best_peak_score"])) 0 else as.numeric(row["best_peak_score"]),
            if(is.na(row["exonic_fraction"])) 0 else as.numeric(row["exonic_fraction"])
    )
  })
)

writeLines(detailed_output, opt$output_summary)
cat("Report written successfully!\n")
