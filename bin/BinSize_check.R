#!/usr/bin/env Rscript
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: BinSize_check.R <bigwig_file>")
bw_file <- args[1]
seqs <- seqnames(seqinfo(BigWigFile(bw_file)))
found <- FALSE
for (chr in seqs) {
  gr <- GRanges(chr, IRanges(1, 1e6))
  bw <- import(bw_file, which = gr)
  if (length(bw) >= 2) {
    bin_sizes <- diff(start(bw))
    cat("Chromosome:", chr, "\n")
    cat("Bin size(s):", paste(unique(bin_sizes), collapse = ", "), "\n")
    found <- TRUE
    break
  }
}
if (!found) stop("No chromosome in BigWig file contains enough data to detect bin size.")
