#!/usr/bin/env Rscript
library(Biostrings)
library(Rsamtools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: extract_seq.R <fasta> <chrom> <start> <end>")
fasta <- args[1]
chrom <- args[2]
start <- as.integer(args[3])
end <- as.integer(args[4])
fa <- FaFile(fasta)
open(fa)
seq <- scanFa(fa, GRanges(chrom, IRanges(start, end)))
close(fa)
cat(as.character(seq), "\n")