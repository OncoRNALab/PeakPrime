#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
})

option_list <- list(
  make_option("--gtf", type="character", help="Input GTF file"),
  make_option("--genes", type="character", help="Gene list file (one gene per line)"),
  make_option("--out", type="character", help="Output filtered GTF file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$gtf) || is.null(opt$genes) || is.null(opt$out)) {
  stop("All parameters (--gtf, --genes, --out) are required")
}

if (!file.exists(opt$gtf)) {
  stop(sprintf("GTF file not found: %s", opt$gtf))
}

if (!file.exists(opt$genes)) {
  stop(sprintf("Gene list file not found: %s", opt$genes))
}

# Load GTF (this may take 30-60 seconds for full genome)
cat("============================================================\n")
cat("GTF FILTERING FOR PLOTTING OPTIMIZATION\n")
cat("============================================================\n")
cat(sprintf("Loading GTF: %s\n", basename(opt$gtf)))
start_time <- Sys.time()
gtf <- import(opt$gtf)
load_time <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat(sprintf("✓ Loaded %d features in %.1f seconds\n", length(gtf), load_time))

# Read gene list
genes <- readLines(opt$genes)
genes <- genes[nchar(genes) > 0]  # Remove empty lines
genes <- unique(genes)  # Remove duplicates
cat(sprintf("✓ Target genes: %d\n", length(genes)))

# Determine which gene column to use
gene_col <- NULL
if ("gene_id" %in% names(mcols(gtf))) {
  gene_col <- "gene_id"
} else if ("gene_name" %in% names(mcols(gtf))) {
  gene_col <- "gene_name"
} else if ("gene" %in% names(mcols(gtf))) {
  gene_col <- "gene"
} else {
  stop("GTF must have one of these columns: gene_id, gene_name, or gene")
}
cat(sprintf("✓ Using column: %s\n", gene_col))

# Filter GTF to target genes
cat("Filtering GTF...\n")
filtered <- gtf[mcols(gtf)[[gene_col]] %in% genes]

# Check if we got any results
if (length(filtered) == 0) {
  warning("No features matched! Check that your gene IDs match the GTF format.")
  cat("\nDiagnostic information:\n")
  cat(sprintf("  First gene in your list: %s\n", genes[1]))
  cat(sprintf("  First gene_id in GTF: %s\n", as.character(mcols(gtf)[[gene_col]])[1]))
  cat(sprintf("  Gene ID format mismatch? Check if you're using symbols vs Ensembl IDs\n"))
}

# Export filtered GTF
export(filtered, opt$out, format="gtf")
cat("============================================================\n")
cat(sprintf("✓ Filtered GTF: %d -> %d features (%.1f%% reduction)\n", 
            length(gtf), length(filtered), 
            100 * (1 - length(filtered)/length(gtf))))
cat(sprintf("✓ Output file: %s\n", opt$out))
cat("============================================================\n")

# Calculate expected speedup
if (length(filtered) > 0) {
  speedup <- length(gtf) / length(filtered)
  cat(sprintf("\n🚀 Expected GTF loading speedup: %.0f×\n\n", speedup))
}
