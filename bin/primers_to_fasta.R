#!/usr/bin/env Rscript

# Convert primer TSV to FASTA format for transcriptome alignment
suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

opt_list <- list(
  make_option("--primers_tsv", type="character", help="Input primers TSV file (cdna_primers.tsv)"),
  make_option("--out_fasta", type="character", default="primers_for_alignment.fa", help="Output FASTA file"),
  make_option("--max_primers_per_gene", type="integer", default=3, help="Maximum number of primers per gene to include (default: 3 best)")
)

opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$primers_tsv)) {
  stop("--primers_tsv is required")
}

# Read primers data
primers <- read.delim(opt$primers_tsv, stringsAsFactors = FALSE)

cat("Input primers:", nrow(primers), "from", length(unique(primers$gene_id)), "genes\n")

# Limit to max_primers_per_gene (take first N primers per gene, assuming they're sorted by quality)
if (!is.na(opt$max_primers_per_gene) && opt$max_primers_per_gene > 0) {
  primers_filtered <- do.call(rbind, lapply(split(primers, primers$gene_id), function(gene_primers) {
    head(gene_primers, opt$max_primers_per_gene)
  }))
  primers <- primers_filtered
  cat("Filtered to:", nrow(primers), "primers (max", opt$max_primers_per_gene, "per gene)\n")
}

# Create FASTA sequences
sequences <- DNAStringSet(primers$primer_sequence)

# Create descriptive FASTA headers
# Format: >gene_id|primer_index|primer_type|strand
headers <- paste0(
  primers$gene_id, "|",
  "idx", primers$primer_index, "|",
  primers$primer_type, "|",
  "strand", primers$gene_strand
)

names(sequences) <- headers

# Write FASTA file
writeXStringSet(sequences, filepath = opt$out_fasta)

cat("FASTA file written:", opt$out_fasta, "\n")
cat("Ready for transcriptome alignment with bowtie2\n")

# Print summary
strand_summary <- table(primers$gene_strand, primers$primer_type)
cat("\nPrimer summary:\n")
print(strand_summary)

# Print example command for bowtie2 alignment
cat("\nExample bowtie2 alignment command:\n")
cat("# First, build transcriptome index (if not already done):\n")
cat("# bowtie2-build transcriptome.fa transcriptome_idx\n\n")
cat("# Then align primers:\n")
cat("bowtie2 -x transcriptome_idx -f -U", opt$out_fasta, "-S primers_alignment.sam --very-sensitive-local -k 10\n\n")
cat("# Convert to sorted BAM for analysis:\n")
cat("samtools view -bS primers_alignment.sam | samtools sort -o primers_alignment.bam\n")
cat("samtools index primers_alignment.bam\n")
