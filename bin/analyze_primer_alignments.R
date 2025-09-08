#!/usr/bin/env Rscript

# Analyze primer alignment results from bowtie2
suppressPackageStartupMessages({
  library(optparse)
  library(Rsamtools)
  library(GenomicAlignments)
  library(Biostrings)
})

opt_list <- list(
  make_option("--alignment_bam", type="character", help="Input BAM file from bowtie2 alignment"),
  make_option("--primers_tsv", type="character", help="Original primers TSV file"),
  make_option("--transcriptome_fasta", type="character", default=NULL, help="Transcriptome FASTA file for gene name mapping"),
  make_option("--transcript_mapping", type="character", default=NULL, help="Pre-built transcript-to-gene mapping TSV file (faster alternative to FASTA)"),
  make_option("--out_report", type="character", default="primer_alignment_report.tsv", help="Output alignment report"),
  make_option("--out_summary", type="character", default="primer_alignment_summary.tsv", help="Detailed alignment summary per primer")
)

opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$alignment_bam) || is.null(opt$primers_tsv)) {
  stop("Both --alignment_bam and --primers_tsv are required")
}

# Function to parse gene names from transcriptome FASTA headers
parse_gene_names <- function(transcriptome_fasta = NULL, transcript_mapping = NULL) {
  
  if (!is.null(transcript_mapping) && file.exists(transcript_mapping)) {
    # Use pre-built mapping file (much faster)
    cat("Loading pre-built transcript-to-gene mapping...\n")
    mapping_data <- read.table(transcript_mapping, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    
    if (!"transcript_id" %in% colnames(mapping_data) || !"gene_name" %in% colnames(mapping_data)) {
      stop("Mapping file must contain 'transcript_id' and 'gene_name' columns")
    }
    
    # Convert to list format for compatibility
    transcript_to_gene <- as.list(setNames(mapping_data$gene_name, mapping_data$transcript_id))
    cat(sprintf("Loaded mapping for %d transcripts\n", length(transcript_to_gene)))
    return(transcript_to_gene)
  }
  
  # Fallback to FASTA parsing if no mapping file
  if (is.null(transcriptome_fasta) || !file.exists(transcriptome_fasta)) {
    cat("Warning: Transcriptome FASTA not provided or not found. Gene names will not be available.\n")
    return(NULL)
  }
  
  cat("Parsing gene names from transcriptome FASTA (this may take time)...\n")
  tryCatch({
    # Read FASTA sequences
    fasta_seqs <- readDNAStringSet(transcriptome_fasta)
    transcript_to_gene <- list()
    
    # Parse each header to extract gene name
    for (i in seq_along(fasta_seqs)) {
      full_header <- names(fasta_seqs)[i]
      
      # Extract transcript ID (first part before space)
      transcript_id <- strsplit(full_header, " ")[[1]][1]
      
      # Extract gene name from header (assuming format like "ENST00000456328 gene=DDX11L2")
      gene_name <- "Unknown"
      if (grepl("gene=", full_header)) {
        # Use more robust regex to extract gene name
        gene_match <- regmatches(full_header, regexpr("gene=([A-Za-z0-9_.-]+)", full_header))
        if (length(gene_match) > 0) {
          gene_name <- sub("gene=", "", gene_match)
        }
      }
      
      # Store mapping using clean transcript ID
      transcript_to_gene[[transcript_id]] <- gene_name
      
      # Progress indicator for large files
      if (i %% 10000 == 0) {
        cat(sprintf("Processed %d/%d sequences\n", i, length(fasta_seqs)))
      }
    }
    
    cat("Parsed", length(transcript_to_gene), "transcript-to-gene mappings\n")
    
    # Debug: check if a few specific transcripts are in the mapping
    test_transcripts <- c("ENST00000521262", "ENST00000524349", "ENST00000371437")
    found_count <- 0
    for (test_id in test_transcripts) {
      if (test_id %in% names(transcript_to_gene)) {
        found_count <- found_count + 1
        if (found_count <= 2) {
          cat("  Example mapping:", test_id, "->", transcript_to_gene[[test_id]], "\n")
        }
      }
    }
    
    return(transcript_to_gene)
    
  }, error = function(e) {
    cat("Error parsing transcriptome FASTA:", e$message, "\n")
    return(NULL)
  })
}

# Read original primer data
primers <- read.delim(opt$primers_tsv, stringsAsFactors = FALSE)
primers$primer_id <- paste0(
  primers$gene_id, "|",
  "idx", primers$primer_index, "|",
  primers$primer_type, "|",
  "strand", primers$gene_strand
)

# Parse transcriptome FASTA for gene name mapping
transcript_to_gene <- parse_gene_names(opt$transcriptome_fasta, opt$transcript_mapping)

# Read alignment results
bam_file <- BamFile(opt$alignment_bam)

# Try to read alignments with query names
tryCatch({
  # Use scanBam to get alignments (exclude MAPQ since it's meaningless with -a flag)
  param <- ScanBamParam(what = c("qname", "rname", "pos", "qwidth", "strand"))
  bam_data <- scanBam(bam_file, param = param)[[1]]
  
  cat("Read", length(bam_data$qname), "alignments from BAM file\n")
  
  # Handle case where there are no alignments
  if (length(bam_data$qname) == 0) {
    cat("Warning: No alignments found in BAM file\n")
    alignment_df <- data.frame(
      primer_id = character(0),
      target = character(0),
      start = integer(0),
      end = integer(0),
      width = integer(0),
      strand = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    # Create alignment dataframe from scanBam results
    alignment_df <- data.frame(
      primer_id = bam_data$qname,
      target = as.character(bam_data$rname),
      start = bam_data$pos,
      end = bam_data$pos + bam_data$qwidth - 1,
      width = bam_data$qwidth,
      strand = as.character(bam_data$strand),
      stringsAsFactors = FALSE
    )
    
    # Remove any rows with NA values
    alignment_df <- alignment_df[!is.na(alignment_df$primer_id), ]
  }
  
}, error = function(e) {
  cat("Error reading BAM file with scanBam, trying readGAlignments...\n")
  
  # Fallback to original method
  alignments <- readGAlignments(bam_file)
  cat("Read", length(alignments), "alignments from BAM file\n")
  
  # Handle case where there are no alignments
  if (length(alignments) == 0) {
    cat("Warning: No alignments found in BAM file\n")
    alignment_df <<- data.frame(
      primer_id = character(0),
      target = character(0),
      start = integer(0),
      end = integer(0),
      width = integer(0),
      strand = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    # Extract alignment statistics
    primer_names <- names(alignments)
    if (is.null(primer_names)) {
      primer_names <- paste0("read_", seq_along(alignments))
    }
    
    alignment_df <<- data.frame(
      primer_id = primer_names,
      target = as.character(seqnames(alignments)),
      start = start(alignments),
      end = end(alignments),
      width = width(alignments),
      strand = as.character(strand(alignments)),
      stringsAsFactors = FALSE
    )
  }
})

# Count alignments per primer
alignment_counts <- table(alignment_df$primer_id)

# Create comprehensive report
report <- data.frame(
  gene_id = primers$gene_id,
  primer_index = primers$primer_index,
  primer_type = primers$primer_type,
  primer_sequence = primers$primer_sequence,
  gene_strand = primers$gene_strand,
  primer_id = primers$primer_id,
  num_alignments = as.numeric(alignment_counts[primers$primer_id]),
  stringsAsFactors = FALSE
)

# Replace NA with 0 for primers with no alignments
report$num_alignments[is.na(report$num_alignments)] <- 0

# Write report (no alignment quality or MAPQ columns)
write.table(report, opt$out_report, sep = "\t", quote = FALSE, row.names = FALSE)

# Create detailed alignment summary (gene, primer, strand, transcript alignments)
if (nrow(alignment_df) > 0) {
  # Use all alignments since MAPQ is meaningless with bowtie2 -a
  cat("Processing", nrow(alignment_df), "total alignments\n")
  
  # Merge alignment data with primer information
  detailed_summary <- merge(alignment_df, primers[, c("primer_id", "gene_id", "primer_index", 
                                                       "primer_type", "primer_sequence", "gene_strand")], 
                           by = "primer_id", all.x = TRUE)
  
  if (nrow(detailed_summary) > 0) {
    # Add gene names from transcriptome mapping
    if (!is.null(transcript_to_gene)) {
      detailed_summary$aligned_gene_name <- sapply(detailed_summary$target, function(transcript_id) {
        gene_name <- transcript_to_gene[[transcript_id]]
        if (is.null(gene_name)) gene_name <- "Unknown"
        return(gene_name)
      })
    } else {
      detailed_summary$aligned_gene_name <- "Unknown"
    }
    
    # Reorder columns for better readability (remove MAPQ)
    detailed_summary <- detailed_summary[, c("gene_id", "primer_index", "primer_type", 
                                            "primer_sequence", "gene_strand", "target", 
                                            "aligned_gene_name", "start", "end", "width", 
                                            "strand")]
    
    # Sort by gene, primer index, and number of mismatches (best alignments first)
    detailed_summary <- detailed_summary[order(detailed_summary$gene_id, 
                                              detailed_summary$primer_index), ]
    
    # Rename columns for clarity
    names(detailed_summary)[names(detailed_summary) == "target"] <- "aligned_transcript"
    names(detailed_summary)[names(detailed_summary) == "start"] <- "alignment_start"
    names(detailed_summary)[names(detailed_summary) == "end"] <- "alignment_end"
    names(detailed_summary)[names(detailed_summary) == "width"] <- "alignment_length"
    names(detailed_summary)[names(detailed_summary) == "strand"] <- "alignment_strand"
    
    # Write detailed summary
    write.table(detailed_summary, opt$out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("Detailed alignment summary written to:", opt$out_summary, "\n")
  } else {
    # Create empty detailed summary file  
    empty_summary <- data.frame(
      gene_id = character(0),
      primer_index = integer(0),
      primer_type = character(0),
      primer_sequence = character(0),
      gene_strand = character(0),
      aligned_transcript = character(0),
      aligned_gene_name = character(0),
      alignment_start = integer(0),
      alignment_end = integer(0),
      alignment_length = integer(0),
      alignment_strand = character(0)
    )
    write.table(empty_summary, opt$out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("No high-quality alignments found - empty detailed summary written to:", opt$out_summary, "\n")
  }
} else {
  # Create empty detailed summary file
  empty_summary <- data.frame(
    gene_id = character(0),
    primer_index = integer(0),
    primer_type = character(0),
    primer_sequence = character(0),
    gene_strand = character(0),
    aligned_transcript = character(0),
    aligned_gene_name = character(0),
    alignment_start = integer(0),
    alignment_end = integer(0),
    alignment_length = integer(0),
    alignment_strand = character(0)
  )
  write.table(empty_summary, opt$out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("No alignments found - empty detailed summary written to:", opt$out_summary, "\n")
}

cat("Alignment analysis complete\n")
cat("Report written to:", opt$out_report, "\n")

# Print summary statistics
cat("=== PRIMER ALIGNMENT SUMMARY ===\n")
cat("Total primers analyzed:", nrow(report), "\n")

cat("\nAlignment count distribution:\n")
alignment_summary <- table(cut(report$num_alignments, 
                              breaks = c(-1, 0, 1, 5, 20, Inf), 
                              labels = c("No hits", "1 hit", "2-5 hits", "6-20 hits", ">20 hits")))
for (bin in names(alignment_summary)) {
  cat("  ", bin, ":", alignment_summary[bin], "\n")
}

# Identify problematic primers
problematic <- report[report$alignment_quality %in% c("POOR", "FAIL"), ]
if (nrow(problematic) > 0) {
  cat("\n=== PROBLEMATIC PRIMERS ===\n")
  cat("Primers with poor/no alignment (consider excluding):\n")
  for (i in 1:min(10, nrow(problematic))) {  # Show max 10
    p <- problematic[i, ]
    cat("  ", p$gene_id, "primer", p$primer_index, "(", p$primer_type, ") -", p$num_alignments, "hits\n")
  }
  if (nrow(problematic) > 10) {
    cat("  ... and", nrow(problematic) - 10, "more (see full report)\n")
  }
}

# Identify best primers
best_primers <- report[report$alignment_quality == "PERFECT", ]
if (nrow(best_primers) > 0) {
  cat("\n=== BEST PRIMERS (PERFECT ALIGNMENT) ===\n")
  cat("Primers with unique transcriptome hits:\n")
  for (i in 1:min(5, nrow(best_primers))) {  # Show max 5
    p <- best_primers[i, ]
    cat("  ", p$gene_id, "primer", p$primer_index, "(", p$primer_type, ") - MAPQ:", p$best_mapq, "\n")
  }
  if (nrow(best_primers) > 5) {
    cat("  ... and", nrow(best_primers) - 5, "more perfect primers\n")
  }
}
