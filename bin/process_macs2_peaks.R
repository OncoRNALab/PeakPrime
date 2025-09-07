#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  
  # Load required packages with better error messages
  tryCatch({
    library(rtracklayer)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(Biostrings)
    library(data.table)
  }, error = function(e) {
    cat("Error loading required packages:", conditionMessage(e), "\n")
    quit("no", status=1)
  })
})

opt_list <- list(
  make_option("--narrowpeak", type="character", help="MACS2 narrowPeak file"),
  make_option("--gtf", type="character", help="GTF annotation file"),
  make_option("--genes", type="character", help="Target genes file"),
  make_option("--fasta", type="character", default="NO_FILE", help="Genome FASTA file"),
  make_option("--pvalue_threshold", type="double", default=0.05, help="p-value threshold for peak filtering"),
  make_option("--min_peak_score", type="double", default=0, help="Minimum peak score threshold"),
  make_option("--pad", type="integer", default=100, help="[UNUSED] Legacy parameter - target regions now use exact peak boundaries"),
  make_option("--min_exonic_fraction", type="character", default="NA", help="Minimum exonic fraction (or NA to disable)"),
  make_option("--trim_to_exon", type="character", default="false", help="Trim peaks to exon boundaries"),
  make_option("--out_fa", type="character", help="Output FASTA file"),
  make_option("--out_bed", type="character", help="Output BED file"),
  make_option("--out_peaks", type="character", help="Output peaks TSV file"),
  make_option("--out_qc", type="character", help="Output QC summary file")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Parse parameters
min_exonic_fraction <- if(opt$min_exonic_fraction == "NA") NA_real_ else as.numeric(opt$min_exonic_fraction)
trim_to_exon <- opt$trim_to_exon == "true"

cat("Loading annotation and gene list...\n")

# Load GTF and create gene/exon annotations
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
g_all <- genes(txdb)
target_genes <- unique(scan(opt$genes, what=character(), quiet=TRUE))
sel <- names(g_all) %in% target_genes
g_target <- g_all[sel]

# Get exons by gene for overlap calculations
ex_by_gene <- exonsBy(txdb, by="gene")
ex_by_gene <- ex_by_gene[names(ex_by_gene) %in% names(g_target)]
ex_merged <- endoapply(ex_by_gene, function(x) reduce(trim(x)))

# Initialize comprehensive QC tracking for ALL target genes
all_genes_qc <- data.frame(
  gene_id = target_genes,
  gene_chr = as.character(seqnames(g_target))[match(target_genes, names(g_target))],
  gene_start = start(g_target)[match(target_genes, names(g_target))],
  gene_end = end(g_target)[match(target_genes, names(g_target))],
  gene_strand = as.character(strand(g_target))[match(target_genes, names(g_target))],
  has_macs2_peaks = FALSE,
  has_significant_peaks = FALSE,
  has_overlapping_peaks = FALSE,
  is_best_peak = FALSE,
  passes_exonic_filter = FALSE,
  final_selection = FALSE,
  failure_reason = "No MACS2 peaks detected",
  peak_count_raw = 0,
  peak_count_significant = 0,
  peak_count_overlapping = 0,
  best_peak_score = NA_real_,
  best_peak_pvalue = NA_real_,
  best_peak_qvalue = NA_real_,
  exonic_fraction = NA_real_,
  stringsAsFactors = FALSE
)

cat("Loading MACS2 peaks...\n")

# Read MACS2 narrowPeak file
# Columns: chr, start, end, name, score, strand, fold_change, pvalue, qvalue, summit_offset
peaks <- fread(opt$narrowpeak, header=FALSE, col.names=c("chr", "start", "end", "name", "score", "strand", "fold_change", "pvalue", "qvalue", "summit_offset"))

cat("Found", nrow(peaks), "total peaks\n")

# Track genes with any MACS2 peaks
if(nrow(peaks) > 0) {
  peaks_gr_all <- GRanges(
    seqnames = peaks$chr,
    ranges = IRanges(start = peaks$start + 1, end = peaks$end),
    strand = "*"
  )
  
  # Find which genes have any MACS2 peaks
  overlaps_all <- findOverlaps(peaks_gr_all, g_target)
  genes_with_any_peaks <- unique(names(g_target)[subjectHits(overlaps_all)])
  all_genes_qc$has_macs2_peaks[all_genes_qc$gene_id %in% genes_with_any_peaks] <- TRUE
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_with_any_peaks] <- "Has MACS2 peaks"
  
  # Count raw peaks per gene
  peak_counts <- table(names(g_target)[subjectHits(overlaps_all)])
  all_genes_qc$peak_count_raw[match(names(peak_counts), all_genes_qc$gene_id)] <- as.numeric(peak_counts)
}

# Filter peaks by p-value and score
significant_peaks <- peaks[pvalue >= -log10(opt$pvalue_threshold) & score >= opt$min_peak_score]
cat("Retained", nrow(significant_peaks), "significant peaks after filtering\n")

# Track genes with significant peaks
if(nrow(significant_peaks) > 0) {
  peaks_gr_sig <- GRanges(
    seqnames = significant_peaks$chr,
    ranges = IRanges(start = significant_peaks$start + 1, end = significant_peaks$end),
    strand = "*"
  )
  
  overlaps_sig <- findOverlaps(peaks_gr_sig, g_target)
  genes_with_sig_peaks <- unique(names(g_target)[subjectHits(overlaps_sig)])
  all_genes_qc$has_significant_peaks[all_genes_qc$gene_id %in% genes_with_sig_peaks] <- TRUE
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_with_sig_peaks] <- "Has significant peaks"
  
  # Count significant peaks per gene
  peak_counts_sig <- table(names(g_target)[subjectHits(overlaps_sig)])
  all_genes_qc$peak_count_significant[match(names(peak_counts_sig), all_genes_qc$gene_id)] <- as.numeric(peak_counts_sig)
  
  # Update failure reasons for genes that had peaks but lost them
  failed_at_significance <- all_genes_qc$has_macs2_peaks & !all_genes_qc$has_significant_peaks
  all_genes_qc$failure_reason[failed_at_significance] <- paste0("Failed significance filter (p>", opt$pvalue_threshold, " or score<", opt$min_peak_score, ")")
}

if(nrow(significant_peaks) == 0) {
  cat("No significant peaks found. Writing comprehensive QC report and creating empty output files.\n")
  # Write comprehensive QC before exiting
  write.table(all_genes_qc, file = opt$out_qc, sep = "\t", quote = FALSE, row.names = FALSE)
  # Create empty output files
  writeLines("", opt$out_fa)
  writeLines("", opt$out_bed)
  write.table(data.frame(), opt$out_peaks, sep="\t", quote=FALSE, row.names=FALSE)
  quit("no", status=0)
}

# Convert to GRanges
peaks_gr <- GRanges(
  seqnames = significant_peaks$chr,
  ranges = IRanges(start = significant_peaks$start + 1, end = significant_peaks$end), # Convert to 1-based
  strand = "*",
  name = significant_peaks$name,
  score = significant_peaks$score,
  fold_change = significant_peaks$fold_change,
  pvalue = significant_peaks$pvalue,
  qvalue = significant_peaks$qvalue,
  summit_offset = significant_peaks$summit_offset
)

cat("Finding overlaps with target genes...\n")

# Find overlaps between peaks and target genes
overlaps <- findOverlaps(peaks_gr, g_target)
peak_gene_map <- data.frame(
  peak_idx = queryHits(overlaps),
  gene_idx = subjectHits(overlaps),
  gene_id = names(g_target)[subjectHits(overlaps)],
  stringsAsFactors = FALSE
)

if(nrow(peak_gene_map) == 0) {
  cat("No peaks overlap with target genes. Writing comprehensive QC report and creating empty output files.\n")
  # Update failure reason for genes with significant peaks but no overlap
  no_overlap <- all_genes_qc$has_significant_peaks & !all_genes_qc$gene_id %in% peak_gene_map$gene_id
  all_genes_qc$failure_reason[no_overlap] <- "Significant peaks found but no overlap with gene boundaries"
  
  # Write comprehensive QC before exiting
  write.table(all_genes_qc, file = opt$out_qc, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", opt$out_fa)
  writeLines("", opt$out_bed)
  write.table(data.frame(), opt$out_peaks, sep="\t", quote=FALSE, row.names=FALSE)
  quit("no", status=0)
}

# Track genes with overlapping peaks
genes_with_overlaps <- unique(peak_gene_map$gene_id)
all_genes_qc$has_overlapping_peaks[all_genes_qc$gene_id %in% genes_with_overlaps] <- TRUE
all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_with_overlaps] <- "Has overlapping significant peaks"

# Count overlapping peaks per gene
peak_counts_overlap <- table(peak_gene_map$gene_id)
all_genes_qc$peak_count_overlapping[match(names(peak_counts_overlap), all_genes_qc$gene_id)] <- as.numeric(peak_counts_overlap)

# Update failure reason for genes that had significant peaks but no overlaps
failed_at_overlap <- all_genes_qc$has_significant_peaks & !all_genes_qc$has_overlapping_peaks
all_genes_qc$failure_reason[failed_at_overlap] <- "Significant peaks found but no overlap with gene boundaries"

cat("Found", nrow(peak_gene_map), "peak-gene overlaps for", length(unique(peak_gene_map$gene_id)), "genes\n")

# Select best peak per gene (highest score)
best_peaks <- peak_gene_map[order(peak_gene_map$gene_id, -peaks_gr[peak_gene_map$peak_idx]$score), ]
best_peaks <- best_peaks[!duplicated(best_peaks$gene_id), ]

# Track which genes get selected as best peaks
genes_with_best_peaks <- unique(best_peaks$gene_id)
all_genes_qc$is_best_peak[all_genes_qc$gene_id %in% genes_with_best_peaks] <- TRUE
all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_with_best_peaks] <- "Selected as best peak"

# Store best peak statistics
for(i in 1:nrow(best_peaks)) {
  gene_id <- best_peaks$gene_id[i]
  peak_idx <- best_peaks$peak_idx[i]
  qc_idx <- which(all_genes_qc$gene_id == gene_id)
  
  all_genes_qc$best_peak_score[qc_idx] <- peaks_gr[peak_idx]$score
  all_genes_qc$best_peak_pvalue[qc_idx] <- peaks_gr[peak_idx]$pvalue
  all_genes_qc$best_peak_qvalue[qc_idx] <- peaks_gr[peak_idx]$qvalue
}

cat("Selected", nrow(best_peaks), "best peaks (one per gene)\n")

# Process selected peaks
selected_peaks_gr <- peaks_gr[best_peaks$peak_idx]
selected_peaks_gr$gene_id <- best_peaks$gene_id

# Use exact peak boundaries as target regions (no padding)
summit_positions <- start(selected_peaks_gr) + selected_peaks_gr$summit_offset

# Create target regions using exact peak coordinates
target_regions <- GRanges(
  seqnames = seqnames(selected_peaks_gr),
  ranges = IRanges(start = start(selected_peaks_gr), end = end(selected_peaks_gr)),
  strand = strand(g_target[best_peaks$gene_idx]),
  gene_id = selected_peaks_gr$gene_id,
  peak_score = selected_peaks_gr$score,
  peak_pvalue = selected_peaks_gr$pvalue,
  peak_qvalue = selected_peaks_gr$qvalue,
  summit_pos = summit_positions
)

cat("Calculating exonic fractions...\n")

# Calculate exonic fraction for each target region
exonic_fractions <- numeric(length(target_regions))
trimmed_to_exon <- logical(length(target_regions))

for(i in seq_along(target_regions)) {
  gene_id <- target_regions$gene_id[i]
  
  # Get exons for this gene
  if(gene_id %in% names(ex_merged)) {
    exons <- ex_merged[[gene_id]]
    
    # Calculate overlap with exons
    overlaps_exon <- findOverlaps(target_regions[i], exons)
    if(length(overlaps_exon) > 0) {
      intersect_ranges <- intersect(target_regions[i], exons[subjectHits(overlaps_exon)])
      exonic_bases <- sum(width(intersect_ranges))
      exonic_fractions[i] <- exonic_bases / width(target_regions[i])
    } else {
      exonic_fractions[i] <- 0
    }
    
    # Trim to exon if requested
    if(trim_to_exon) {
      summit_pos <- target_regions$summit_pos[i]
      containing_exons <- exons[summit_pos >= start(exons) & summit_pos <= end(exons)]
      if(length(containing_exons) > 0) {
        # Use the first containing exon
        exon_range <- containing_exons[1]
        new_start <- max(start(target_regions[i]), start(exon_range))
        new_end <- min(end(target_regions[i]), end(exon_range))
        
        # Update the target region
        target_regions[i] <- GRanges(
          seqnames = seqnames(target_regions[i]),
          ranges = IRanges(start = new_start, end = new_end),
          strand = strand(target_regions[i]),
          gene_id = gene_id,
          peak_score = target_regions$peak_score[i],
          peak_pvalue = target_regions$peak_pvalue[i],
          peak_qvalue = target_regions$peak_qvalue[i],
          summit_pos = summit_pos
        )
        
        trimmed_to_exon[i] <- TRUE
        
        # Recalculate exonic fraction
        overlaps_exon_new <- findOverlaps(target_regions[i], exons)
        if(length(overlaps_exon_new) > 0) {
          intersect_ranges_new <- intersect(target_regions[i], exons[subjectHits(overlaps_exon_new)])
          exonic_fractions[i] <- sum(width(intersect_ranges_new)) / width(target_regions[i])
        } else {
          exonic_fractions[i] <- 0
        }
      }
    }
  } else {
    exonic_fractions[i] <- 0
  }
}

# Store exonic fractions in QC before filtering
for(i in 1:length(target_regions)) {
  gene_id <- target_regions$gene_id[i]
  qc_idx <- which(all_genes_qc$gene_id == gene_id)
  all_genes_qc$exonic_fraction[qc_idx] <- exonic_fractions[i]
}

target_regions$exonic_fraction <- exonic_fractions
target_regions$trimmed_to_exon <- trimmed_to_exon

# Track genes that pass exonic filter
all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% target_regions$gene_id] <- TRUE

# Filter by minimum exonic fraction if specified
if(!is.na(min_exonic_fraction)) {
  keep_idx <- exonic_fractions >= min_exonic_fraction
  genes_before_filter <- target_regions$gene_id
  target_regions <- target_regions[keep_idx]
  genes_after_filter <- target_regions$gene_id
  
  # Update QC for genes that pass exonic filter
  all_genes_qc$passes_exonic_filter <- FALSE  # Reset
  all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% genes_after_filter] <- TRUE
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_after_filter] <- "Passed all filters"
  
  # Update failure reason for genes that failed exonic filter
  failed_exonic <- setdiff(genes_before_filter, genes_after_filter)
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% failed_exonic] <- paste0("Failed exonic fraction filter (<", min_exonic_fraction, ")")
  
  cat("Retained", length(target_regions), "regions after exonic fraction filtering\n")
} else {
  # No exonic filtering - all genes with best peaks pass
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% target_regions$gene_id] <- "Passed all filters"
}

# Mark final selections
all_genes_qc$final_selection[all_genes_qc$gene_id %in% target_regions$gene_id] <- TRUE

cat("Extracting sequences...\n")

# Extract sequences
if(opt$fasta != "NO_FILE" && file.exists(opt$fasta)) {
  fasta_file <- Rsamtools::FaFile(opt$fasta)
  Rsamtools::open.FaFile(fasta_file)
  on.exit(Rsamtools::close.FaFile(fasta_file), add = TRUE)
  
  sequences <- Rsamtools::getSeq(fasta_file, target_regions)
  names(sequences) <- paste0(
    target_regions$gene_id, "|",
    seqnames(target_regions), ":",
    start(target_regions), "-",
    end(target_regions), "(",
    strand(target_regions), ")"
  )
  
  Biostrings::writeXStringSet(sequences, filepath = opt$out_fa)
} else {
  cat("No FASTA file provided or file does not exist. Skipping sequence extraction.\n")
  writeLines("", opt$out_fa)
}

# Create BED file
bed_df <- data.frame(
  chr = seqnames(target_regions),
  start = start(target_regions) - 1, # BED is 0-based
  end = end(target_regions),
  name = paste0(target_regions$gene_id, "_peak"),
  score = target_regions$peak_score,
  strand = strand(target_regions)
)

write.table(bed_df, file = opt$out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Create output DataFrame
peaks_df <- data.frame(
  gene = target_regions$gene_id,
  chr = seqnames(target_regions),
  start = start(target_regions),
  end = end(target_regions),
  strand = strand(target_regions),
  summit_pos = target_regions$summit_pos,
  peak_score = target_regions$peak_score,
  peak_pvalue = target_regions$peak_pvalue,
  peak_qvalue = target_regions$peak_qvalue,
  exonic_fraction = target_regions$exonic_fraction,
  trimmed_to_exon = target_regions$trimmed_to_exon,
  strategy = "macs2_peak_boundaries",
  stringsAsFactors = FALSE
)

# Write output files
write.table(peaks_df, file = opt$out_peaks, sep = "\t", quote = FALSE, row.names = FALSE)

# Write comprehensive QC report (all genes, not just selected ones)
write.table(all_genes_qc, file = opt$out_qc, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Process completed successfully!\n")
cat("Output files:\n")
cat("- FASTA:", opt$out_fa, "\n")
cat("- BED:", opt$out_bed, "\n")
cat("- Peaks TSV:", opt$out_peaks, "\n")
cat("- QC Summary:", opt$out_qc, "\n")
