#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

opt_list <- list(
  make_option("--bw", type="character"),
  make_option("--gtf", type="character"),
  make_option("--genes", type="character"),
  make_option("--genome", type="character", default="BSgenome.Hsapiens.UCSC.hg38"),
  make_option("--fasta", type="character", default=NULL),
  make_option("--pad", type="integer", default=1000),
  make_option("--smooth_k", type="integer", default=31),
  make_option("--sliding_window", action="store_true", default=FALSE),
  make_option("--min_window_mean", type="double", default=NA_real_),
  make_option("--min_window_mean_pct", type="double", default=NA_real_),
  make_option("--max_gap", type="integer", default=NA_integer_),
  make_option("--search_slop", type="integer", default=1000),
  make_option("--trim_to_exon", action="store_true", default=FALSE),
  make_option("--trim_low_coverage_pct", type="double", default=NA_real_),
  make_option("--min_exonic_fraction", type="double", default=NA_real_),
  make_option("--out_fa", type="character", default="primer_targets.fa"),
  make_option("--out_bed", type="character", default="primer_targets.bed"),
  make_option("--out_peaks", type="character", default="peaks.tsv"),
  make_option("--out_qc", type="character", default="qc_coverage_summary.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Load gene annotation
cat("Loading GTF and gene list...\n")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
g_all <- genes(txdb)
ids <- unique(scan(opt$genes, what=character(), quiet=TRUE))
sel <- names(g_all) %in% ids
g <- g_all[sel]

ex_by_gene <- exonsBy(txdb, by="gene")
ex_by_gene <- ex_by_gene[names(ex_by_gene) %in% names(g)]
ex_merged <- endoapply(ex_by_gene, function(x) reduce(trim(x)))


# Output containers
peak_rows <- list()
qc_rows <- list()

for (gid in names(g)) {
  cat("Processing gene:", gid, "\n")
  blocks <- ex_merged[[gid]]
  gene_chr <- as.character(seqnames(blocks)[1])
  gene_start <- min(start(blocks))
  gene_end <- max(end(blocks))
  pad <- opt$pad
  span_start <- max(1, gene_start - pad)
  span_end <- gene_end + pad
  gr <- GRanges(seqnames=gene_chr, ranges=IRanges(span_start, span_end))

  # Import coverage for the full span
  cov <- import(opt$bw, which=gr)
  region_length <- span_end - span_start + 1
  cov_rle <- Rle(0, region_length)
  for (j in seq_along(cov)) {
    interval <- cov[j]
    rle_start <- start(interval) - span_start + 1
    rle_end <- end(interval) - span_start + 1
    rle_start <- max(1, rle_start)
    rle_end <- min(region_length, rle_end)
    if (rle_start <= rle_end && rle_start <= region_length) {
      cov_rle[rle_start:rle_end] <- interval$score
    }
  }
  cov_vec <- as.numeric(cov_rle)
  if (!length(cov_vec) || max(cov_vec, na.rm=TRUE) == 0) {
    cat("  Skipping gene", gid, "- no coverage\n")
    next
  }

  # Compute dynamic threshold if percentage-based
  dynamic_min_mean <- opt$min_window_mean
  if (!is.na(opt$min_window_mean_pct)) {
    gene_peak_cov <- max(cov_vec, na.rm=TRUE)
    dynamic_min_mean <- gene_peak_cov * (opt$min_window_mean_pct / 100)
    cat("  Gene", gid, "peak coverage:", gene_peak_cov, "-> dynamic threshold:", dynamic_min_mean, "(", opt$min_window_mean_pct, "%)\n")
  }

  # Initial window centered at peak
  sm <- cov_vec
  max_idx <- which.max(sm)
  peak_pos <- span_start + max_idx - 1
  win_start <- max(span_start, peak_pos - pad)
  win_end <- min(span_end, peak_pos + pad)
  wcov <- cov_vec[(win_start - span_start + 1):(win_end - span_start + 1)]
  w_mean <- mean(wcov)
  w_lzr <- {
    r <- rle(wcov == 0)
    zeros <- r$lengths[r$values]
    if (length(zeros) == 0) 0 else max(zeros)
  }
  pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
  pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)
  cat("    window_mean =", w_mean, "(threshold:", dynamic_min_mean, ")\n")
  cat("    longest_zero_run =", w_lzr, "(threshold:", opt$max_gap, ")\n")
  cat("    pass_min_mean =", pass_min_mean, ", pass_max_gap =", pass_max_gap, "\n")

  # Sliding window search if QC fails
  strategy <- "peak_centered"
  if (!pass_min_mean || !pass_max_gap) {
    cat("    Triggering sliding window search...\n")
    best_mean <- -Inf
    best_center <- NA_integer_
    wsize <- 2 * pad + 1
    n <- length(cov_vec)
    for (i in 0:(n - wsize)) {
      j <- i + wsize
      win <- cov_vec[(i+1):j]
      meanv <- mean(win)
      r <- rle(win == 0)
      zeros <- r$lengths[r$values]
      lzr <- if (length(zeros) == 0) 0 else max(zeros)
      if (!is.na(dynamic_min_mean) && meanv < dynamic_min_mean) next
      if (!is.na(opt$max_gap) && lzr > opt$max_gap) next
      if (meanv > best_mean) {
        best_mean <- meanv
        best_center <- i + floor(wsize/2) + 1
      }
    }
    if (!is.infinite(best_mean)) {
      strategy <- "sliding_best"
      peak_pos <- span_start + best_center - 1
      win_start <- max(span_start, peak_pos - pad)
      win_end <- min(span_end, peak_pos + pad)
      wcov <- cov_vec[(win_start - span_start + 1):(win_end - span_start + 1)]
      w_mean <- mean(wcov)
      r <- rle(wcov == 0)
      zeros <- r$lengths[r$values]
      w_lzr <- if (length(zeros) == 0) 0 else max(zeros)
      pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
      pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)
      cat("    Found better window at position", peak_pos, "with mean coverage", w_mean, "\n")
    } else {
      cat("    No better window found, keeping original\n")
    }
  }

  # Trimming window ends based on coverage threshold
  trim_applied <- FALSE
  if (!is.na(opt$trim_low_coverage_pct)) {
    window_peak <- max(wcov)
    threshold <- window_peak * (opt$trim_low_coverage_pct / 100)
    start_trim <- 0
    for (i in seq_along(wcov)) {
      if (wcov[i] >= threshold) break
      start_trim <- start_trim + 1
    }
    end_trim <- 0
    for (i in rev(seq_along(wcov))) {
      if (wcov[i] >= threshold) break
      end_trim <- end_trim + 1
    }
    trimmed_length <- length(wcov) - start_trim - end_trim
    if (trimmed_length > 0 && (start_trim > 0 || end_trim > 0)) {
      trim_applied <- TRUE
      win_start <- win_start + start_trim
      win_end <- win_end - end_trim
      peak_pos <- as.integer((win_start + win_end) / 2)
      wcov <- cov_vec[(win_start - span_start + 1):(win_end - span_start + 1)]
      w_mean <- mean(wcov)
      r <- rle(wcov == 0)
      zeros <- r$lengths[r$values]
      w_lzr <- if (length(zeros) == 0) 0 else max(zeros)
      pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
      pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)
      strategy <- paste0(strategy, "+trimmed")
      cat("    Trimmed window:", win_start, "-", win_end, "(", win_end - win_start + 1, "bp)\n")
    }
  }

  # Exonic fraction calculation
  window_gr <- GRanges(seqnames=gene_chr, ranges=IRanges(win_start, win_end), strand=strand(g[gid]))
  cat("  Window GRanges: chr=", as.character(seqnames(window_gr)), " start=", start(window_gr), " end=", end(window_gr), " strand=", as.character(strand(window_gr)), "\n")
  cat("  Exonic blocks: chr=", as.character(seqnames(blocks)), " start=", start(blocks), " end=", end(blocks), " strand=", as.character(strand(blocks)), "\n")
  hits <- findOverlaps(window_gr, blocks)
  exonic_bases <- 0
  if (length(hits) > 0) {
    overlapping_blocks <- blocks[subjectHits(hits)]
    ov <- intersect(window_gr, overlapping_blocks)
    exonic_bases <- sum(width(ov))
  }
  cat("  Exonic bases found:", exonic_bases, " out of window width ", width(window_gr), "\n")
  exonic_fraction <- exonic_bases / width(window_gr)
  fail_exonic_fraction <- !is.na(opt$min_exonic_fraction) && exonic_fraction < opt$min_exonic_fraction

  # Exon boundary trimming
  trimmed_to_exon <- FALSE
  if (opt$trim_to_exon) {
    containing <- blocks[peak_pos >= start(blocks) & peak_pos <= end(blocks)]
    if (length(containing) >= 1) {
      exon_range <- containing[1]
      new_start <- max(win_start, start(exon_range))
      new_end <- min(win_end, end(exon_range))
      window_gr <- GRanges(seqnames=gene_chr, ranges=IRanges(new_start, new_end))
      trimmed_to_exon <- TRUE
      # Recompute exonic fraction
      hits <- findOverlaps(window_gr, blocks)
      if (length(hits) > 0) {
        overlapping_blocks <- blocks[subjectHits(hits)]
        ov <- intersect(window_gr, overlapping_blocks)
        exonic_bases <- sum(width(ov))
      } else {
        exonic_bases <- 0
      }
      exonic_fraction <- exonic_bases / width(window_gr)
    }
  }

  # Output rows (without exonic_fraction, trimmed_to_exon, fail_exonic_fraction)
  peak_rows[[gid]] <- data.frame(
    gene = gid,
    chr = gene_chr,
    pos = peak_pos,
    start = start(window_gr),
    end = end(window_gr),
    strand = as.character(strand(g[gid])),
    strategy = strategy,
    stringsAsFactors = FALSE
  )
  qc_rows[[gid]] <- data.frame(
    gene = gid,
    max_cov = max(cov_vec),
    window_mean = w_mean,
    window_zeros = sum(wcov == 0),
    longest_zero_run = w_lzr,
    pass_min_mean = pass_min_mean,
    pass_max_gap = pass_max_gap,
    strategy = strategy,
    stringsAsFactors = FALSE
  )
}



# After all genes processed, recalculate exonic_fraction, trimmed_to_exon, fail_exonic_fraction for all windows
peaks_df <- do.call(rbind, peak_rows)
qc_df <- do.call(rbind, qc_rows)
gr_peaks <- GRanges(peaks_df$chr, IRanges(peaks_df$start, peaks_df$end), strand = peaks_df$strand, gene = peaks_df$gene)
exonic_fraction <- rep(NA_real_, length(gr_peaks))
trimmed <- rep(FALSE, length(gr_peaks))
fail_exonic_fraction <- rep(FALSE, length(gr_peaks))
for (i in seq_along(gr_peaks)) {
  gid <- peaks_df$gene[i]
  blocks <- ex_merged[[gid]]
  w <- gr_peaks[i]
  hits <- findOverlaps(w, blocks)
  if (length(hits) > 0) {
    overlapping_blocks <- blocks[subjectHits(hits)]
    ov <- intersect(w, overlapping_blocks)
    exonic_bases <- sum(width(ov))
    exonic_fraction[i] <- exonic_bases / width(w)
  } else {
    exonic_fraction[i] <- 0
  }
  if (!is.na(opt$min_exonic_fraction) && exonic_fraction[i] < opt$min_exonic_fraction) {
    fail_exonic_fraction[i] <- TRUE
  }
  if (opt$trim_to_exon) {
    peak_pos <- peaks_df$pos[i]
    containing <- blocks[peak_pos >= start(blocks) & peak_pos <= end(blocks)]
    if (length(containing) >= 1) {
      exon_range <- containing[1]
      current_start <- start(gr_peaks[i])
      current_end <- end(gr_peaks[i])
      new_start <- max(current_start, start(exon_range))
      new_end <- min(current_end, end(exon_range))
      gr_peaks[i] <- GRanges(seqnames(exon_range), IRanges(new_start, new_end), strand = strand(exon_range), gene = gid)
      trimmed[i] <- TRUE
      w_new <- gr_peaks[i]
      hits_new <- findOverlaps(w_new, blocks)
      if (length(hits_new) > 0) {
        overlapping_blocks_new <- blocks[subjectHits(hits_new)]
        ov_new <- intersect(w_new, overlapping_blocks_new)
        exonic_bases_new <- sum(width(ov_new))
        exonic_fraction[i] <- exonic_bases_new / width(w_new)
      } else {
        exonic_fraction[i] <- 0
      }
    }
  }
}
peaks_df$exonic_fraction <- exonic_fraction
peaks_df$trimmed_to_exon <- trimmed
peaks_df$fail_exonic_fraction <- fail_exonic_fraction
qc_df$exonic_fraction <- exonic_fraction[match(qc_df$gene, peaks_df$gene)]
qc_df$fail_exonic_fraction <- fail_exonic_fraction[match(qc_df$gene, peaks_df$gene)]
qc_df$trimmed_to_exon <- trimmed[match(qc_df$gene, peaks_df$gene)]
write.table(peaks_df, file = opt$out_peaks, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc_df, file = opt$out_qc, sep = "\t", quote = FALSE, row.names = FALSE)

# Extract sequences from FASTA or BSgenome
library(Biostrings)
library(Rsamtools)
seqs <- NULL
fasta_path <- if (!is.null(opt$fasta)) opt$fasta else NULL
if (!is.null(fasta_path)) {
  cat("Extracting sequences from FASTA...\n")
  faf <- Rsamtools::FaFile(fasta_path)
  Rsamtools::open.FaFile(faf)
  on.exit(Rsamtools::close.FaFile(faf), add = TRUE)
  gr_peaks <- GRanges(peaks_df$chr, IRanges(peaks_df$start, peaks_df$end), strand = peaks_df$strand)
  seqs <- Rsamtools::getSeq(faf, gr_peaks)
} else {
  cat("No FASTA provided, skipping sequence extraction.\n")
}

# Append universal sequence to 3' end
universal_seq <- ""
if (!is.null(seqs)) {
  seqs_appended <- paste0(as.character(seqs), universal_seq)
  seqs_appended <- DNAStringSet(seqs_appended)
  nm <- paste0(peaks_df$gene, "|", peaks_df$chr, ":", peaks_df$start, "-", peaks_df$end, "(", peaks_df$strand, ")")
  names(seqs_appended) <- nm
  writeXStringSet(seqs_appended, filepath = opt$out_fa)
}

# Write BED file
bed_df <- data.frame(
  chr = peaks_df$chr,
  start = peaks_df$start - 1, # BED is 0-based
  end = peaks_df$end,
  name = peaks_df$gene,
  score = 0,
  strand = peaks_df$strand
)
write.table(bed_df, file = opt$out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
