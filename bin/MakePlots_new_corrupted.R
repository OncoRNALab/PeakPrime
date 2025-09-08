#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  l  # lock to first exon's seqname/strand for safety
  gene_chr    <- as.character(seqnames(gene_exons_all)[1])
  gene_strand <- as.character(strand(gene_exons_all)[1])
  gene_exons  <- gene_exons_all[seqnames(gene_exons_all) == gene_chr & strand(gene_exons_all) == gene_strand]
  gene_start  <- min(start(gene_exons))
  gene_end    <- max(end(gene_exons))
  
  # Filter UTRs to same chromosome/strand and expand gene boundaries if needed
  if (length(gene_utrs_all)) {
    gene_utrs <- gene_utrs_all[seqnames(gene_utrs_all) == gene_chr & strand(gene_utrs_all) == gene_strand]
    if (length(gene_utrs)) {
      utr_start <- min(start(gene_utrs))
      utr_end   <- max(end(gene_utrs))
      gene_start <- min(gene_start, utr_start)
      gene_end   <- max(gene_end, utr_end)
    }
  } else {
    gene_utrs <- GRanges()
  }
  
  gene_gr     <- GRanges(gene_chr, IRanges(gene_start, gene_end), strand = gene_strand)rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(data.table)
  library(grid)   # for unit() used in arrow()
})

# ---------- helpers ----------
.normalize_chr <- function(x) {
  # Normalize chromosome names by removing 'chr' prefix for consistent comparison
  sub("^chr", "", x)
}

.pick_gene_col <- function(gr) {
  cand <- c("gene_id","gene","gene_name","geneid","geneId","GeneID")
  hit <- cand[cand %in% colnames(mcols(gr))]
  if (length(hit)) hit[1] else stop("No gene-like column (gene_id/gene_name) found in GTF metadata.")
}

.pick_tx_col <- function(gr) {
  cand <- c("transcript_id","transcript","tx_id","transcriptId","TxID")
  hit <- cand[cand %in% colnames(mcols(gr))]
  if (length(hit)) hit[1] else stop("No transcript-like column (transcript_id) found in GTF metadata.")
}

# Length-weighted binning to reduce coverage points for very long genes
.bin_intervals <- function(dt, region_start, region_end, target_points = 6000L) {
  # dt: data.table(start, end, score) with 1-based inclusive coords
  if (!nrow(dt)) return(dt)
  span <- region_end - region_start + 1L
  if (span <= target_points || target_points <= 0L) return(dt)

  bin_w <- ceiling(span / target_points)
  b_start <- seq(region_start, region_end, by = bin_w)
  b_end   <- pmin(b_start + bin_w - 1L, region_end)
  bins <- data.table(bstart = b_start, bend = b_end)

  iv <- dt[, .(start, end, score)]
  res <- bins[, {
    s <- bstart; e <- bend
    ov <- iv[start <= e & end >= s]
    if (!nrow(ov)) .(start = s, end = e, score = 0)
    else {
      w <- pmin(ov$end, e) - pmax(ov$start, s) + 1L
      .(start = s, end = e, score = sum(ov$score * w) / sum(w))
    }
  }, by = .(bstart, bend)][, .(start, end, score)]

  res[]
}

# ---------- core plotting ----------
plot_gene_with_window <- function(
  gene_id, bw, gtf, peaks_tsv,
  primer_bed = NULL, qc_tsv = NULL, narrowpeak_file = NULL,
  yaxis_mode = c("percent","depth","log10"),
  max_isoforms = Inf, cov_target_points = 6000L
) {
  yaxis_mode <- match.arg(yaxis_mode)

  # --- GTF & gene selection ---
  gtf_gr <- import(gtf)
  gcol <- .pick_gene_col(gtf_gr)
  tcol <- .pick_tx_col(gtf_gr)

  is_exon <- gtf_gr$type %in% "exon"
  is_utr <- gtf_gr$type %in% c("three_prime_utr", "five_prime_utr")
  gene_mask <- is_exon & (as.character(mcols(gtf_gr)[[gcol]]) == gene_id)
  gene_exons_all <- gtf_gr[gene_mask]
  if (!length(gene_exons_all)) stop(sprintf("Gene '%s' not found in GTF (column %s).", gene_id, gcol))
  
  # Also get UTR regions for this gene
  utr_mask <- is_utr & (as.character(mcols(gtf_gr)[[gcol]]) == gene_id)
  gene_utrs_all <- gtf_gr[utr_mask]

  # lock to first exon’s seqname/strand for safety
  gene_chr    <- as.character(seqnames(gene_exons_all)[1])
  gene_strand <- as.character(strand(gene_exons_all)[1])
  gene_exons  <- gene_exons_all[seqnames(gene_exons_all) == gene_chr & strand(gene_exons_all) == gene_strand]
  gene_start  <- min(start(gene_exons))
  gene_end    <- max(end(gene_exons))
  gene_gr     <- GRanges(gene_chr, IRanges(gene_start, gene_end), strand = gene_strand)

  # --- window from peaks.tsv (single read) ---
  peaks <- fread(peaks_tsv)
  if (!"gene" %in% names(peaks)) stop("peaks.tsv must have a 'gene' column.")
  prow <- peaks[gene == gene_id]
  
  # Handle missing genes gracefully
  if (!nrow(prow)) {
    warning(sprintf("Gene '%s' not present in %s - likely filtered out by peak quality thresholds. Creating gene-only plot.", gene_id, basename(peaks_tsv)))
    # Use entire gene span as window when no peaks are available
    window_start <- gene_start
    window_end   <- gene_end
    has_peaks    <- FALSE
  } else {
    window_start <- as.integer(prow$start[1])
    window_end   <- as.integer(prow$end[1])
    if (is.na(window_start) || is.na(window_end)) stop("peaks.tsv must have integer 'start' and 'end' columns.")
    has_peaks    <- TRUE
  }

  # --- coverage over entire gene (interval ribbons; optional binning) ---
  cov_gr <- import(bw, which = gene_gr)

  # TRUE (unbinned) max from BigWig over the gene span
  max_depth_true <- if (length(cov_gr)) {
    x <- as.numeric(mcols(cov_gr)$score)
    x <- x[is.finite(x)]
    if (length(x)) max(x) else 1
  } else 1

  # Interval table clipped to gene span
  cov_df <- if (length(cov_gr)) {
    data.table(
      start = as.integer(pmax(start(cov_gr), gene_start)),
      end   = as.integer(pmin(end(cov_gr),   gene_end)),
      score = as.numeric(mcols(cov_gr)$score)
    )[start <= end]
  } else data.table(start = integer(0), end = integer(0), score = numeric(0))

  # Optional binning for speed; then scale percent by TRUE max
  if (nrow(cov_df)) {
    cov_df <- .bin_intervals(cov_df, region_start = gene_start, region_end = gene_end,
                             target_points = as.integer(cov_target_points))
    cov_df[, score_pct := 100 * score / max_depth_true]
    cov_df[, score_pct := pmin(score_pct, 100)]  # clamp
  }

  # --- build exon/intron tracks for ALL isoforms (full gene view) ---
  transcripts_all <- unique(as.character(mcols(gene_exons)[[tcol]]))
  if (is.finite(max_isoforms) && length(transcripts_all) > max_isoforms) {
    transcripts_all <- head(transcripts_all, max_isoforms)
  }

  exon_df <- data.table()
  intron_df <- data.table()
  utr_df <- data.table()

  for (tx in transcripts_all) {
    tx_exons <- gene_exons[as.character(mcols(gene_exons)[[tcol]]) == tx]
    if (!length(tx_exons)) next
    tx_exons <- tx_exons[order(start(tx_exons))]
    exon_df <- rbind(exon_df, data.table(
      transcript = tx,
      xmin = as.integer(start(tx_exons)),
      xmax = as.integer(end(tx_exons))
    ))
    
    # Add UTRs for this transcript
    if (length(gene_utrs)) {
      tx_utrs <- gene_utrs[as.character(mcols(gene_utrs)[[tcol]]) == tx]
      if (length(tx_utrs)) {
        utr_df <- rbind(utr_df, data.table(
          transcript = tx,
          xmin = as.integer(start(tx_utrs)),
          xmax = as.integer(end(tx_utrs)),
          utr_type = as.character(tx_utrs$type)
        ))
      }
    }
    
    if (length(tx_exons) > 1) {
      intron_starts <- end(tx_exons)[-length(tx_exons)] + 1L
      intron_ends   <- start(tx_exons)[-1] - 1L
      intron_df <- rbind(intron_df, data.table(
        transcript = tx,
        xstart = as.integer(intron_starts),
        xend   = as.integer(intron_ends)
      ))
    }
  }

  # Stable isoform order by genomic start, then end
  if (nrow(exon_df)) {
    tx_order <- unique(exon_df[order(xmin, xmax), transcript])
  } else {
    tx_order <- transcripts_all
  }
  exon_df[,  y := match(transcript, tx_order)]
  intron_df[, y := match(transcript, tx_order)]
  utr_df[, y := match(transcript, tx_order)]
  valid_transcripts <- tx_order
  primer_lane <- length(valid_transcripts) + 1
  peaks_lane <- primer_lane + 1

  cat("Valid transcripts:", length(valid_transcripts), "\n")
  cat("Primer lane:", primer_lane, "\n") 
  cat("Peaks lane:", peaks_lane, "\n")

  # --- QC subtitle (optional) + annotate what 100% means ---
  qc_lab <- NULL
  if (!is.null(qc_tsv) && file.exists(qc_tsv)) {
    qc <- fread(qc_tsv)
    if ("gene" %in% names(qc)) {
      qrow <- qc[gene == gene_id]
      if (nrow(qrow)) {
        # Handle both old (coverage-based) and new (MACS2-based) QC formats
        if ("window_mean" %in% names(qrow)) {
          # Old format: coverage-based metrics
          qc_lab <- sprintf("window_mean=%.1f; longest_zero_run=%d; strategy=%s",
                            qrow$window_mean[1], qrow$longest_zero_run[1], qrow$strategy[1])
        } else if ("peak_score" %in% names(qrow)) {
          # New format: MACS2-based metrics
          qc_lab <- sprintf("peak_score=%.1f; p-value=%.2f; q-value=%.2f; exonic_frac=%.2f; strategy=%s",
                            qrow$peak_score[1], qrow$peak_pvalue[1], qrow$peak_qvalue[1], 
                            qrow$exonic_fraction[1], qrow$strategy[1])
        }
      }
    }
  }
  if (yaxis_mode == "percent") {
    msg <- sprintf("100%% = %.0f (BigWig units)", max_depth_true)
    qc_lab <- if (is.null(qc_lab)) msg else paste(qc_lab, "|", msg)
  }

  # --- primers (optional) ---
  primer_df <- NULL
  if (!is.null(primer_bed) && file.exists(primer_bed)) {
    pb <- fread(primer_bed, header = FALSE)
    setnames(pb, 1:3, c("chr","start","end"))
    if (ncol(pb) >= 4) setnames(pb, 4, "name")  else pb[, name := "primer"]
    if (ncol(pb) >= 6) setnames(pb, 6, "strand") else pb[, strand := gene_strand]
    # Normalize chromosome names for consistent comparison
    pb <- pb[.normalize_chr(pb$chr) == .normalize_chr(gene_chr) & pb$start <= gene_end & pb$end >= gene_start]
    if (nrow(pb)) {
      # Convert from 0-based BED coordinates to 1-based before clipping
      pb$start <- as.integer(pb$start) + 1L
      pb$end <- as.integer(pb$end)
      pb$start <- pmax(pb$start, gene_start)
      pb$end   <- pmin(pb$end,   gene_end)
      pb$xa    <- ifelse(pb$strand == "-", pb$end,   pb$start)  # 5′ →
      pb$xa_e  <- ifelse(pb$strand == "-", pb$start, pb$end)    # → 3′
      primer_df <- pb
    }
  }

    # --- all MACS2 peaks (optional) ---
  all_peaks_df <- NULL
  if (!is.null(narrowpeak_file) && file.exists(narrowpeak_file)) {
    # Read MACS2 narrowPeak file
    # Columns: chr, start, end, name, score, strand, fold_change, pvalue, qvalue, summit_offset
    np <- fread(narrowpeak_file, header = FALSE,
                col.names = c("chr", "start", "end", "name", "score", "strand", 
                              "fold_change", "pvalue", "qvalue", "summit_offset"))
    
    # Remove any track lines or invalid rows (alternative approach for older data.table)
    if (nrow(np) > 0) {
      np <- np[!grepl("^track", chr) & !is.na(chr) & chr != "" & !grepl("^#", chr)]
      # Ensure numeric columns are properly converted
      np$start <- as.integer(np$start)
      np$end <- as.integer(np$end)
      np$score <- as.numeric(np$score)
      # Remove rows with invalid coordinates
      np <- np[!is.na(start) & !is.na(end) & start >= 0 & end >= 0]
    }
    
    cat("Read", nrow(np), "total peaks from narrowPeak file\n")
    
    # Convert from 0-based to 1-based coordinates BEFORE filtering
    if (nrow(np) > 0) {
      np$start <- as.integer(np$start) + 1L
      np$end <- as.integer(np$end)
    }
    
    cat("Gene chromosome:", gene_chr, "\n")
    cat("Gene region:", gene_start, "-", gene_end, "\n")
    if (nrow(np) > 0) {
      cat("Peak chromosomes in file:", unique(np$chr), "\n")
      cat("Normalized gene chr:", .normalize_chr(gene_chr), "\n")
      cat("Normalized peak chrs:", unique(.normalize_chr(np$chr)), "\n")
    }
    
    # Filter peaks overlapping with current gene (using normalized chromosome names)
    np <- np[.normalize_chr(np$chr) == .normalize_chr(gene_chr) & np$start <= gene_end & np$end >= gene_start]
    cat("Found", nrow(np), "peaks overlapping with gene", gene_id, "\n")
    if (nrow(np)) {
      # Clip to gene boundaries after coordinate conversion
      np$start <- pmax(np$start, gene_start)
      np$end   <- pmin(np$end, gene_end)
      
      # Calculate peak intensity (use score, but could also use fold_change or -log10(pvalue))
      np$peak_intensity <- np$score
      
      all_peaks_df <- np[np$start <= np$end]  # Remove any invalid ranges after clipping
      
      cat("Found", nrow(all_peaks_df), "peaks overlapping with gene", gene_id, "\n")
      if (nrow(all_peaks_df) > 0) {
        cat("Peak coordinates after clipping:\n")
        for (i in 1:min(5, nrow(all_peaks_df))) {
          cat("  Peak", i, ": chr", all_peaks_df$chr[i], ":", all_peaks_df$start[i], "-", all_peaks_df$end[i], "\n")
        }
      }
    } else {
      all_peaks_df <- NULL
    }
  }

  # ---------- plots ----------
  # Update title to indicate if peaks are missing
  if (has_peaks) {
    cov_title <- sprintf("%s  %s:%d-%d (%s)", gene_id, gene_chr, window_start, window_end, gene_strand)
  } else {
    cov_title <- sprintf("%s  %s:%d-%d (%s) [NO PEAKS - GENE VIEW ONLY]", gene_id, gene_chr, window_start, window_end, gene_strand)
  }
  win_highlight <- data.frame(xmin = window_start, xmax = window_end,
                              ymin = 0.5, ymax = peaks_lane + 0.5)

  # Coverage plot (mode-dependent; percent uses TRUE max)
  p_cov <- ggplot()
  if (nrow(cov_df)) {
    if (yaxis_mode == "percent") {
      # Use bin midpoints for area plot
      cov_df_area <- cov_df[, .(pos = round((start + end)/2), score_pct = score_pct)]
      p_cov <- p_cov +
        geom_area(data = cov_df_area, aes(x = pos, y = score_pct), fill = "steelblue", alpha = 0.5) +
        scale_y_continuous(labels = function(x) sprintf("%d%%", round(x)),
                           limits = c(0, 100), name = "Coverage (% of peak)")
    } else if (yaxis_mode == "depth") {
      cov_df_area <- cov_df[, .(pos = round((start + end)/2), score = score)]
      p_cov <- p_cov +
        geom_area(data = cov_df_area, aes(x = pos, y = score), fill = "steelblue", alpha = 0.5) +
        labs(y = "Coverage (depth)")
    } else if (yaxis_mode == "log10") {
      cov_df_area <- cov_df[, .(pos = round((start + end)/2), score = score)]
      p_cov <- p_cov +
        geom_area(data = cov_df_area, aes(x = pos, y = score), fill = "steelblue", alpha = 0.5) +
        scale_y_continuous(trans = "log10", name = "Coverage (log10 depth)")
    }
  } else {
    # empty coverage; still set axis titles
    if (yaxis_mode == "percent") {
      p_cov <- p_cov + scale_y_continuous(labels = function(x) sprintf("%d%%", round(x)),
                                           limits = c(0, 100), name = "Coverage (% of peak)")
    } else if (yaxis_mode == "depth") {
      p_cov <- p_cov + labs(y = "Coverage (depth)")
    } else {
      p_cov <- p_cov + scale_y_continuous(trans = "log10", name = "Coverage (log10 depth)")
    }
  }
  p_cov <- p_cov +
    scale_x_continuous(limits = c(gene_start, gene_end)) +
    labs(x = NULL, title = cov_title, subtitle = qc_lab) +
    theme_bw(base_size = 11) +
    theme(plot.title.position = "plot",
          plot.subtitle = element_text(size = 9),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  # Add peak highlight only if peaks are available
  if (has_peaks) {
    p_cov <- p_cov + annotate("rect", xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf, alpha = 0.12)
  }

  # Feature track: all isoforms + primer lane
  feat <- ggplot()
  
  # Add window outline only if peaks are available
  if (has_peaks) {
    feat <- feat + geom_rect(data = win_highlight,
                             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                             fill = "red", alpha = 0.18, inherit.aes = FALSE)
  }
  
  feat <- feat +
    # exons (standard height)
    { if (nrow(exon_df))
        geom_rect(data = exon_df,
                  aes(xmin = xmin, xmax = xmax, ymin = y, ymax = y + 0.8),
                  fill = "steelblue", alpha = 0.35) else NULL } +
    # UTRs (narrower height - centered on transcript line)
    { if (nrow(utr_df)) {
        # 5' UTRs in lighter blue, 3' UTRs in orange
        utr_5prime <- utr_df[utr_type == "five_prime_utr"]
        utr_3prime <- utr_df[utr_type == "three_prime_utr"]
        
        list(
          if (nrow(utr_5prime))
            geom_rect(data = utr_5prime,
                      aes(xmin = xmin, xmax = xmax, ymin = y + 0.25, ymax = y + 0.55),
                      fill = "lightblue", alpha = 0.7, color = "darkblue", linewidth = 0.3),
          if (nrow(utr_3prime))
            geom_rect(data = utr_3prime,
                      aes(xmin = xmin, xmax = xmax, ymin = y + 0.25, ymax = y + 0.55),
                      fill = "orange", alpha = 0.7, color = "darkorange", linewidth = 0.3)
        )
      } else NULL } +
    # introns
    { if (nrow(intron_df))
        geom_segment(data = intron_df,
                     aes(x = xstart, xend = xend, y = y + 0.4, yend = y + 0.4),
                     linewidth = 0.7) else NULL } +
    # primer arrows on their own lane
    { if (!is.null(primer_df) && nrow(primer_df))
        geom_segment(data = primer_df,
                     aes(x = xa, xend = xa_e, y = primer_lane + 0.2, yend = primer_lane + 0.2),
                     arrow = arrow(length = unit(3, "mm"), ends = "last"),
                     linewidth = 0.9) else NULL } +
    # all MACS2 peaks as horizontal bars  
    { if (!is.null(all_peaks_df) && nrow(all_peaks_df)) {
        cat("Drawing", nrow(all_peaks_df), "peaks at y-level", peaks_lane, "\n")
        cat("Y coordinates: ymin =", peaks_lane + 0.1, "ymax =", peaks_lane + 0.7, "\n")
        cat("X-axis limits: gene_start =", gene_start, "gene_end =", gene_end, "\n")
        for (i in 1:min(3, nrow(all_peaks_df))) {
          cat("  Peak", i, "x-coords: xmin =", all_peaks_df$start[i], "xmax =", all_peaks_df$end[i], "\n")
        }
        geom_rect(data = all_peaks_df,
                  aes(xmin = start, xmax = end, 
                      ymin = peaks_lane + 0.1, ymax = peaks_lane + 0.7),
                  fill = "darkred", alpha = 0.7, color = "black", linewidth = 0.2)
      } else NULL } +
    scale_x_continuous(limits = c(gene_start, gene_end)) +
    scale_y_continuous(
      limits = c(0.5, peaks_lane + 1.0),  # Increase upper limit to accommodate peak rectangles
      breaks = c(seq_along(valid_transcripts), primer_lane, peaks_lane),
      labels = c(valid_transcripts, "Primer", "All Peaks")
    ) +
    labs(x = sprintf("%s:%d-%d", gene_chr, gene_start, gene_end), y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7),
          legend.position = "none")

  # Stack tracks
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p_cov, feat, ncol = 1, heights = c(3, 2))  # Increase feature track height
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(p_cov, feat, ncol = 1, rel_heights = c(3, 2))  # Increase feature track height
  } else {
    # fallback if patchwork/cowplot not installed
    list(coverage = p_cov, features = feat)
  }
}

# ---------- CLI ----------
opt_list <- list(
  make_option("--gene", type="character", help="Gene identifier (must match 'gene' in peaks.tsv and the GTF column)"),
  make_option("--bw", type="character", help="Coverage BigWig file"),
  make_option("--gtf", type="character", help="Gene annotation GTF/GFF"),
  make_option("--peaks", type="character", default="peaks.tsv", help="peaks.tsv from your pipeline [default %default]"),
  make_option("--qc", type="character", default="qc_coverage_summary.tsv", help="QC table (optional) [default %default]"),
  make_option("--primer", type="character", default=NULL, help="Primer BED (optional)"),
  make_option("--narrowpeak", type="character", default=NULL, help="MACS2 narrowPeak file to show all peaks (optional)"),
  make_option("--out", type="character", default=NULL, help="Output image path (PNG/PDF). If omitted, uses plot_<gene>.png"),
  make_option("--yaxis", type="character", default="percent",
              help="Y-axis mode: 'percent' (relative to TRUE BigWig peak), 'depth' (raw), or 'log10' (raw log10) [default %default]"),
  make_option("--max_isoforms", type="integer", default=Inf, help="Limit number of isoforms for readability (default: all)"),
  make_option("--cov_points", type="integer", default=6000, help="Target #points for coverage binning across gene [default %default]"),
  make_option("--width", type="double", default=12, help="Figure width inches [default %default]"),
  make_option("--height", type="double", default=6.5, help="Figure height inches [default %default]"),
  make_option("--dpi", type="double", default=300, help="DPI for raster outputs [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$gene) || is.null(opt$bw) || is.null(opt$gtf))
  stop("Required: --gene --bw --gtf (and ensure --peaks points to your peaks.tsv)")

plt <- plot_gene_with_window(
  gene_id          = opt$gene,
  bw               = opt$bw,
  gtf              = opt$gtf,
  peaks_tsv        = opt$peaks,
  primer_bed       = opt$primer,
  qc_tsv           = if (file.exists(opt$qc)) opt$qc else NULL,
  narrowpeak_file  = opt$narrowpeak,
  yaxis_mode       = opt$yaxis,
  max_isoforms     = opt$max_isoforms,
  cov_target_points= as.integer(opt$cov_points)
)

outfile <- if (!is.null(opt$out)) opt$out else sprintf("plot_%s.png", opt$gene)
ext <- tolower(tools::file_ext(outfile))

if (inherits(plt, "list")) {
  # fallback: use gridExtra to combine panels when patchwork/cowplot unavailable
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plt <- gridExtra::grid.arrange(plt$coverage, plt$features, ncol = 1, heights = c(3, 2))
    ggsave(outfile, combined_plt, width = opt$width, height = opt$height, dpi = opt$dpi)
  } else {
    warning("Neither patchwork, cowplot, nor gridExtra available. Saving coverage panel only.")
    ggsave(outfile, plt$coverage, width = opt$width, height = opt$height, dpi = opt$dpi)
  }
} else {
  if (ext %in% c("pdf","svg")) {
    ggsave(outfile, plt, width = opt$width, height = opt$height, limitsize = FALSE)
  } else {
    ggsave(outfile, plt, width = opt$width, height = opt$height, dpi = opt$dpi)
  }
}
cat("Wrote", outfile, "\n")
