#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(derfinder)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
  library(Biostrings)
  library(Rsamtools)
})

opt_list <- list(
  make_option("--bw", type="character"),
  make_option("--gtf", type="character"),
  make_option("--genes", type="character"),
  make_option("--genome", type="character", default="BSgenome.Hsapiens.UCSC.hg38"),
  make_option("--fasta", type="character", default=NULL),
  make_option("--pad", type="integer", default=60),
  make_option("--smooth_k", type="integer", default=31),
  make_option("--sliding_window", action="store_true", default=FALSE, help="Enable sliding-window max-average selection if default window fails QC"),
  make_option("--min_window_mean", type="double", default=NA_real_, help="Minimum mean coverage across the final window"),
  make_option("--min_window_mean_pct", type="double", default=NA_real_, help="Minimum window mean as percentage of gene's peak coverage (0-100, overrides --min_window_mean)"),
  make_option("--max_gap", type="integer", default=NA_integer_, help="Maximum allowed longest zero-coverage run within the final window"),
  make_option("--search_slop", type="integer", default=1000, help="Extra bases around exon span for coverage import"),
  make_option("--trim_to_exon", action="store_true", default=FALSE, help="Trim final window to boundaries of exon containing the peak (prevents intronic spill-over)"),
  make_option("--trim_low_coverage_pct", type="double", default=NA_real_, help="After window selection, trim ends with coverage below X% of window peak (0-100, e.g. 10 for 10%)"),
  make_option("--min_exonic_fraction", type="double", default=NA_real_, help="If set and window exonic fraction < value, mark as failing (reported in QC); no automatic rescue unless --trim_to_exon used"),
  make_option("--out_fa", type="character", default="primer_targets.fa"),
  make_option("--out_bed", type="character", default="primer_targets.bed"),
  make_option("--out_peaks", type="character", default="peaks.tsv"),
  make_option("--out_qc", type="character", default="qc_coverage_summary.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Load BSgenome dynamically (only if needed)
if (is.null(opt$fasta)) {
  suppressWarnings(suppressMessages({
    genome_pkg <- opt$genome
    library(genome_pkg, character.only=TRUE)
  }))
}

# Helpers
runmean_rle <- function(x, k=31L) {
  if (k < 1) return(x)
  if (k %% 2 == 0) k <- k + 1L
  as.numeric(S4Vectors::runmean(Rle(x), k=k, endrule="constant"))
}

# Inputs
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
g_all <- genes(txdb)
ids <- unique(scan(opt$genes, what=character(), quiet=TRUE))
sel <- names(g_all) %in% ids
if (!any(sel)) stop("None of the provided gene IDs were found in the GTF.")
g <- g_all[sel]

# Union of exons per gene (isoform‑agnostic exonic parts)
ex_by_gene <- exonsBy(txdb, by="gene")
ex_by_gene <- ex_by_gene[names(ex_by_gene) %in% names(g)]
ex_merged  <- endoapply(ex_by_gene, function(x) reduce(trim(x)))

# Union all exon ranges for coverage query
ex_all <- unlist(ex_merged, use.names=TRUE)
if (length(ex_all) == 0) stop("No exons for selected genes.")

# Debug: print chromosome names
cat("GTF chromosome names:", head(unique(as.character(seqnames(ex_all)))), "\n")

# Filter to only the chromosomes where our genes are located
target_chroms <- unique(as.character(seqnames(ex_all)))
cat("Target chromosomes:", target_chroms, "\n")

# Filter to standard chromosomes only (exclude scaffolds/patches)
standard_chroms <- c(as.character(1:22), "X", "Y", "MT")
target_chroms_filtered <- intersect(target_chroms, standard_chroms)

cat("Standard chromosomes needed:", target_chroms_filtered, "\n")

if (length(target_chroms_filtered) == 0) {
  stop("No genes found on standard chromosomes (1-22, X, Y, MT). Your genes may be on alternative scaffolds not supported by the BigWig.")
}

# Create new GRanges with only standard chromosomes
ex_std <- ex_all[seqnames(ex_all) %in% target_chroms_filtered]

# Create a completely new GRanges object to avoid seqlevel issues
ex_clean <- GRanges(
  seqnames = Rle(as.character(seqnames(ex_std))),
  ranges = ranges(ex_std),
  strand = strand(ex_std)
)
names(ex_clean) <- names(ex_std)

cat("Clean regions created:", length(ex_clean), "regions on", length(unique(seqnames(ex_clean))), "chromosomes\n")

# Debug: Check chromosome names in detail
cat("Clean GRanges seqlevels:", paste(seqlevels(ex_clean), collapse=", "), "\n")
cat("Clean GRanges seqnames:", paste(unique(as.character(seqnames(ex_clean))), collapse=", "), "\n")

# Try alternative approach using rtracklayer directly instead of derfinder
tryCatch({
  library(rtracklayer)
  
  cat("Using simple BigWig reading approach...\n")
  
  # Read coverage data chromosome by chromosome with a simpler approach
  reg_cov <- list()
  cov_cache <- list()
  
  for (chr_name in target_chroms_filtered) {
    cat("Reading BigWig data for chromosome", chr_name, "\n")
    
    # Get all regions for this chromosome
    chr_regions <- ex_clean[seqnames(ex_clean) == chr_name]
    
    if (length(chr_regions) > 0) {
      # Create a range spanning all regions on this chromosome (with some padding)
      chr_start <- min(start(chr_regions)) - 1000
      chr_end <- max(end(chr_regions)) + 1000
      chr_start <- max(1, chr_start)  # Don't go below 1
      
      cat("  Querying chr", chr_name, "from", chr_start, "to", chr_end, "\n")
      
      # Import coverage for this region
      query_range <- GRanges(seqnames = chr_name, 
                           ranges = IRanges(start = chr_start, end = chr_end))
      
      chr_coverage <- import(opt$bw, which = query_range)
      cat("  Found", length(chr_coverage), "coverage intervals in BigWig\n")
      
      # Debug: show some example coverage intervals
      if (length(chr_coverage) > 0) {
        cat("  Example coverage intervals:\n")
        for (k in 1:min(3, length(chr_coverage))) {
          cat("    ", start(chr_coverage[k]), "-", end(chr_coverage[k]), " score:", chr_coverage[k]$score, "\n")
        }
      } else {
        cat("  WARNING: No coverage intervals found in BigWig for this region!\n")
        # Try querying a larger region to see if there's any data
        larger_range <- GRanges(seqnames = chr_name, 
                               ranges = IRanges(start = 1, end = 1000000))
        test_cov <- import(opt$bw, which = larger_range)
        cat("  Testing larger region (1-1M): found", length(test_cov), "intervals\n")
      }
      
  # Convert to Rle for easier processing
      if (length(chr_coverage) > 0) {
        # Create an RLE manually from the BigWig intervals
        # Initialize with zeros
        chr_rle <- Rle(0, chr_end - chr_start + 1)
        
        # Fill in the coverage values from BigWig intervals
        for (j in seq_along(chr_coverage)) {
          interval <- chr_coverage[j]
          # Convert genomic coordinates to RLE coordinates
          rle_start <- start(interval) - chr_start + 1
          rle_end <- end(interval) - chr_start + 1
          
          # Ensure we're within bounds
          rle_start <- max(1, rle_start)
          rle_end <- min(length(chr_rle), rle_end)
          
          if (rle_start <= rle_end && rle_start <= length(chr_rle)) {
            # Set the coverage values
            chr_rle[rle_start:rle_end] <- interval$score
          }
        }
      } else {
        # Create empty coverage if no data
        chr_rle <- Rle(0, chr_end - chr_start + 1)
      }
      
  cat("  Coverage RLE length:", length(chr_rle), "values\n")
      cat("  Coverage RLE summary: min =", min(chr_rle), "max =", max(chr_rle), "mean =", mean(chr_rle), "\n")
  # Save cache for this chromosome
  cov_cache[[chr_name]] <- list(rle = chr_rle, offset = chr_start)
      
      # Extract coverage for each region on this chromosome
      for (i in which(seqnames(ex_clean) == chr_name)) {
        region <- ex_clean[i]
        region_start <- start(region) - chr_start + 1
        region_end <- end(region) - chr_start + 1
        
        cat("  Region", i, ":", as.character(seqnames(region)), ":", start(region), "-", end(region), 
            "(RLE coords:", region_start, "-", region_end, ")\n")
        
        # Ensure we don't go out of bounds
        if (region_start > 0 && region_start <= length(chr_rle)) {
          region_end <- min(region_end, length(chr_rle))
          region_start <- max(1, region_start)
          
          if (region_start <= region_end) {
            cov_values <- as.numeric(chr_rle[region_start:region_end])
            cat("    Extracted", length(cov_values), "coverage values, max =", max(cov_values), "\n")
          } else {
            cov_values <- numeric(width(region))
            cat("    Invalid region coordinates, using zeros\n")
          }
        } else {
          cov_values <- numeric(width(region))
          cat("    Region out of bounds, using zeros\n")
        }
        
        # Ensure the coverage vector has the right length
        if (length(cov_values) < width(region)) {
          cov_values <- c(cov_values, rep(0, width(region) - length(cov_values)))
        } else if (length(cov_values) > width(region)) {
          cov_values <- cov_values[1:width(region)]
        }
        
        reg_cov[[i]] <- list(value = cov_values)
      }
      
      cat("Successfully processed", sum(seqnames(ex_clean) == chr_name), "regions on chromosome", chr_name, "\n")
    }
  }
  
  cat("Successfully got coverage for", length(reg_cov), "regions using simple approach\n")
  
}, error = function(e) {
  cat("Error with rtracklayer approach:", e$message, "\n")
  
  # Fallback: try the original derfinder approach with more specific debugging
  tryCatch({
    cat("Falling back to derfinder approach...\n")
    cat("Attempting to load coverage for chromosomes:", paste(target_chroms_filtered, collapse=", "), "\n")
    
    # Get coverage chromosome by chromosome to isolate any issues
    reg_cov <- list()
    for (chr_name in target_chroms_filtered) {
      cat("Processing chromosome", chr_name, "\n")
      chr_regions <- ex_clean[seqnames(ex_clean) == chr_name]
      if (length(chr_regions) > 0) {
        chr_cov <- getRegionCoverage(regions = chr_regions, files = opt$bw, verbose = FALSE)
        reg_cov <- c(reg_cov, chr_cov)
      }
    }
    
    cat("Successfully got coverage for", length(reg_cov), "regions\n")
  }, error = function(e2) {
    cat("Error getting coverage:", e2$message, "\n")
    
    # Try alternative: use rtracklayer to check BigWig contents
    library(rtracklayer)
    tryCatch({
      cat("Checking BigWig file contents...\n")
      bw_info <- seqinfo(BigWigFile(opt$bw))
      cat("BigWig seqnames:", paste(seqnames(bw_info), collapse=", "), "\n")
      cat("BigWig seqlengths for target chromosomes:\n")
      for (chr in target_chroms_filtered) {
        if (chr %in% seqnames(bw_info)) {
          cat("  ", chr, ":", seqlengths(bw_info)[chr], "bp\n")
        } else {
          cat("  ", chr, ": NOT FOUND\n")
        }
      }
    }, error = function(e3) {
      cat("Could not check BigWig contents:", e3$message, "\n")
    })
    
    stop("Failed to get coverage data. Check that the BigWig file contains the required chromosomes: ", paste(target_chroms_filtered, collapse=", "))
  })
})

# Create region2gene mapping for the clean regions
region2gene <- character()
clean_names <- names(ex_clean)
for (name in clean_names) {
  gene_id <- strsplit(name, "\\.")[[1]][1]
  region2gene <- c(region2gene, gene_id)
}

# Iterate genes
peak_rows <- list()
qc_rows <- list()
genes_processed <- 0
genes_with_coverage <- 0

cat("Processing", length(names(g)), "genes for peak detection\n")

get_window_cov <- function(chr, start_pos, end_pos) {
  if (!(chr %in% names(cov_cache))) return(rep(NA_real_, end_pos - start_pos + 1))
  cc <- cov_cache[[chr]]
  start_idx <- start_pos - cc$offset + 1
  end_idx <- end_pos - cc$offset + 1
  if (is.na(start_idx) || is.na(end_idx) || end_idx < 1 || start_idx > length(cc$rle)) {
    return(rep(NA_real_, end_pos - start_pos + 1))
  }
  start_idx <- max(1, start_idx)
  end_idx <- min(length(cc$rle), end_idx)
  vec <- as.numeric(cc$rle[start_idx:end_idx])
  left_pad <- (start_pos - cc$offset + 1) - start_idx
  right_pad <- (end_pos - cc$offset + 1) - end_idx
  if (left_pad > 0) vec <- c(rep(NA_real_, left_pad), vec)
  if (right_pad > 0) vec <- c(vec, rep(NA_real_, right_pad))
  need <- end_pos - start_pos + 1
  if (length(vec) < need) vec <- c(vec, rep(NA_real_, need - length(vec)))
  if (length(vec) > need) vec <- vec[seq_len(need)]
  vec
}

longest_zero_run <- function(x) {
  if (!length(x)) return(0L)
  y <- ifelse(is.na(x), 0, x)
  r <- rle(y == 0)
  if (!length(r$lengths)) return(0L)
  zero_runs <- r$lengths[r$values]
  if (length(zero_runs) == 0) return(0L)
  max(zero_runs)
}

# Trim window ends with coverage below threshold percentage of window peak
trim_low_coverage <- function(coverage_vec, threshold_pct) {
  if (is.na(threshold_pct) || threshold_pct <= 0) {
    return(list(start_trim = 0, end_trim = 0, trimmed_length = length(coverage_vec)))
  }
  
  cov_nona <- replace(coverage_vec, is.na(coverage_vec), 0)
  window_peak <- max(cov_nona)
  threshold <- window_peak * (threshold_pct / 100)
  
  if (window_peak == 0) {
    return(list(start_trim = 0, end_trim = 0, trimmed_length = length(coverage_vec)))
  }
  
  # Find first position >= threshold from start
  start_trim <- 0
  for (i in seq_along(cov_nona)) {
    if (cov_nona[i] >= threshold) break
    start_trim <- start_trim + 1
  }
  
  # Find first position >= threshold from end
  end_trim <- 0
  for (i in rev(seq_along(cov_nona))) {
    if (cov_nona[i] >= threshold) break
    end_trim <- end_trim + 1
  }
  
  # Ensure we don't trim everything
  trimmed_length <- length(cov_nona) - start_trim - end_trim
  if (trimmed_length <= 0) {
    # If trimming would remove everything, keep the original window
    start_trim <- 0
    end_trim <- 0
    trimmed_length <- length(coverage_vec)
  }
  
  list(start_trim = start_trim, end_trim = end_trim, trimmed_length = trimmed_length)
}

find_best_window <- function(chr, search_start, search_end, wsize, max_gap = NA_integer_) {
  if (!(chr %in% names(cov_cache))) return(NULL)
  cc <- cov_cache[[chr]]
  s_idx <- max(1, as.integer(search_start - cc$offset + 1))
  e_idx <- min(length(cc$rle), as.integer(search_end - cc$offset + 1))
  if (e_idx - s_idx + 1 < wsize) return(NULL)
  vec <- as.numeric(cc$rle[s_idx:e_idx])
  n <- length(vec)
  v0 <- replace(vec, is.na(vec), 0)
  cs <- c(0, cumsum(v0))
  best_mean <- -Inf
  best_center <- NA_integer_
  for (i in 0:(n - wsize)) {
    j <- i + wsize
    sumv <- cs[j] - cs[i]
    meanv <- sumv / wsize
    if (!is.na(max_gap)) {
      sub <- vec[(i + 1):j]
      if (longest_zero_run(sub) > max_gap) next
    }
    if (meanv > best_mean) {
      best_mean <- meanv
      center_idx <- floor((i + 1 + j) / 2)
      best_center <- (s_idx + center_idx - 1) + cc$offset - 1
    }
  }
  if (is.infinite(best_mean)) return(NULL)
  list(center = best_center, mean = best_mean)
}

for (gid in names(g)) {
  idx <- which(region2gene == gid)
  if (!length(idx)) {
    cat("Gene", gid, "has no matching regions\n")
    next
  }
  
  genes_processed <- genes_processed + 1

  blocks <- ex_merged[[gid]]
  cov_vec <- unlist(lapply(reg_cov[idx], function(rr) rr$value), use.names = FALSE)
  
  cat("Gene", gid, "- regions:", length(idx), "coverage length:", length(cov_vec), "max coverage:", max(cov_vec, na.rm = TRUE), "\n")
  
  if (!length(cov_vec) || max(cov_vec) == 0) {
    cat("  Skipping gene", gid, "- no coverage\n")
    next
  }
  
  genes_with_coverage <- genes_with_coverage + 1

  sm <- runmean_rle(cov_vec, k = opt$smooth_k)
  max_idx <- which.max(sm)

  cumw <- cumsum(width(blocks))
  block_i <- which(max_idx <= cumw)[1]
  offset_in_block <- max_idx - ifelse(block_i == 1, 0, cumw[block_i - 1])
  peak_pos <- start(blocks[block_i]) + offset_in_block - 1

  chr_g <- as.character(seqnames(blocks[block_i]))
  strand_g <- as.character(strand(g[gid]))
  pad <- opt$pad
  win_start <- max(1L, as.integer(peak_pos - pad))
  win_end <- as.integer(peak_pos + pad)
  strategy <- "peak_centered"

  # Compute dynamic threshold if percentage-based
  dynamic_min_mean <- opt$min_window_mean
  if (!is.na(opt$min_window_mean_pct)) {
    gene_peak_cov <- max(cov_vec)
    dynamic_min_mean <- gene_peak_cov * (opt$min_window_mean_pct / 100)
    cat("  Gene", gid, "peak coverage:", gene_peak_cov, "-> dynamic threshold:", dynamic_min_mean, "(", opt$min_window_mean_pct, "%)\n")
  }

  wcov <- get_window_cov(chr_g, win_start, win_end)
  wcov_nona <- replace(wcov, is.na(wcov), 0)
  w_mean <- mean(wcov_nona)
  w_median <- stats::median(wcov_nona)
  w_min <- min(wcov_nona)
  w_max <- max(wcov_nona)
  w_zeros <- sum(wcov_nona == 0)
  w_lzr <- longest_zero_run(wcov_nona)
  pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
  pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)

  # Debug: print sliding window decision details
  cat("  Gene", gid, "initial window QC:\n")
  cat("    window_mean =", w_mean, "(threshold:", dynamic_min_mean, ")\n")
  cat("    longest_zero_run =", w_lzr, "(threshold:", opt$max_gap, ")\n")
  cat("    pass_min_mean =", pass_min_mean, ", pass_max_gap =", pass_max_gap, "\n")
  
  if (opt$sliding_window && (!pass_min_mean || !pass_max_gap)) {
    cat("    Triggering sliding window search...\n")
    search_start <- max(1L, min(start(blocks)) - pad)
    search_end <- max(end(blocks)) + pad
    cat("    Search range:", search_start, "to", search_end, "(window size:", 2*pad+1, "bp)\n")
    best <- find_best_window(chr_g, search_start, search_end, wsize = (2 * pad + 1L), max_gap = opt$max_gap)
    if (!is.null(best)) {
      cat("    Found better window at position", best$center, "with mean coverage", best$mean, "\n")
      strategy <- "sliding_best"
      peak_pos <- as.integer(best$center)
      win_start <- max(1L, as.integer(peak_pos - pad))
      win_end <- as.integer(peak_pos + pad)
      wcov <- get_window_cov(chr_g, win_start, win_end)
      wcov_nona <- replace(wcov, is.na(wcov), 0)
      w_mean <- mean(wcov_nona)
      w_median <- stats::median(wcov_nona)
      w_min <- min(wcov_nona)
      w_max <- max(wcov_nona)
      w_zeros <- sum(wcov_nona == 0)
      w_lzr <- longest_zero_run(wcov_nona)
      pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
      pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)
    } else {
      cat("    No better window found, keeping original\n")
    }
  }

  # Apply coverage-based trimming if requested
  original_start <- win_start
  original_end <- win_end
  trim_applied <- FALSE
  if (!is.na(opt$trim_low_coverage_pct)) {
    trim_result <- trim_low_coverage(wcov, opt$trim_low_coverage_pct)
    if (trim_result$start_trim > 0 || trim_result$end_trim > 0) {
      trim_applied <- TRUE
      new_start <- win_start + trim_result$start_trim
      new_end <- win_end - trim_result$end_trim
      
      cat("  Gene", gid, "coverage trimming (", opt$trim_low_coverage_pct, "% threshold):\n")
      cat("    Original window:", win_start, "-", win_end, "(", win_end - win_start + 1, "bp)\n")
      cat("    Trimmed window: ", new_start, "-", new_end, "(", new_end - new_start + 1, "bp)\n")
      cat("    Trimmed:", trim_result$start_trim, "bp from start,", trim_result$end_trim, "bp from end\n")
      
      # Update window coordinates
      win_start <- new_start
      win_end <- new_end
      peak_pos <- as.integer((win_start + win_end) / 2)  # Recenter peak
      strategy <- paste0(strategy, "+trimmed")
      
      # Recalculate coverage statistics for trimmed window
      wcov <- get_window_cov(chr_g, win_start, win_end)
      wcov_nona <- replace(wcov, is.na(wcov), 0)
      w_mean <- mean(wcov_nona)
      w_median <- stats::median(wcov_nona)
      w_min <- min(wcov_nona)
      w_max <- max(wcov_nona)
      w_zeros <- sum(wcov_nona == 0)
      w_lzr <- longest_zero_run(wcov_nona)
      pass_min_mean <- if (is.na(dynamic_min_mean)) TRUE else (w_mean >= dynamic_min_mean)
      pass_max_gap <- if (is.na(opt$max_gap)) TRUE else (w_lzr <= opt$max_gap)
    }
  }

  peak_rows[[gid]] <- data.frame(
    gene = gid,
    chr = chr_g,
    pos = peak_pos,
    start = win_start,
    end = win_end,
    strand = strand_g,
    stringsAsFactors = FALSE
  )

  qc_rows[[gid]] <- data.frame(
    gene = gid,
    total_exonic_bases = sum(width(blocks)),
    max_cov = max(cov_vec),
    mean_cov = mean(cov_vec),
    median_cov = stats::median(cov_vec),
    peak_pos = peak_pos,
    window_start = win_start,
    window_end = win_end,
    window_mean = w_mean,
    window_median = w_median,
    window_min = w_min,
    window_max = w_max,
    window_zeros = w_zeros,
    longest_zero_run = w_lzr,
    pass_min_mean = pass_min_mean,
    pass_max_gap = pass_max_gap,
    strategy = strategy,
    stringsAsFactors = FALSE
  )
}

cat("Summary: Processed", genes_processed, "genes,", genes_with_coverage, "had coverage\n")

peaks_df <- do.call(rbind, peak_rows)
if (is.null(peaks_df) || nrow(peaks_df) == 0) stop("No peaks found — check coverage and gene IDs.")
qc_df <- do.call(rbind, qc_rows)

# Build windows (start/end already computed with pad)
gr_peaks <- GRanges(peaks_df$chr,
                    IRanges(pmax(1, as.integer(peaks_df$start)), as.integer(peaks_df$end)),
                    strand = peaks_df$strand,
                    gene = peaks_df$gene)

# Annotate exon overlap / trimming if requested
exonic_fraction <- rep(NA_real_, length(gr_peaks))
trimmed <- rep(FALSE, length(gr_peaks))
fail_exonic_fraction <- rep(FALSE, length(gr_peaks))

for (i in seq_along(gr_peaks)) {
  gid <- peaks_df$gene[i]
  blocks <- ex_merged[[gid]]
  w <- gr_peaks[i]
  # Use findOverlaps + intersect instead of pintersect to avoid length issues
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
    # Find exon that contains peak_pos (peak stored as pos)
    peak_pos <- peaks_df$pos[i]
    containing <- blocks[peak_pos >= start(blocks) & peak_pos <= end(blocks)]
    if (length(containing) >= 1) {
      # Trim the window to the exon boundaries, not use the entire exon
      exon_range <- containing[1]
      current_start <- start(gr_peaks[i])
      current_end <- end(gr_peaks[i])
      new_start <- max(current_start, start(exon_range))
      new_end <- min(current_end, end(exon_range))
      gr_peaks[i] <- GRanges(seqnames(exon_range), IRanges(new_start, new_end), strand = strand(exon_range), gene = gid)
      trimmed[i] <- TRUE
      # Recompute exonic fraction after trimming
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

# Update QC with exonic fraction info (align rows by gene)
qc_df$exonic_fraction <- exonic_fraction[match(qc_df$gene, peaks_df$gene)]
qc_df$fail_exonic_fraction <- fail_exonic_fraction[match(qc_df$gene, peaks_df$gene)]
qc_df$trimmed_to_exon <- trimmed[match(qc_df$gene, peaks_df$gene)]

# Extract sequences
# Choose sequence source: FASTA (preferred if provided) or BSgenome
if (!is.null(opt$fasta)) {
  faf <- Rsamtools::FaFile(opt$fasta)
  Rsamtools::open.FaFile(faf)
  on.exit(Rsamtools::close.FaFile(faf), add = TRUE)
  seqs <- Rsamtools::getSeq(faf, gr_peaks)   # names taken from GRanges
} else {
  # Try to get BSgenome object
  tryCatch({
    bs <- get(opt$genome)                      # BSgenome fallback
    seqs <- Biostrings::getSeq(bs, gr_peaks)
  }, error = function(e) {
    stop("No FASTA file provided and BSgenome package '", opt$genome, "' is not available. Please provide a FASTA file with --fasta or install the correct BSgenome package.")
  })
}

nm <- paste0(peaks_df$gene, "|", as.character(seqnames(gr_peaks)), ":",
             start(gr_peaks), "-", end(gr_peaks), "(", peaks_df$strand, ")")
names(seqs) <- nm

writeXStringSet(seqs, filepath = opt$out_fa)
export(gr_peaks, con = opt$out_bed, format = "BED")

peaks_out <- data.frame(
  gene = peaks_df$gene,
  chr = as.character(seqnames(gr_peaks)),
  start = start(gr_peaks),
  end = end(gr_peaks),
  strand = peaks_df$strand,
  exonic_fraction = peaks_df$exonic_fraction,
  trimmed_to_exon = peaks_df$trimmed_to_exon,
  fail_exonic_fraction = peaks_df$fail_exonic_fraction
)
write.table(peaks_out, file = opt$out_peaks, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(qc_df, file = opt$out_qc, sep = "\t", quote = FALSE, row.names = FALSE)