#!/usr/bin/env Rscript

# PeakPrime Standalone Explorer - Run from anywhere, point to any results directory
# Full plotting features with proper separated layout (coverage + gene structure)
#
# Usage:
#   Rscript app_standalone.R /path/to/results/directory
#   OR
#   Rscript app_standalone.R  # Will prompt for directory

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(grid)  # for arrow units
})

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Function to setup logo with automatic copying and resource path
setup_logo_resources <- function(results_dir) {
  # Check if logo exists in app www directory (primary source)
  app_logo <- "www/logo2.png"
  
  if (file.exists(app_logo)) {
    # Option 1: Try to set up resource path to results directory
    results_www <- file.path(results_dir, "www")
    
    # Create results www directory if it doesn't exist
    if (!dir.exists(results_www)) {
      dir.create(results_www, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Copy logo to results directory for proper serving
    results_logo <- file.path(results_www, "logo2.png")
    if (!file.exists(results_logo) || 
        file.info(app_logo)$mtime > file.info(results_logo)$mtime) {
      file.copy(app_logo, results_logo, overwrite = TRUE)
      cat("📁 Logo copied to results directory:", results_www, "\n")
    }
    
    # Set up resource path
    addResourcePath("logo", results_www)
    cat("📁 Logo resources configured from results directory\n")
    return("logo/logo2.png")
    
  } else if (file.exists("logo2.png")) {
    # Fallback: logo in app root directory
    dir.create("www", showWarnings = FALSE)
    file.copy("logo2.png", "www/logo2.png", overwrite = TRUE)
    cat("📁 Logo copied from app root to www directory\n")
    return("logo2.png")
    
  } else {
    cat("⚠️ Logo not found, using text fallback\n")
    return(NULL)
  }
}

# Function to get logo element
get_logo_element <- function(logo_src = NULL) {
  if (!is.null(logo_src)) {
    return(tags$img(src = logo_src, height = "80px", style = "margin-right: 25px;"))
  } else {
    # Fallback to text/emoji logo
    return(tags$span("🧬 PeakPrime", 
                     style = "font-size: 24px; font-weight: bold; color: #2c3e50; margin-right: 15px;"))
  }
}

# ==============================================================================
# DIRECTORY SELECTION
# ==============================================================================

# Get results directory from command line or user input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  results_dir <- args[1]
} else {
  # Interactive directory selection
  cat("=== PeakPrime Standalone Explorer ===\n")
  cat("Enter the path to your results directory:\n")
  cat("(e.g., /path/to/results/testboth or results/Class1)\n")
  results_dir <- readline("Results directory: ")
}

# Validate directory
if (!dir.exists(results_dir)) {
  stop("❌ Directory does not exist: ", results_dir)
}

# Convert to absolute path and store original working directory
results_dir <- normalizePath(results_dir, mustWork = TRUE)
original_wd <- getwd()

cat("📁 Using results directory:", results_dir, "\n")

# Check for preprocessed data
manifest_file <- file.path(results_dir, "data_manifest.rds")
if (!file.exists(manifest_file)) {
  cat("❌ No preprocessed data found in:", results_dir, "\n")
  cat("💡 Run preprocess_for_standalone.R first:\n")
  cat("   source('preprocess_for_standalone.R')\n")
  cat("   preprocess_peakprime_standalone('", results_dir, "', '/path/to/gtf')\n")
  stop("Preprocessing required")
}

# ==============================================================================
# DATA LOADING FROM RESULTS DIRECTORY
# ==============================================================================

# Load all preprocessed data from results directory
cat("🚀 Loading preprocessed data from:", results_dir, "\n")
start_time <- Sys.time()

manifest <- readRDS(file.path(results_dir, "data_manifest.rds"))
qc_data_raw <- readRDS(file.path(results_dir, "qc_data.rds"))
selected_peaks <- if (file.exists(file.path(results_dir, "peaks_data.rds"))) {
  readRDS(file.path(results_dir, "peaks_data.rds"))
} else NULL
gtf_data <- if (file.exists(file.path(results_dir, "gtf_data.rds"))) {
  readRDS(file.path(results_dir, "gtf_data.rds"))
} else NULL
coverage_index <- if (file.exists(file.path(results_dir, "coverage_index.rds"))) {
  readRDS(file.path(results_dir, "coverage_index.rds"))
} else NULL

# Load alignment summary data
alignment_summary <- if (file.exists(file.path(results_dir, "alignment_summary.rds"))) {
  readRDS(file.path(results_dir, "alignment_summary.rds"))
} else data.frame()

# Load all peaks data if available
all_peaks_data <- NULL
narrowpeak_file <- list.files(file.path(results_dir, "macs2_peaks"), 
                             pattern = "_peaks\\.narrowPeak$", 
                             full.names = TRUE)[1]
if (!is.na(narrowpeak_file) && file.exists(narrowpeak_file)) {
  cat("📈 Loading all detected peaks data...\n")
  all_peaks_raw <- fread(narrowpeak_file, 
                        header = FALSE, 
                        col.names = c("chr", "start", "end", "name", "score", 
                                     "strand", "fold_change", "pvalue", "qvalue", "summit_offset"))
  
  # Process peaks into gene-indexed structure
  if (nrow(all_peaks_raw) > 0) {
    # Use QC data to map peaks to genes
    qc_summary <- qc_data_raw
    
    # Create a mapping structure: gene_id -> list of all peaks
    all_peaks_data <- list()
    
    for (i in 1:nrow(qc_summary)) {
      gene_id <- qc_summary$gene_id[i]
      gene_chr <- qc_summary$gene_chr[i]
      gene_start <- qc_summary$gene_start[i]
      gene_end <- qc_summary$gene_end[i]
      
      # Find all peaks overlapping this gene
      gene_peaks <- all_peaks_raw[
        chr == gene_chr & 
        end >= gene_start & 
        start <= gene_end
      ]
      
      if (nrow(gene_peaks) > 0) {
        # Sort peaks by score (best first)
        gene_peaks <- gene_peaks[order(-score)]
        all_peaks_data[[gene_id]] <- gene_peaks
      }
    }
    
    cat("   Processed peaks for", length(all_peaks_data), "genes\n")
  }
}

# Convert data.table to regular data.frame to avoid filtering issues
qc_data <- as.data.frame(qc_data_raw)
cat("📋 Converted data.table to data.frame for reliable filtering\n")

load_time <- Sys.time() - start_time
cat("⚡ Loaded", nrow(qc_data), "genes in", format(load_time, digits = 2), "\n")

# Report alignment data availability
if (nrow(alignment_summary) > 0) {
  cat("🔍 Alignment data loaded:", nrow(alignment_summary), "records\n")
} else {
  cat("ℹ️ No alignment data available\n")
}

# Store source information for UI display
source_info <- paste0("📂 Source: ", basename(results_dir), 
                     "\n📄 GTF: ", basename(manifest$gtf_path %||% "Unknown"),
                     "\n⏱️ Processed: ", format(manifest$processed_time, "%Y-%m-%d %H:%M"))

# ==============================================================================
# ENHANCED PLOTTING FUNCTIONS
# ==============================================================================

#' Get gene coverage from preprocessed data
get_gene_coverage <- function(gene_id, coverage_index) {
  if (is.null(coverage_index) || !gene_id %in% names(coverage_index)) {
    return(data.table(position = numeric(0), coverage = numeric(0)))
  }
  
  # Return preprocessed coverage data
  coverage_index[[gene_id]]
}

#' Get gene structure from preprocessed GTF
extract_gene_info <- function(gene_id, gtf_data, max_transcripts = 20) {
  # Find all features for this gene (handle versioned gene IDs)
  gene_col <- as.character(mcols(gtf_data)$gene_id)
  gene_mask <- !is.na(gene_col) & 
                (gene_col == gene_id |
                 # Handle cases where one has version and other doesn't
                 sub("\\.\\d+$", "", gene_col) == gene_id |
                 sub("\\.\\d+$", "", gene_id) == gene_col)
  
  gene_features <- gtf_data[gene_mask]
  
  if (length(gene_features) == 0) {
    cat("Warning: No GTF features found for gene", gene_id, "\n")
    return(list(exons = GRanges(), cds = GRanges(), utrs = GRanges(), transcripts = character(0)))
  }
  
  # Split by feature type - following MakePlots_new.R approach
  type_col <- mcols(gene_features)$type
  if (is.null(type_col)) {
    cat("Warning: No type column found in GTF data\n")
    return(list(exons = gene_features, cds = GRanges(), utrs = GRanges(), transcripts = "Gene_Region"))
  }
  
  # Extract all feature types like MakePlots_new.R
  exons <- gene_features[type_col == "exon"]
  cds <- gene_features[type_col == "CDS"]
  utrs <- gene_features[type_col %in% c("three_prime_utr", "five_prime_utr")]
  
  # Get unique transcripts from exons (includes both coding and non-coding)
  if (length(exons) > 0) {
    tx_col <- mcols(exons)$transcript_id
    if (!is.null(tx_col)) {
      transcripts <- unique(na.omit(as.character(tx_col)))
      # Note: transcript limiting will be handled by max_isoforms parameter in plotting function
      # Keep all transcripts available here for flexibility
    } else {
      transcripts <- "Transcript_1"
    }
  } else {
    transcripts <- "No_Exons"
  }
  
  if (length(transcripts) == 0) transcripts <- "Unknown_Transcript"
  
  # Apply transcript limit if specified
  if (length(transcripts) > max_transcripts) {
    cat("Limiting from", length(transcripts), "to", max_transcripts, "transcripts for visualization\n")
    transcripts <- head(transcripts, max_transcripts)
  }
  
  # Get gene strand information
  gene_strand <- if (length(gene_features) > 0) {
    strand_info <- as.character(strand(gene_features))[1]
    if (is.na(strand_info) || strand_info == "*") "+" else strand_info
  } else {
    "+"  # Default to positive strand
  }
  
  cat("Found", length(transcripts), "transcripts for gene", gene_id, "(including non-coding)\n")
  
  list(exons = exons, cds = cds, utrs = utrs, transcripts = transcripts, gene_strand = gene_strand)
}

#' Get all detected peaks for a gene
get_all_peaks <- function(gene_id, all_peaks_data, max_peaks = 5) {
  if (is.null(all_peaks_data) || !gene_id %in% names(all_peaks_data)) {
    return(data.frame())
  }
  
  # Return top peaks (already sorted by score)
  peaks <- all_peaks_data[[gene_id]]
  head(peaks, max_peaks)
}

#' Create comprehensive gene plot with full features
plot_gene_comprehensive <- function(gene_id, qc_data, all_peaks_data, yaxis_mode = "depth", 
                                  show_primers = TRUE, max_isoforms = 5,
                                  show_peak_region = TRUE, show_all_peaks = FALSE, 
                                  max_peaks = 3) {
  
  # Get gene info from QC data - ensure fresh filtering
  gene_info <- qc_data[qc_data$gene_id == gene_id, , drop = FALSE]
  if (nrow(gene_info) == 0) {
    return(ggplot() + 
           annotate("text", x = 0.5, y = 0.5, label = "Gene not found", size = 6) +
           theme_void() + xlim(0, 1) + ylim(0, 1))
  }
  
  # Ensure we get exactly one row
  if (nrow(gene_info) > 1) {
    warning("Multiple rows found for gene ", gene_id, ". Using first row.")
    gene_info <- gene_info[1, , drop = FALSE]
  }
  
  # Extract coordinates safely
  peak_start_coord <- as.numeric(gene_info$peak_start[1])
  peak_end_coord <- as.numeric(gene_info$peak_end[1])
  gene_start_coord <- as.numeric(gene_info$gene_start[1])
  gene_end_coord <- as.numeric(gene_info$gene_end[1])
  
  # Debug: print gene coordinates for troubleshooting
  cat("Gene:", gene_id, "\n")
  cat("QC data peak coords:", peak_start_coord, "-", peak_end_coord, "\n")
  cat("Gene coords:", gene_start_coord, "-", gene_end_coord, "\n")
  
  # Store the trimmed coordinates from QC data (for highlighting)
  trimmed_peak_start <- peak_start_coord
  trimmed_peak_end <- peak_end_coord
  
  # Always get original MACS2 peaks for display in peak bars
  original_peaks <- get_all_peaks(gene_id, all_peaks_data, max_peaks = 10)
  
  # Determine what to use for highlighting vs peak display
  peak_was_trimmed <- (peak_start_coord != gene_start_coord || peak_end_coord != gene_end_coord)
  
  if (peak_was_trimmed) {
    cat("Trimming detected - Using separate coordinates:\n")
    cat("  Trimmed region (for highlighting):", trimmed_peak_start, "-", trimmed_peak_end, "\n")
    cat("  Gene boundaries:", gene_start_coord, "-", gene_end_coord, "\n")
    if (nrow(original_peaks) > 0) {
      cat("  Original MACS2 peak (for bars):", original_peaks$start[1], "-", original_peaks$end[1], "\n")
    }
  } else {
    cat("No trimming detected - peak coords match gene coords or using MACS2 coords\n")
    if (nrow(original_peaks) > 0) {
      cat("Using MACS2 peak coords:", original_peaks$start[1], "-", original_peaks$end[1], "\n")
      # For highlighting, use MACS2 coordinates when no trimming was applied
      trimmed_peak_start <- original_peaks$start[1]  
      trimmed_peak_end <- original_peaks$end[1]
    } else {
      cat("No MACS2 peaks found, using QC data coordinates for both\n")
    }
  }
  
  # Note: peak_start_coord and peak_end_coord now represent the actual peak region
  # These should be used for highlighting (primer design region)
  # All MACS2 peaks will be displayed via the all_peaks system
  
  # Get coverage data
  coverage_df <- get_gene_coverage(gene_id, coverage_index)
  
  # Extract gene structure
  structure <- extract_gene_info(gene_id, gtf_data, max_transcripts = max_isoforms)
  
  # Determine plot coordinates - ensure peak is always included
  if (nrow(coverage_df) > 0) {
    plot_start <- min(coverage_df$position)
    plot_end <- max(coverage_df$position)
  } else if (length(structure$exons) > 0) {
    plot_start <- min(start(structure$exons))
    plot_end <- max(end(structure$exons))
  } else {
    # Fallback to QC data coordinates
    plot_start <- gene_start_coord
    plot_end <- gene_end_coord
  }
  
  # CRITICAL FIX: Find original MACS2 peak coordinates and expand plot window
  original_macs2_start <- peak_start_coord  # Default to trimmed coordinates
  original_macs2_end <- peak_end_coord
  
  # Find the original MACS2 peak that overlaps with the selected region
  if (!is.na(peak_start_coord) && !is.na(peak_end_coord)) {
    all_peaks_for_gene <- get_all_peaks(gene_id, all_peaks_data, max_peaks)
    
    if (nrow(all_peaks_for_gene) > 0) {
      for (i in 1:nrow(all_peaks_for_gene)) {
        macs2_start <- all_peaks_for_gene$start[i]
        macs2_end <- all_peaks_for_gene$end[i]
        
        # Check for overlap between MACS2 peak and selected region
        if (!(macs2_end < peak_start_coord || macs2_start > peak_end_coord)) {
          original_macs2_start <- macs2_start
          original_macs2_end <- macs2_end
          cat("Found overlapping MACS2 peak for plot expansion:", macs2_start, "-", macs2_end, "\n")
          break
        }
      }
    }
    
    # Expand plot window to include the ORIGINAL MACS2 peak (not just trimmed region)
    plot_start <- min(plot_start, original_macs2_start)
    plot_end <- max(plot_end, original_macs2_end)
    cat("Expanded plot window to include MACS2 peak: ", plot_start, "-", plot_end, "\n")
  }
  
  # CRITICAL FIX: Expand plot window to include ALL UTRs (fixes missing UTR display)
  if (length(structure$utrs) > 0) {
    utr_start_min <- min(start(structure$utrs))
    utr_end_max <- max(end(structure$utrs))
    original_start <- plot_start
    original_end <- plot_end
    plot_start <- min(plot_start, utr_start_min)
    plot_end <- max(plot_end, utr_end_max)
    if (plot_start != original_start || plot_end != original_end) {
      cat("Expanded plot window to include UTRs: ", plot_start, "-", plot_end, "\\n")
      cat("UTR range: ", utr_start_min, "-", utr_end_max, "\\n")
    }
  }
  
  # Create main title and subtitle  
  gene_chr <- as.character(gene_info$gene_chr[1])
  gene_strand <- as.character(gene_info$gene_strand[1])
  
  main_title <- sprintf("%s  %s:%d-%d (%s)", 
                       gene_id, 
                       gene_chr,
                       plot_start, 
                       plot_end, 
                       ifelse(is.na(gene_strand), "*", gene_strand))
  
  qc_subtitle <- sprintf("Peak Score: %.2f | Exonic Fraction: %.3f | Selected: %s",
                        ifelse(is.na(gene_info$best_peak_score[1]), 0, gene_info$best_peak_score[1]),
                        ifelse(is.na(gene_info$exonic_fraction[1]), 0, gene_info$exonic_fraction[1]),
                        ifelse(is.na(gene_info$final_selection[1]) || !gene_info$final_selection[1], "NO", "YES"))
  
  # ==== COVERAGE PLOT ====
  p_cov <- ggplot()
  
  if (nrow(coverage_df) > 0) {
    # Prepare coverage data based on y-axis mode
    if (yaxis_mode == "percent") {
      max_cov <- max(coverage_df$coverage, na.rm = TRUE)
      coverage_df$y_value <- (coverage_df$coverage / max_cov) * 100
      y_label <- "Coverage (% of max)"
      y_limits <- c(0, 100)
    } else if (yaxis_mode == "log10") {
      coverage_df$y_value <- log10(pmax(coverage_df$coverage, 1))
      y_label <- "Coverage (log10)"
      y_limits <- NULL
    } else {  # depth
      coverage_df$y_value <- coverage_df$coverage
      y_label <- "Coverage (depth)"
      y_limits <- NULL
    }
    
    # Create coverage area plot
    p_cov <- p_cov +
      geom_area(data = coverage_df, 
               aes(x = position, y = y_value), 
               fill = "steelblue", alpha = 0.6) +
      geom_line(data = coverage_df,
               aes(x = position, y = y_value),
               color = "darkblue", size = 0.5)
    
    if (!is.null(y_limits)) {
      p_cov <- p_cov + scale_y_continuous(name = y_label, limits = y_limits)
    } else {
      p_cov <- p_cov + scale_y_continuous(name = y_label)
    }
  } else {
    # No coverage data - create empty plot
    p_cov <- p_cov +
      annotate("text", x = (plot_start + plot_end) / 2, y = 0.5, 
              label = "No coverage data available", size = 4, color = "gray50") +
      scale_y_continuous(name = "Coverage", limits = c(0, 1))
  }
  
  # Add TRIMMED peak region highlight (primer design region)
  # This highlights only the selected/trimmed peak, not the full MACS2 peak
  if (show_peak_region && !is.na(trimmed_peak_start) && !is.na(trimmed_peak_end)) {
    p_cov <- p_cov + 
      annotate("rect", 
              xmin = trimmed_peak_start, 
              xmax = trimmed_peak_end,
              ymin = -Inf, ymax = Inf, 
              alpha = 0.3, fill = "yellow", color = "orange", linewidth = 0.8) +
      annotate("text", x = (trimmed_peak_start + trimmed_peak_end) / 2, y = Inf, 
              label = "TRIMMED PEAK\n(Primer Region)", vjust = 1.2, size = 3, 
              color = "darkred", fontface = "bold")
  }
  
  p_cov <- p_cov +
    scale_x_continuous(limits = c(plot_start, plot_end)) +
    labs(title = main_title, subtitle = qc_subtitle, x = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title.position = "plot",
          axis.title.x = element_blank())
  
  # ==== GENE STRUCTURE PLOT ====
  p_feat <- ggplot()
  
  # Track layout - ensure we always have at least one track for primers and peaks
  transcript_tracks <- head(structure$transcripts, max_isoforms)
  n_transcript_tracks <- length(transcript_tracks)
  
  # GTF data should be mandatory now - but provide fallback
  if (n_transcript_tracks == 0) {
    warning("No transcript data found for gene ", gene_id, ". GTF annotations may be missing.")
    n_transcript_tracks <- 1
    transcript_tracks <- paste("Gene_Region (", gene_id, ")")
  }
  
  primer_track <- n_transcript_tracks + 1
  peak_track <- primer_track + 1
  
  # Add TRIMMED peak region highlight to gene structure (primer design region)
  if (show_peak_region && !is.na(trimmed_peak_start) && !is.na(trimmed_peak_end)) {
    p_feat <- p_feat + 
      annotate("rect",
              xmin = trimmed_peak_start,
              xmax = trimmed_peak_end, 
              ymin = 0.5, ymax = peak_track + 0.5,
              alpha = 0.3, fill = "yellow", color = "orange", linewidth = 0.8)
  }
  
  # STEP 1: Draw ALL intron lines first (background layer)
  cat("Drawing intron lines for all transcripts...\n")
  for (i in seq_along(transcript_tracks)) {
    tx_id <- transcript_tracks[i]
    y_pos <- i
    
    if (!(tx_id %in% c("Gene_Region", "No_Exons", "Unknown_Transcript")) && 
        length(structure$exons) > 0) {
      
      # Get exons for this transcript
      tx_col <- mcols(structure$exons)$transcript_id
      tx_mask <- !is.na(tx_col) & as.character(tx_col) == tx_id
      tx_exons <- structure$exons[tx_mask]
      
      if (length(tx_exons) > 0) {
        # Use distinct colors for each transcript
        tx_colors <- c("black", "darkblue", "darkgreen", "purple", "darkred")
        tx_color <- tx_colors[((i - 1) %% length(tx_colors)) + 1]
        
        cat("  Drawing introns for", tx_id, "at y =", y_pos, "color =", tx_color, "\n")
        
        # Draw true intron segments between consecutive exons
        if (length(tx_exons) > 1) {
          # Sort exons by position
          exon_df <- data.frame(start = start(tx_exons), end = end(tx_exons))
          exon_df <- exon_df[order(exon_df$start, exon_df$end), , drop = FALSE]
          
          # Create intron segments between consecutive exons
          introns <- data.frame(
            x    = head(exon_df$end, -1),
            xend = tail(exon_df$start, -1),
            y    = y_pos,
            yend = y_pos
          )
          
          # Only draw introns where there's actual gap between exons
          introns <- introns[introns$x < introns$xend, , drop = FALSE]
          
          if (nrow(introns) > 0) {
            cat("    Drawing", nrow(introns), "intron segments\n")
            p_feat <- p_feat +
              geom_segment(
                data = introns,
                aes(x = x, xend = xend, y = y, yend = yend),
                inherit.aes = FALSE,
                linewidth = 1.2,
                color = tx_color
              )
          }
        } else if (length(tx_exons) == 1) {
          # Single-exon transcript: draw a short baseline for visual consistency
          tx_start <- start(tx_exons[1])
          tx_end <- end(tx_exons[1])
          cat("    Single exon transcript - drawing baseline\n")
          p_feat <- p_feat +
            annotate(
              "segment",
              x = tx_start, xend = tx_end,
              y = y_pos, yend = y_pos,
              color = tx_color, linewidth = 1.0
            )
        }
      }
    }
  }
  
  # STEP 2: Draw transcript structures (exons/CDS on top)
  cat("Drawing exon structures for all transcripts...\n")
  for (i in seq_along(transcript_tracks)) {
    tx_id <- transcript_tracks[i]
    y_pos <- i
    
    cat("Drawing transcript", i, ":", tx_id, "at position", y_pos, "\n")
    
    if (tx_id %in% c("Gene_Region", "No_Exons", "Unknown_Transcript") || 
        length(structure$exons) == 0) {
      # No GTF data - draw a simple gene region line
      p_feat <- p_feat + 
        geom_segment(aes(x = plot_start, xend = plot_end, y = y_pos, yend = y_pos), 
                    color = "gray40", linewidth = 1) +
        geom_rect(aes(xmin = plot_start, xmax = plot_end, 
                     ymin = y_pos - 0.1, ymax = y_pos + 0.1),
                 fill = "lightgray", color = "black", linewidth = 0.1)
    } else {
      # GTF data available - draw transcript structures
      # Get exons for this specific transcript
      if (length(structure$exons) > 0) {
        tx_col <- mcols(structure$exons)$transcript_id
        tx_mask <- !is.na(tx_col) & as.character(tx_col) == tx_id
        tx_exons <- structure$exons[tx_mask]
        
        # Debug transcript filtering
        cat("  Transcript", tx_id, ": matching exons =", length(tx_exons), "\n")
      } else {
        tx_exons <- GRanges()
      }
      
      # Get CDS for this transcript
      if (length(structure$cds) > 0) {
        cds_col <- mcols(structure$cds)$transcript_id  
        cds_mask <- !is.na(cds_col) & as.character(cds_col) == tx_id
        tx_cds <- structure$cds[cds_mask]
      } else {
        tx_cds <- GRanges()
      }
      
      if (length(tx_exons) > 0) {
        # Exons and CDS are drawn below - intron lines already drawn above
        
        # FIRST: Draw all exons as light gray background (ALL transcribed regions)
        # This includes both coding and non-coding transcript exons
        if (length(tx_exons) > 0) {
          exon_df <- data.frame(
            xmin = start(tx_exons),
            xmax = end(tx_exons),
            ymin = y_pos - 0.15,
            ymax = y_pos + 0.15
          )
          p_feat <- p_feat + 
            geom_rect(data = exon_df,
                     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     fill = "lightgray", color = "gray", linewidth = 0.2, alpha = 0.7)
          
          cat("    Drew", length(tx_exons), "exon regions (transcribed areas)\\n")
        }
        
        # SECOND: Draw CDS (coding sequences) on top - ONLY for protein-coding transcripts
        if (length(tx_cds) > 0) {
          cds_df <- data.frame(
            xmin = start(tx_cds),
            xmax = end(tx_cds),
            ymin = y_pos - 0.3,
            ymax = y_pos + 0.3
          )
          p_feat <- p_feat +
            geom_rect(data = cds_df,
                     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     fill = "steelblue", color = "black", linewidth = 0.3)
          
          cat("    Drew", length(tx_cds), "CDS regions (protein-coding areas)\\n")
        } else {
          cat("    No CDS regions - this is a non-coding transcript\\n")
        }
        
        # THIRD: Draw UTR regions (on top of exons, narrower than CDS)
        if (length(structure$utrs) > 0) {
          utr_col <- mcols(structure$utrs)$transcript_id  
          utr_mask <- !is.na(utr_col) & as.character(utr_col) == tx_id
          tx_utrs <- structure$utrs[utr_mask]
          
          if (length(tx_utrs) > 0) {
            utr_types <- mcols(tx_utrs)$type
            
            # 5' UTRs in light blue
            utr_5prime <- tx_utrs[utr_types == "five_prime_utr"]
            if (length(utr_5prime) > 0) {
              # Enhance visibility for small UTRs
              utr5_starts <- start(utr_5prime)
              utr5_ends <- end(utr_5prime)
              utr5_widths <- utr5_ends - utr5_starts + 1
              
              # Add minimum visual width for very small UTRs (< 100bp)
              min_visual_width <- 50
              for (i in seq_along(utr5_starts)) {
                if (utr5_widths[i] < min_visual_width) {
                  center <- (utr5_starts[i] + utr5_ends[i]) / 2
                  utr5_starts[i] <- center - min_visual_width/2
                  utr5_ends[i] <- center + min_visual_width/2
                }
              }
              
              utr5_df <- data.frame(
                xmin = utr5_starts,
                xmax = utr5_ends,
                ymin = y_pos - 0.25,  # Make thicker for visibility
                ymax = y_pos + 0.25
              )
              p_feat <- p_feat +
                geom_rect(data = utr5_df,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         fill = "lightblue", alpha = 0.95, color = "darkblue", linewidth = 0.5)
              
              cat("    Drew", length(utr_5prime), "5' UTR regions (light blue) - Width:", paste(utr5_widths, "bp"), "\\n")
            }
            
            # 3' UTRs in orange
            utr_3prime <- tx_utrs[utr_types == "three_prime_utr"]
            if (length(utr_3prime) > 0) {
              # Enhance visibility for small UTRs
              utr_starts <- start(utr_3prime)
              utr_ends <- end(utr_3prime)
              utr_widths <- utr_ends - utr_starts + 1
              
              # Add minimum visual width for very small UTRs (< 100bp)
              min_visual_width <- 50  # minimum pixels for visibility
              for (i in seq_along(utr_starts)) {
                if (utr_widths[i] < min_visual_width) {
                  center <- (utr_starts[i] + utr_ends[i]) / 2
                  utr_starts[i] <- center - min_visual_width/2
                  utr_ends[i] <- center + min_visual_width/2
                }
              }
              
              utr3_df <- data.frame(
                xmin = utr_starts,
                xmax = utr_ends,
                ymin = y_pos - 0.25,  # Make thicker for visibility
                ymax = y_pos + 0.25
              )
              p_feat <- p_feat +
                geom_rect(data = utr3_df,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         fill = "orange", alpha = 0.95, color = "darkorange", linewidth = 0.5)
              
              cat("    Drew", length(utr_3prime), "3' UTR regions (orange) - Original width:", paste(utr_widths, "bp"), "\\n")
            }
          }
        }
        
        cat("Drew", length(tx_exons), "exons and", length(tx_cds), "CDS for", tx_id, "\n")
      } else {
        cat("No exons found for transcript", tx_id, "\n")
        # Draw fallback line
        p_feat <- p_feat + 
          geom_segment(aes(x = plot_start, xend = plot_end, y = y_pos, yend = y_pos), 
                      color = "gray60", linewidth = 0.5)
      }
    }
  }
  
  # Add primer annotations - use trimmed coordinates for primer positioning
  peak_start_val <- trimmed_peak_start
  peak_end_val <- trimmed_peak_end
  
  cat("Checking primers: show_primers =", show_primers, "peak_start =", peak_start_val, "peak_end =", peak_end_val, "\n")
  
  if (show_primers && !is.na(peak_start_val) && !is.na(peak_end_val)) {
    # Calculate primer position at TRIMMED peak center (primer design region)
    selected_peak_center <- (peak_start_val + peak_end_val) / 2
    
    # Determine primer direction based on gene strand (towards 3' end)
    primer_direction <- if (!is.null(gene_strand) && length(gene_strand) > 0 && gene_strand == "-") {
      "left"  # - strand: 3' is left
    } else {
      "right"  # + strand: 3' is right (default)
    }
    
    cat("Drawing single primer at SELECTED peak position:", selected_peak_center, "pointing", primer_direction, "(towards 3' end) on track", primer_track, "\n")
    
    # Single primer arrow pointing towards 3' end - ALWAYS use selected peak
    if (primer_direction == "left") {
      # Gene on minus strand: primer points left (towards 3' end)
      p_feat <- p_feat +
        geom_segment(aes(x = selected_peak_center + 30, xend = selected_peak_center - 30,
                        y = primer_track, yend = primer_track), 
                    arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
                    color = "red", linewidth = 2)
    } else {
      # Gene on plus strand: primer points right (towards 3' end)
      p_feat <- p_feat +
        geom_segment(aes(x = selected_peak_center - 30, xend = selected_peak_center + 30, 
                        y = primer_track, yend = primer_track),
                    arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
                    color = "red", linewidth = 2)
    }
  } else {
    cat("Primers not drawn: show_primers =", show_primers, 
        "coordinates valid =", !is.na(peak_start_val) && !is.na(peak_end_val), "\n")
  }
  
  # Draw peak regions on dedicated peak track
  if (show_all_peaks) {
    # Draw all detected peaks as rectangles on bottom track only
    all_peaks <- get_all_peaks(gene_id, all_peaks_data, max_peaks)
    
    # Always ensure selected peak is included if it exists
    selected_peak_drawn <- FALSE
    
    if (nrow(all_peaks) > 0) {
      cat("Drawing", nrow(all_peaks), "detected peaks on track", peak_track, "\n")
      
      # First pass: Draw all detected peaks
      for (i in 1:nrow(all_peaks)) {
        # Use immediate evaluation to avoid ggplot lazy evaluation issues
        peak_start_coord <- all_peaks$start[i]
        peak_end_coord <- all_peaks$end[i]
        peak_score_val <- all_peaks$score[i]
        
        # Check if this MACS2 peak overlaps with the selected/trimmed peak region
        # This handles cases where the original MACS2 peak is larger than the trimmed region
        is_selected_peak <- (!is.na(peak_start_val) && !is.na(peak_end_val) &&
                            !(peak_end_coord < peak_start_val || peak_start_coord > peak_end_val))
        
        if (is_selected_peak) {
          selected_peak_drawn <- TRUE
          # Selected peak - prominent style
          peak_fill <- "darkred"
          alpha_value <- 0.9
          border_col <- "black"
          border_size <- 0.5
        } else {
          # Non-selected peaks - muted style
          peak_fill <- "lightgray"
          alpha_value <- 0.6
          border_col <- "gray"
          border_size <- 0.2
        }
        
        # Y-position with slight offset for multiple peaks
        y_pos_offset <- (i - 1) * 0.1 - 0.2
        
        # Use annotate to avoid lazy evaluation issues
        p_feat <- p_feat +
          annotate("rect", 
                  xmin = peak_start_coord, xmax = peak_end_coord,
                  ymin = peak_track - 0.4 + y_pos_offset, 
                  ymax = peak_track + 0.4 + y_pos_offset,
                  fill = peak_fill, alpha = alpha_value, 
                  color = border_col, linewidth = border_size)
      }
    }
    
    # If selected peak wasn't found in detected peaks, draw the trimmed region separately  
    if (!selected_peak_drawn && !is.na(peak_start_val) && !is.na(peak_end_val)) {
      cat("Drawing trimmed peak region separately (not overlapping with MACS2 peaks)\n")
      p_feat <- p_feat +
        annotate("rect", 
                xmin = peak_start_val, xmax = peak_end_val,
                ymin = peak_track - 0.4, ymax = peak_track + 0.4,
                fill = "darkred", alpha = 0.9, 
                color = "black", linewidth = 0.5)
    }
  } else {
    # Draw only selected peak using ORIGINAL MACS2 coordinates (not trimmed)
    # Find the original MACS2 peak that overlaps with our selected region
    all_peaks <- get_all_peaks(gene_id, all_peaks_data, max_peaks)
    selected_macs2_peak <- NULL
    
    if (nrow(all_peaks) > 0 && !is.na(peak_start_val) && !is.na(peak_end_val)) {
      # Find MACS2 peak that overlaps with selected region
      for (i in 1:nrow(all_peaks)) {
        peak_start_coord <- all_peaks$start[i]
        peak_end_coord <- all_peaks$end[i]
        
        # Check for overlap between MACS2 peak and selected region
        if (!(peak_end_coord < peak_start_val || peak_start_coord > peak_end_val)) {
          selected_macs2_peak <- list(start = peak_start_coord, end = peak_end_coord)
          break
        }
      }
    }
    
    if (!is.null(selected_macs2_peak)) {
      cat("Drawing selected MACS2 peak (original coordinates) from", 
          selected_macs2_peak$start, "to", selected_macs2_peak$end, "on track", peak_track, "\n")
      
      p_feat <- p_feat +
        annotate("rect", 
                xmin = selected_macs2_peak$start, xmax = selected_macs2_peak$end,
                ymin = peak_track - 0.4, ymax = peak_track + 0.4,
                fill = "darkred", alpha = 0.9, 
                color = "black", linewidth = 0.5)
    } else if (!is.na(peak_start_val) && !is.na(peak_end_val)) {
      # Fallback: if no MACS2 peak found, draw trimmed region
      cat("No overlapping MACS2 peak found, drawing trimmed region from", 
          peak_start_val, "to", peak_end_val, "on track", peak_track, "\n")
      
      p_feat <- p_feat +
        annotate("rect", 
                xmin = peak_start_val, xmax = peak_end_val,
                ymin = peak_track - 0.4, ymax = peak_track + 0.4,
                fill = "darkred", alpha = 0.9, 
                color = "black", linewidth = 0.5)
    } else {
      cat("Peak region not drawn: coordinates invalid\n")
    }
  }
  
  # Format gene structure plot
  peak_label <- if (show_all_peaks) "All MACS2 Peaks" else "Trimmed Peak"
  track_labels <- c(transcript_tracks, "Primers", peak_label)
  track_positions <- seq_len(peak_track)
  
  p_feat <- p_feat +
    scale_x_continuous(limits = c(plot_start, plot_end)) +
    scale_y_continuous(limits = c(0.5, peak_track + 0.5),
                      breaks = track_positions,
                      labels = track_labels) +
    labs(x = sprintf("Position (%s)", ifelse(is.na(gene_chr), "chr?", gene_chr)), y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 8))
  
  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p_cov, p_feat, ncol = 1, heights = c(2, 1)))
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(p_cov, p_feat, ncol = 1, rel_heights = c(2, 1)))
  } else {
    # Fallback - return just coverage plot with note
    p_cov <- p_cov + 
      labs(caption = "Install 'patchwork' or 'cowplot' for combined gene structure view")
    return(p_cov)
  }
}

# ==============================================================================
# SHINY UI SETUP
# ==============================================================================

# Setup logo resources early (before UI definition)
logo_src <- setup_logo_resources(results_dir)

ui <- fluidPage(
  titlePanel(div(
    style = "display: flex; align-items: center; margin-bottom: 20px;",
    get_logo_element(logo_src),  # Use configured logo source
    span("PeakPrime Explorer", style = "font-size: 28px; font-weight: bold;")
  )),
  
  tags$head(
    tags$style(HTML("
      .sidebar { background-color: #f8f9fa; padding: 15px; }
      .performance-badge { 
        background: linear-gradient(45deg, #28a745, #20c997);
        color: white; padding: 5px 10px; border-radius: 15px;
        font-size: 0.8em; font-weight: bold;
      }
      .source-info {
        background: linear-gradient(45deg, #6f42c1, #007bff);
        color: white; padding: 8px 12px; border-radius: 10px;
        font-size: 0.85em; margin-bottom: 15px;
      }
      .feature-badge {
        background: linear-gradient(45deg, #007bff, #6f42c1);
        color: white; padding: 5px 10px; border-radius: 15px;
        font-size: 0.8em; font-weight: bold;
      }
      .title-panel .container-fluid {
        padding: 15px 15px 0px 15px;
      }
      .title-panel h1 {
        margin-bottom: 10px;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      width = 3,
      
      # Performance indicator
      div(class = "performance-badge", 
          paste("⚡ Loaded in", format(load_time, digits = 2))),
      br(),
      div(class = "source-info", source_info),
      br(),
      
      # Gene selection
      h4("🧬 Gene Selection"),
      selectInput("gene_id", "Choose Gene:",
                 choices = qc_data$gene_id,
                 selected = qc_data$gene_id[1]),
      
      checkboxInput("filter_selected", 
                   "Show only selected genes", 
                   value = TRUE),
      
      textOutput("gene_count_info"),
      br(),
      
      # Plot controls
      h4("📊 Plot Controls"),
      
      radioButtons("yaxis_mode", "Y-axis Scale:",
                  choices = list(
                    "Absolute depth" = "depth",
                    "Percentage of max" = "percent",
                    "Log10 scale" = "log10"
                  ),
                  selected = "depth"),
      
      checkboxInput("show_primers", "Show primer positions", value = TRUE),
      checkboxInput("show_peak_region", "Highlight selected peak", value = TRUE),
      
      # Peak display options
      checkboxInput("show_all_peaks", "Show all MACS2 detected peaks", value = FALSE),
      helpText("Note: Yellow highlighting shows trimmed peak (primer region)"),
      conditionalPanel(
        condition = "input.show_all_peaks",
        numericInput("max_peaks", "Max peaks to display:",
                    value = 3, min = 1, max = 8, step = 1)
      ),
      
      numericInput("max_isoforms", "Max transcript isoforms:",
                  value = 15, min = 1, max = 25, step = 1)
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        tabPanel("Gene Plot",
                plotOutput("comprehensive_plot", height = "650px"),
                br(),
                downloadButton("download_plot_png", "Download Plot as PNG")),
        
        tabPanel("Gene Data",
                h4("Gene Quality Control Summary"),
                DT::dataTableOutput("gene_table")),
        
        tabPanel("Performance Info",
                h4("Performance & Data Summary"),
                verbatimTextOutput("performance_info"),
                
                h4("Available Data"),
                verbatimTextOutput("data_summary")),
        
        tabPanel("Primer Alignment",
                h4("🔍 Primer Alignment Summary"),
                br(),
                conditionalPanel(
                  condition = "output.alignment_available",
                  p("Showing alignment results for the selected gene from the main panel."),
                  fluidRow(
                    column(6,
                      h5("Filter by Primer Index:"),
                      selectInput("alignment_primer_filter",
                                 label = NULL,
                                 choices = c("All Primers" = "all"),
                                 selected = "all",
                                 multiple = FALSE)
                    ),
                    column(6,
                      h5("Show only perfect matches:"),
                      checkboxInput("alignment_perfect_only", 
                                   label = "No mismatches", 
                                   value = FALSE),
                      br(),
                      downloadButton("download_alignment_csv", "Download Alignment CSV")
                    )
                  ),
                  br(),
                  h5("Alignment Summary Statistics:"),
                  verbatimTextOutput("alignment_stats"),
                  br(),
                  h5("Detailed Alignment Results:"),
                  DT::dataTableOutput("alignment_table")
                ),
                
                conditionalPanel(
                  condition = "!output.alignment_available",
                  wellPanel(
                    h5("⚠️ No Alignment Data Available"),
                    p("Primer alignment analysis was not performed or results are not available for this dataset."),
                    p("To generate alignment data, run the PeakPrime pipeline with primer alignment enabled."),
                    tags$code("nextflow run main.nf --bam sample.bam --genes targets.txt --transcriptome_qc")
                  )
                ))
      )
    )
  )
)

# ==============================================================================
# SHINY SERVER  
# ==============================================================================

server <- function(input, output, session) {
  
  # Reactive gene filtering
  filtered_genes <- reactive({
    if (input$filter_selected) {
      qc_data[qc_data$final_selection == TRUE, ]
    } else {
      qc_data
    }
  })
  
  # Update gene choices when filter changes
  observe({
    genes <- filtered_genes()
    updateSelectInput(session, "gene_id", 
                     choices = setNames(genes$gene_id, genes$gene_id),
                     selected = input$gene_id)
  })
  
  # Gene count info
  output$gene_count_info <- renderText({
    total <- nrow(qc_data)
    selected <- sum(qc_data$final_selection, na.rm = TRUE)
    showing <- nrow(filtered_genes())
    
    paste("Showing", showing, "of", total, "genes",
          sprintf("(%d selected)", selected))
  })
  
  # Main comprehensive plot
  output$comprehensive_plot <- renderPlot({
    req(input$gene_id)
    
    plot_gene_comprehensive(
      gene_id = input$gene_id,
      qc_data = qc_data,
      all_peaks_data = all_peaks_data,
      yaxis_mode = input$yaxis_mode,
      show_primers = input$show_primers,
      max_isoforms = input$max_isoforms,
      show_peak_region = input$show_peak_region,
      show_all_peaks = input$show_all_peaks,
      max_peaks = if(!is.null(input$max_peaks)) input$max_peaks else 3
    )
  })
  
  # Gene data table
  output$gene_table <- DT::renderDataTable({
    req(input$gene_id)
    
    gene_data <- qc_data[qc_data$gene_id == input$gene_id, ]
    
    # Transpose for better readability
    if (nrow(gene_data) > 0) {
      gene_transposed <- data.frame(
        Metric = names(gene_data),
        Value = as.character(unlist(gene_data[1, ]))
      )
      
      DT::datatable(gene_transposed, 
                   options = list(pageLength = 15, dom = 't'),
                   rownames = FALSE)
    } else {
      DT::datatable(data.frame(Message = "Gene not found"))
    }
  })
  
  # Performance information
  output$performance_info <- renderText({
    paste(
      "=== PERFORMANCE METRICS ===",
      sprintf("App startup time: %s", format(load_time, digits = 2)),
      sprintf("Preprocessed at: %s", format(manifest$processed_time)),
      sprintf("Total genes loaded: %d", nrow(qc_data)),
      "",
      "=== SPEED OPTIMIZATIONS ===", 
      "✓ All data preloaded into memory",
      "✓ Binary RDS format for instant access",
      "✓ Preprocessed coverage data",
      "✓ Optimized GTF structures", 
      "✓ Reactive filtering with cached data",
      "",
      "=== FEATURE COMPLETENESS ===",
      "✓ Multi-transcript gene models", 
      "✓ Coverage tracks with multiple scaling options",
      "✓ Peak region highlighting",
      "✓ Primer position annotations",
      "✓ Interactive gene filtering",
      sep = "\n"
    )
  })
  
  # Data summary
  output$data_summary <- renderText({
    paste(
      "=== LOADED DATA COMPONENTS ===",
      sprintf("QC Summary: %d genes", nrow(qc_data)),
      sprintf("Selected Peaks: %s", ifelse(is.null(selected_peaks), "Not available", 
                                          paste(nrow(selected_peaks), "peaks"))),
      sprintf("GTF Annotations: %s", ifelse(is.null(gtf_data), "Not available",
                                           paste(length(gtf_data), "features"))),
      sprintf("Coverage Index: %s", ifelse(is.null(coverage_index), "Not available",
                                          paste(length(coverage_index), "genes"))),
      "",
      "=== FILE SOURCES ===",
      sprintf("QC file: %s", basename(manifest$files$qc_summary %||% "missing")),
      sprintf("BigWig file: %s", basename(manifest$files$bigwig %||% "missing")),
      sprintf("GTF file: %s", basename(manifest$files$gtf %||% "missing")),
      sprintf("Peak file: %s", basename(manifest$files$selected_peaks %||% "missing")),
      "",
      ifelse(is.null(coverage_index), 
            "⚠️  Run preprocess_data_hybrid.R for coverage data",
            "✅ All data components loaded successfully"),
      sep = "\n"
    )
  })
  
  # ==============================================================================
  # PRIMER ALIGNMENT TAB SERVER LOGIC
  # ==============================================================================
  
  # Check if alignment data is available
  output$alignment_available <- reactive({
    nrow(alignment_summary) > 0
  })
  outputOptions(output, "alignment_available", suspendWhenHidden = FALSE)
  
  # Update primer index filter choices based on selected gene
  observe({
    if (nrow(alignment_summary) > 0 && !is.null(input$gene_id)) {
      gene_data <- alignment_summary[alignment_summary$gene_id == input$gene_id, ]
      if (nrow(gene_data) > 0) {
        primer_indices <- sort(unique(gene_data$primer_index))
        primer_choices <- c("All Primers" = "all", setNames(primer_indices, paste("Primer", primer_indices)))
        updateSelectInput(session, "alignment_primer_filter", 
                         choices = primer_choices, 
                         selected = "all")
      } else {
        updateSelectInput(session, "alignment_primer_filter", 
                         choices = c("No primers for this gene" = "none"), 
                         selected = "none")
      }
    }
  })
  
  # Filtered alignment data
  filtered_alignment <- reactive({
    if (nrow(alignment_summary) == 0 || is.null(input$gene_id)) return(data.frame())
    
    # Filter by selected gene from main panel
    data <- alignment_summary[alignment_summary$gene_id == input$gene_id, ]
    
    if (nrow(data) == 0) return(data.frame())
    
    # Filter by primer index
    if (!is.null(input$alignment_primer_filter) && input$alignment_primer_filter != "all" && input$alignment_primer_filter != "none") {
      data <- data[data$primer_index == as.numeric(input$alignment_primer_filter), ]
    }
    
    # Filter for perfect matches only
    if (!is.null(input$alignment_perfect_only) && input$alignment_perfect_only) {
      data <- data[data$mismatches == 0, ]
    }
    
    return(data)
  })

  # Download handler for plot PNG
  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste0("gene_plot_", input$gene_id, ".png")
    },
    content = function(file) {
      # Use the same plotting function and arguments as renderPlot
      plot_obj <- plot_gene_comprehensive(
        gene_id = input$gene_id,
        qc_data = qc_data,
        all_peaks_data = all_peaks_data,
        yaxis_mode = input$yaxis_mode,
        show_primers = input$show_primers,
        max_isoforms = input$max_isoforms,
        show_peak_region = input$show_peak_region,
        show_all_peaks = input$show_all_peaks,
        max_peaks = if(!is.null(input$max_peaks)) input$max_peaks else 3
      )
      ggsave(file, plot = plot_obj, width = 12, height = 8, dpi = 100)
    }
  )

  # Download handler for alignment summary CSV
  output$download_alignment_csv <- downloadHandler(
    filename = function() {
      paste0("alignment_summary_", input$gene_id, ".csv")
    },
    content = function(file) {
      data <- filtered_alignment()
      if (nrow(data) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
      } else {
        write.csv(data, file, row.names = FALSE)
      }
    }
  )
  
  # Alignment statistics
  output$alignment_stats <- renderText({
    if (nrow(alignment_summary) == 0) return("No alignment data available")
    if (is.null(input$gene_id)) return("Please select a gene")
    
    filtered_data <- filtered_alignment()
    
    if (nrow(filtered_data) == 0) {
      return(paste("No alignment data available for gene:", input$gene_id))
    }
    
    # Calculate statistics for selected gene
    total_alignments <- nrow(filtered_data)
    unique_primers <- length(unique(filtered_data$primer_index))
    perfect_matches <- sum(filtered_data$mismatches == 0, na.rm = TRUE)
    unique_transcripts <- length(unique(filtered_data$aligned_transcript))
    
    # Mismatch distribution
    mismatch_table <- table(filtered_data$mismatches)
    mismatch_summary <- paste(names(mismatch_table), "mismatches:", mismatch_table, collapse = " | ")
    
    paste(
      paste("Gene:", input$gene_id),
      paste("Total alignments:", total_alignments),
      paste("Primers analyzed:", unique_primers),
      paste("Aligned transcripts:", unique_transcripts),
      paste("Perfect matches:", perfect_matches, paste0("(", round(100 * perfect_matches / total_alignments, 1), "%)")),
      paste("Mismatch distribution:", mismatch_summary),
      sep = "\n"
    )
  })
  
  # Alignment data table
  output$alignment_table <- DT::renderDataTable({
    filtered_data <- filtered_alignment()
    
    if (nrow(filtered_data) == 0) {
      return(data.frame(Message = "No data available or no matches for current filters"))
    }
    
    # Select and rename columns for display
    display_data <- filtered_data[, c("gene_id", "primer_index", "primer_type", 
                                     "primer_sequence", "aligned_transcript", 
                                     "aligned_gene_name", "alignment_start", 
                                     "alignment_end", "alignment_length", 
                                     "mismatches", "distance_to_end")]
    
    colnames(display_data) <- c("Gene ID", "Primer Index", "Type", "Sequence", 
                               "Transcript", "Gene Name", "Start", "End", 
                               "Length", "Mismatches", "Distance to End")
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = c(1, 2, 6, 7, 8, 9, 10)),
          list(width = "120px", targets = 3),  # Sequence column
          list(width = "80px", targets = c(1, 2, 8, 9, 10))  # Numeric columns
        )
      ),
      class = "cell-border stripe",
      rownames = FALSE
    ) %>%
      DT::formatStyle("Mismatches",
        backgroundColor = DT::styleInterval(c(0, 1), c("#d4edda", "#fff3cd", "#f8d7da"))
      )
  })
}

# ==============================================================================
# LAUNCH APP
# ==============================================================================

# Restore original working directory on exit
on.exit(setwd(original_wd))

cat("\n=== PeakPrime Standalone Explorer ===\n")
cat("🚀 Ultra-fast loading + ⭐ Full features\n")
cat("📁 Source:", results_dir, "\n")
cat("Ready to explore", nrow(qc_data), "genes!\n\n")

# Logo already configured in UI setup section

# Create the Shiny app object
app <- shinyApp(ui = ui, server = server)

# Always make app available in global environment
assign("peakprime_app", app, envir = globalenv())

# Check if we should auto-launch
auto_launch <- exists(".standalone_launch", envir = globalenv()) && 
               get(".standalone_launch", envir = globalenv()) == TRUE

if (auto_launch || !interactive()) {
  cat("🚀 Launching Shiny app...\n")
  cat("💡 The app will open in your browser or RStudio Viewer\n")
  cat("🛑 Press Ctrl+C (or Escape in RStudio) to stop the app\n\n")
  
  # Run the app
  runApp(app, launch.browser = TRUE)
} else {
  cat("📝 App loaded successfully! Use one of these commands to launch:\n")
  cat("   runApp(peakprime_app)                    # Launch immediately\n")
  cat("   .standalone_launch <- TRUE; source(...)  # Auto-launch on next source\n\n")
}