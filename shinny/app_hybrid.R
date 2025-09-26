#!/usr/bin/env Rscript

# PeakPrime Hybrid Explorer - Ultra-fast loading with full plotting features
# Combines the speed of ultra-fast mode with comprehensive plot_gene_with_window functionality
#
# Usage:
#   1. Run preprocess_data_hybrid.R once to create enhanced preprocessed files
#   2. Launch this app with: shiny::runApp("app_hybrid.R")

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
# DATA LOADING AND INITIALIZATION
# ==============================================================================

# Check for preprocessed data
if (!file.exists("data_manifest.rds")) {
  stop("❌ No preprocessed data found!\n",
       "Run preprocess_data_hybrid.R first to generate required data files.\n",
       "This will include GTF annotations which are mandatory for proper visualization.")
}

# Load all preprocessed data instantly
cat("🚀 Loading preprocessed data...\n")
start_time <- Sys.time()

manifest <- readRDS("data_manifest.rds")
qc_data_raw <- readRDS("qc_data.rds")
selected_peaks <- if (file.exists("peaks_data.rds")) readRDS("peaks_data.rds") else NULL
gtf_data <- if (file.exists("gtf_data.rds")) readRDS("gtf_data.rds") else NULL
coverage_index <- if (file.exists("coverage_index.rds")) readRDS("coverage_index.rds") else NULL

# Load all peaks data if available
all_peaks_data <- NULL
if (file.exists("macs2_peaks/Merged_S7_S12_peaks.narrowPeak")) {
  cat("📈 Loading all detected peaks data...\n")
  all_peaks_raw <- fread("macs2_peaks/Merged_S7_S12_peaks.narrowPeak", 
                        header = FALSE, 
                        col.names = c("chr", "start", "end", "name", "score", 
                                     "strand", "fold_change", "pvalue", "qvalue", "summit_offset"))
  
  # Process peaks into gene-indexed structure
  if (nrow(all_peaks_raw) > 0) {
    # Load GTF to map peaks to genes
    if (file.exists("processed_peaks/peaks_qc_summary.tsv")) {
      qc_summary <- fread("processed_peaks/peaks_qc_summary.tsv")
      
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
}

# Convert data.table to regular data.frame to avoid filtering issues
qc_data <- as.data.frame(qc_data_raw)
cat("📋 Converted data.table to data.frame for reliable filtering\n")

load_time <- Sys.time() - start_time
cat("⚡ Loaded", nrow(qc_data), "genes in", format(load_time, digits = 2), "\n")

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
  
  cat("Found", length(transcripts), "transcripts for gene", gene_id, "(including non-coding)\n")
  
  list(exons = exons, cds = cds, utrs = utrs, transcripts = transcripts)
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
plot_gene_comprehensive <- function(gene_id, qc_data, yaxis_mode = "depth", 
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
  cat("Selected/Trimmed peak coords:", peak_start_coord, "-", peak_end_coord, "\n")
  cat("Gene coords:", gene_start_coord, "-", gene_end_coord, "\n")
  
  # Note: peak_start_coord and peak_end_coord represent the TRIMMED/SELECTED peak
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
  
  # CRITICAL FIX: Expand plot window to always include peak coordinates
  if (!is.na(peak_start_coord) && !is.na(peak_end_coord)) {
    plot_start <- min(plot_start, peak_start_coord)
    plot_end <- max(plot_end, peak_end_coord)
    cat("Expanded plot window to include peak: ", plot_start, "-", plot_end, "\\n")
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
  if (show_peak_region && !is.na(peak_start_coord) && !is.na(peak_end_coord)) {
    p_cov <- p_cov + 
      annotate("rect", 
              xmin = peak_start_coord, 
              xmax = peak_end_coord,
              ymin = -Inf, ymax = Inf, 
              alpha = 0.3, fill = "yellow", color = "orange", linewidth = 0.8) +
      annotate("text", x = (peak_start_coord + peak_end_coord) / 2, y = Inf, 
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
  if (show_peak_region && !is.na(peak_start_coord) && !is.na(peak_end_coord)) {
    p_feat <- p_feat + 
      annotate("rect",
              xmin = peak_start_coord,
              xmax = peak_end_coord, 
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
  
  # Add primer annotations - ALWAYS draw if coordinates exist
  peak_start_val <- peak_start_coord
  peak_end_val <- peak_end_coord
  
  cat("Checking primers: show_primers =", show_primers, "peak_start =", peak_start_val, "peak_end =", peak_end_val, "\n")
  
  if (show_primers && !is.na(peak_start_val) && !is.na(peak_end_val)) {
    # Calculate primer position at SELECTED peak center (always use selected peak)
    selected_peak_center <- (peak_start_val + peak_end_val) / 2
    
    # Determine primer direction based on gene strand (towards 3' end)
    primer_direction <- ifelse(gene_strand == "-", "left", "right")  # - strand: 3' is left, + strand: 3' is right
    
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
        
        # Check if this matches the selected peak (exact coordinate match)
        is_selected_peak <- (!is.na(peak_start_val) && !is.na(peak_end_val) &&
                            abs(peak_start_coord - peak_start_val) < 50 && 
                            abs(peak_end_coord - peak_end_val) < 50)
        
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
    
    # If selected peak wasn't found in detected peaks, draw it separately
    if (!selected_peak_drawn && !is.na(peak_start_val) && !is.na(peak_end_val)) {
      cat("Drawing selected peak separately (not in detected peaks list)\n")
      p_feat <- p_feat +
        annotate("rect", 
                xmin = peak_start_val, xmax = peak_end_val,
                ymin = peak_track - 0.4, ymax = peak_track + 0.4,
                fill = "darkred", alpha = 0.9, 
                color = "black", linewidth = 0.5)
    }
  } else {
    # Draw only selected peak (original behavior)
    if (!is.na(peak_start_val) && !is.na(peak_end_val)) {
      cat("Drawing selected peak region from", peak_start_val, "to", peak_end_val, "on track", peak_track, "\n")
      
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
# SHINY UI
# ==============================================================================

ui <- fluidPage(
  titlePanel(
    div(style = "display: flex; align-items: center; margin-bottom: 20px;",
        img(src = "logo2.png", height = "60px", style = "margin-right: 15px;"),
        span("PeakPrime Explorer", style = "font-size: 28px; font-weight: bold;")
    )
  ),
  
  tags$head(
    tags$style(HTML("
      .sidebar { background-color: #f8f9fa; padding: 15px; }
      .performance-badge { 
        background: linear-gradient(45deg, #28a745, #20c997);
        color: white; padding: 5px 10px; border-radius: 15px;
        font-size: 0.8em; font-weight: bold;
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
      br(), br(),
      
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
                  value = 15, min = 1, max = 25, step = 1),
      
      br(),
      div(class = "feature-badge", "Full plot_gene_with_window features")
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        tabPanel("Gene Plot",
                plotOutput("comprehensive_plot", height = "650px")),
        
        tabPanel("Gene Data",
                h4("Gene Quality Control Summary"),
                DT::dataTableOutput("gene_table")),
        
        tabPanel("Performance Info",
                h4("Performance & Data Summary"),
                verbatimTextOutput("performance_info"),
                
                h4("Available Data"),
                verbatimTextOutput("data_summary"))
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
      "✓ Full plot_gene_with_window functionality",
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
}

# ==============================================================================
# APP LAUNCH
# ==============================================================================

# Run the application
if (interactive()) {
  cat("\n=== PeakPrime Hybrid Explorer ===\n")
  cat("🚀 Ultra-fast loading + ⭐ Full features\n")
  cat("Ready to explore", nrow(qc_data), "genes!\n\n")
  
  shinyApp(ui = ui, server = server)
} else {
  shinyApp(ui = ui, server = server)
}