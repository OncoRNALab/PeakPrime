#!/usr/bin/env Rscript

# PeakPrime Standalone Data Preprocessor
# Creates fast-loading binary data files for standalone app
#
# This version includes:
# - BigWig coverage data extraction for all genes
# - Complete GTF processing with transcript structures  
# - Optimized data structures for instant loading
# - Flexible GTF path parameter
#
# Usage: preprocess_peakprime_standalone("results/directory", "path/to/gtf")

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

#' Enhanced preprocessing for standalone app with flexible GTF path
#' @param results_dir Path to results directory containing PeakPrime outputs
#' @param gtf_path Path to GTF annotation file
#' @param max_coverage_points Maximum coverage points per gene for performance
#' @param coverage_window_pad Window padding for coverage extraction
preprocess_peakprime_standalone <- function(results_dir, 
                                           gtf_path = NULL,
                                           max_coverage_points = 1000,
                                           coverage_window_pad = 500) {
  
  # Validate and normalize results directory
  if (!dir.exists(results_dir)) {
    stop("❌ Results directory does not exist: ", results_dir)
  }
  
  results_dir <- normalizePath(results_dir, mustWork = TRUE)
  cat("🔍 Discovering PeakPrime files in:", results_dir, "\n")
  
  # Auto-detect GTF if not provided
  if (is.null(gtf_path)) {
    # Try common GTF locations
    gtf_candidates <- c(
      "/data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf",
      "~/annotations/Homo_sapiens.GRCh38.109.gtf",
      "../annotations/Homo_sapiens.GRCh38.109.gtf"
    )
    
    for (candidate in gtf_candidates) {
      if (file.exists(candidate)) {
        gtf_path <- candidate
        cat("📄 Auto-detected GTF:", gtf_path, "\n")
        break
      }
    }
    
    if (is.null(gtf_path)) {
      stop("❌ No GTF file found. Please provide gtf_path parameter.")
    }
  }
  
  if (!file.exists(gtf_path)) {
    stop("❌ GTF file does not exist: ", gtf_path)
  }
  
  # File discovery with robust pattern matching in results directory
  files <- list(
    qc_summary = list.files(results_dir, pattern = "peaks_qc_summary\\.tsv$", recursive = TRUE, full.names = TRUE)[1],
    selected_peaks = list.files(results_dir, pattern = "selected_peaks\\.tsv$", recursive = TRUE, full.names = TRUE)[1],
    bigwig = list.files(results_dir, pattern = "\\.bw$", recursive = TRUE, full.names = TRUE)[1],
    gtf = gtf_path,
    narrowpeak = list.files(results_dir, pattern = "_peaks\\.narrowPeak$", recursive = TRUE, full.names = TRUE)[1]
  )
  
  # Report file discovery status
  missing_files <- c()
  for (name in names(files)) {
    if (!is.na(files[[name]]) && file.exists(files[[name]])) {
      cat("✓ Found", name, ":", basename(files[[name]]), "\n")
    } else {
      cat("✗ Missing", name, "\n")
      files[[name]] <- NA
      missing_files <- c(missing_files, name)
    }
  }
  
  # Check for mandatory files
  mandatory_files <- c("qc_summary", "bigwig", "gtf")
  missing_mandatory <- intersect(missing_files, mandatory_files)
  
  if (length(missing_mandatory) > 0) {
    cat("\n❌ ERROR: Missing mandatory files:", paste(missing_mandatory, collapse = ", "), "\n")
    cat("GTF file is required for proper gene annotation visualization.\n")
    cat("Expected GTF: /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf\n")
    stop("Cannot proceed without mandatory input files.")
  }
  
  # ==== PROCESS QC DATA ====
  cat("\n📊 Processing QC summary data...\n")
  qc_raw <- fread(files$qc_summary)
  
  # Extract coordinates and enhance QC data - process row by row
  qc_enhanced <- qc_raw[, {
    # Process each row individually to preserve gene-specific coordinates
    gene_start_val <- pmax(1, gene_start - coverage_window_pad)
    gene_end_val <- gene_end + coverage_window_pad
    
    # Create peak region string for display
    peak_region_str <- ifelse(!is.na(peak_chr) & !is.na(peak_start) & !is.na(peak_end),
                             paste0(peak_chr, ":", peak_start, "-", peak_end),
                             NA_character_)
    
    # Return data preserving all original values per row
    list(
      gene_id = gene_id,
      final_selection = final_selection,
      best_peak_score = best_peak_score,
      best_peak_pvalue = best_peak_pvalue, 
      best_peak_qvalue = best_peak_qvalue,
      exonic_fraction = exonic_fraction,
      gene_chr = gene_chr,
      gene_start = gene_start_val,
      gene_end = gene_end_val,
      gene_strand = gene_strand,
      peak_chr = peak_chr,      # Use original column directly
      peak_start = peak_start,  # Use original column directly
      peak_end = peak_end,      # Use original column directly
      peak_region = peak_region_str,
      failure_reason = failure_reason
    )
  }, by = seq_len(nrow(qc_raw))]  # Process by row index to preserve individual values
  
  # Remove the grouping column
  qc_enhanced[, seq_len := NULL]
  
  cat("   Processed", nrow(qc_enhanced), "genes\n")
  
  # ==== PROCESS SELECTED PEAKS ====
  selected_data <- NULL
  if (!is.na(files$selected_peaks) && file.exists(files$selected_peaks)) {
    cat("🎯 Processing selected peaks data...\n")
    selected_data <- fread(files$selected_peaks)
    cat("   Loaded", nrow(selected_data), "selected peaks\n")
  }
  
  # ==== PROCESS GTF ANNOTATIONS ====
  gtf_processed <- NULL
  if (!is.na(files$gtf) && file.exists(files$gtf)) {
    cat("🧬 Processing GTF annotations...\n")
    
    tryCatch({
      gtf_full <- import(files$gtf)
      
      # Keep only essential feature types and columns for speed
      essential_types <- c("gene", "transcript", "exon", "CDS", "start_codon", "stop_codon", 
                          "five_prime_utr", "three_prime_utr", "UTR")
      
      gtf_filtered <- gtf_full[mcols(gtf_full)$type %in% essential_types]
      
      # Simplify metadata to essential columns
      essential_cols <- intersect(colnames(mcols(gtf_filtered)), 
                                 c("type", "gene_id", "gene_name", "transcript_id", 
                                   "exon_number", "protein_id"))
      
      if (length(essential_cols) > 0) {
        mcols(gtf_filtered) <- mcols(gtf_filtered)[essential_cols]
        gtf_processed <- gtf_filtered
        cat("   Processed", length(gtf_processed), "GTF features\n")
      } else {
        cat("   Warning: No essential GTF columns found\n")
      }
      
    }, error = function(e) {
      cat("   Error processing GTF:", e$message, "\n")
    })
  }
  
  # ==== PROCESS BIGWIG COVERAGE DATA ====
  coverage_index <- list()
  
  if (!is.na(files$bigwig) && file.exists(files$bigwig)) {
    cat("📈 Extracting BigWig coverage data for all genes...\n")
    cat("   This may take a few minutes for", nrow(qc_enhanced), "genes...\n")
    
    # Process genes with valid coordinates
    valid_genes <- qc_enhanced[!is.na(gene_chr) & !is.na(gene_start) & !is.na(gene_end)]
    
    tryCatch({
      for (i in 1:nrow(valid_genes)) {
        gene_info <- valid_genes[i]
        
        if (i %% 10 == 0) {
          cat("   Processing gene", i, "of", nrow(valid_genes), 
              paste0("(", round(100 * i / nrow(valid_genes)), "%)...\\r"))
        }
        
        # Create genomic range for this gene
        gene_gr <- GRanges(
          seqnames = gene_info$gene_chr,
          ranges = IRanges(start = gene_info$gene_start, end = gene_info$gene_end)
        )
        
        # Extract coverage
        cov_data <- tryCatch({
          import(files$bigwig, which = gene_gr)
        }, error = function(e) {
          GRanges()  # Return empty if extraction fails
        })
        
        if (length(cov_data) > 0) {
          # Convert to simplified data table
          cov_df <- data.table(
            position = as.integer((start(cov_data) + end(cov_data)) / 2),
            coverage = as.numeric(mcols(cov_data)$score)
          )
          
          # Remove invalid values
          cov_df <- cov_df[is.finite(coverage) & coverage >= 0]
          
          # Downsample if too many points (for performance)
          if (nrow(cov_df) > max_coverage_points) {
            keep_idx <- seq(1, nrow(cov_df), length.out = max_coverage_points)
            cov_df <- cov_df[keep_idx]
          }
          
          # Store in index
          if (nrow(cov_df) > 0) {
            coverage_index[[gene_info$gene_id]] <- cov_df
          }
        }
      }
      
      cat("\\n   Extracted coverage for", length(coverage_index), "genes\\n")
      
    }, error = function(e) {
      cat("\\n   Error during BigWig processing:", e$message, "\\n")
      cat("   Continuing without coverage data...\\n")
    })
  } else {
    cat("📈 No BigWig file found - coverage plots will be unavailable\\n")
  }
  
  # ==== SAVE PROCESSED DATA ====
  cat("\\n💾 Saving processed data files...\\n")
  
  # Save QC data to results directory
  saveRDS(qc_enhanced, file.path(results_dir, "qc_data.rds"))
  cat("   ✓ qc_data.rds\\n")
  
  # Save selected peaks
  if (!is.null(selected_data)) {
    saveRDS(selected_data, file.path(results_dir, "peaks_data.rds"))
    cat("   ✓ peaks_data.rds\\n")
  }
  
  # Save GTF data
  if (!is.null(gtf_processed)) {
    saveRDS(gtf_processed, file.path(results_dir, "gtf_data.rds"))
    cat("   ✓ gtf_data.rds\\n")
  }
  
  # Save coverage index
  if (length(coverage_index) > 0) {
    saveRDS(coverage_index, file.path(results_dir, "coverage_index.rds"))
    cat("   ✓ coverage_index.rds\\n")
  }
  
  # Create enhanced manifest
  manifest <- list(
    files = files,
    processed_time = Sys.time(),
    processing_version = "standalone_v1.0",
    results_dir = results_dir,
    gtf_path = gtf_path,
    qc_genes = nrow(qc_enhanced),
    selected_genes = sum(qc_enhanced$final_selection, na.rm = TRUE),
    coverage_genes = length(coverage_index),
    gtf_features = if(!is.null(gtf_processed)) length(gtf_processed) else 0,
    settings = list(
      max_coverage_points = max_coverage_points,
      coverage_window_pad = coverage_window_pad
    )
  )
  
  saveRDS(manifest, file.path(results_dir, "data_manifest.rds"))
  cat("   ✓ data_manifest.rds\\n")
  
  # ==== SUMMARY REPORT ====
  cat("\\n🎉 Enhanced preprocessing complete!\\n")
  cat("\\n=== SUMMARY ===\\n")
  cat("Total genes:", manifest$qc_genes, "\\n")
  cat("Selected genes:", manifest$selected_genes, "\\n")  
  cat("Genes with coverage:", manifest$coverage_genes, "\\n")
  cat("GTF features:", manifest$gtf_features, "\\n")
  cat("\\n=== FILES CREATED ===\\n")
  cat("📊 qc_data.rds - Enhanced QC data with coordinates\\n")
  if (!is.null(selected_data)) cat("🎯 peaks_data.rds - Selected peaks data\\n")
  if (!is.null(gtf_processed)) cat("🧬 gtf_data.rds - Processed gene annotations\\n") 
  if (length(coverage_index) > 0) cat("📈 coverage_index.rds - Preprocessed coverage data\\n")
  cat("📋 data_manifest.rds - Processing metadata\\n")
  cat("\\n🚀 Now launch app_standalone.R with this results directory for ultra-fast exploration!\\n")
  cat("   Usage: Rscript app_standalone.R", results_dir, "\\n")
  
  return(invisible(manifest))
}

# ==== COMMAND LINE USAGE ====
# This script can be sourced as a library or run directly with parameters
# Example usage:
#   source("preprocess_for_standalone.R")  
#   preprocess_peakprime_standalone("results/testboth", "/path/to/gtf")

cat("📦 Standalone preprocessing functions loaded.\\n")
cat("   Use: preprocess_peakprime_standalone(results_dir, gtf_path)\\n")