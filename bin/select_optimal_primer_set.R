#!/usr/bin/env Rscript
# Select optimal primer combination that maximizes isoform coverage
#
# This script implements a greedy algorithm to select primers that:
# 1. Map perfectly (0 mismatches) to only one gene (target specificity)
# 2. Are within a distance threshold from the transcript 3' end
# 3. Maximize the number of distinct isoforms covered

library(data.table)

#' Filter primers by specificity (no off-target perfect matches)
filter_specific_primers <- function(dt, target_gene) {
  cat(sprintf("\nFiltering primers for gene: %s\n", target_gene))
  
  # Get all alignments for this gene's primers
  gene_primers <- dt[gene_id == target_gene]
  cat(sprintf("  Total alignments for %s: %d\n", target_gene, nrow(gene_primers)))
  
  primer_indices <- unique(gene_primers$primer_index)
  cat(sprintf("  Unique primers: %d\n", length(primer_indices)))
  
  # Check each primer for specificity
  specific_primers <- sapply(primer_indices, function(idx) {
    # Get all alignments for this primer
    primer_aligns <- gene_primers[primer_index == idx]
    
    # Check perfect matches
    perfect <- primer_aligns[mismatches == 0]
    
    # Get unique genes with perfect matches
    perfect_genes <- unique(perfect$aligned_gene_name)
    
    # Keep if all perfect matches are to the same gene (target gene)
    return(length(perfect_genes) == 1)
  })
  
  specific_indices <- primer_indices[specific_primers]
  cat(sprintf("  Specific primers (no off-target): %d\n", length(specific_indices)))
  
  return(gene_primers[primer_index %in% specific_indices])
}

#' Filter primers by distance to transcript end
filter_by_distance <- function(dt, max_distance) {
  cat(sprintf("\nFiltering by distance to end (max=%dnt)\n", max_distance))
  
  before <- length(unique(dt$primer_index))
  result <- dt[distance_to_end <= max_distance]
  after <- length(unique(result$primer_index))
  
  cat(sprintf("  Primers before: %d, after: %d, removed: %d\n", 
              before, after, before - after))
  
  return(result)
}

#' Greedy algorithm to select optimal primer set
select_optimal_primers <- function(dt, max_primers = 5) {
  cat(sprintf("\nSelecting optimal primer combination (max=%d primers)\n", max_primers))
  
  # Build primer -> isoforms mapping (only perfect matches)
  perfect_matches <- dt[mismatches == 0]
  
  primer_to_isoforms <- split(
    perfect_matches$aligned_transcript,
    perfect_matches$primer_index
  )
  primer_to_isoforms <- lapply(primer_to_isoforms, unique)
  
  # Get primer metadata
  primer_info <- perfect_matches[, .(
    sequence = first(primer_sequence),
    gene_id = first(gene_id),
    gene_name = first(aligned_gene_name),
    gene_strand = first(gene_strand),
    min_distance = min(distance_to_end)
  ), by = primer_index]
  
  cat(sprintf("  Total candidate primers: %d\n", length(primer_to_isoforms)))
  
  all_isoforms <- unique(unlist(primer_to_isoforms))
  cat(sprintf("  Total isoforms to cover: %d\n", length(all_isoforms)))
  
  # Greedy selection
  selected_primers <- integer(0)
  covered_isoforms <- character(0)
  
  for (iteration in seq_len(max_primers)) {
    if (length(primer_to_isoforms) == 0) {
      cat(sprintf("  Iteration %d: No more candidate primers\n", iteration))
      break
    }
    
    # Calculate marginal coverage for each remaining primer
    marginal_coverage <- sapply(primer_to_isoforms, function(isoforms) {
      new_isoforms <- setdiff(isoforms, covered_isoforms)
      length(new_isoforms)
    })
    
    # Find best primer
    best_idx <- which.max(marginal_coverage)
    best_new_coverage <- marginal_coverage[best_idx]
    
    if (best_new_coverage == 0) {
      cat(sprintf("  Iteration %d: No additional isoforms can be covered\n", iteration))
      break
    }
    
    # Get primer index and new isoforms
    best_primer <- as.integer(names(primer_to_isoforms)[best_idx])
    new_isoforms <- setdiff(primer_to_isoforms[[best_idx]], covered_isoforms)
    
    # Update selection
    selected_primers <- c(selected_primers, best_primer)
    covered_isoforms <- c(covered_isoforms, new_isoforms)
    primer_to_isoforms[[best_idx]] <- NULL
    
    coverage_pct <- (length(covered_isoforms) / length(all_isoforms)) * 100
    
    cat(sprintf("  Iteration %d: Selected primer %d\n", iteration, best_primer))
    cat(sprintf("    - Covers %d new isoforms: %s\n", 
                best_new_coverage, paste(sort(new_isoforms), collapse=", ")))
    cat(sprintf("    - Total coverage: %d/%d (%.1f%%)\n", 
                length(covered_isoforms), length(all_isoforms), coverage_pct))
    
    min_dist <- primer_info[primer_index == best_primer, min_distance]
    cat(sprintf("    - Distance to end: %dnt\n", min_dist))
  }
  
  # Prepare results
  if (length(selected_primers) == 0) {
    return(NULL)
  }
  
  results <- lapply(seq_along(selected_primers), function(i) {
    primer_idx <- selected_primers[i]
    info <- primer_info[primer_index == primer_idx]
    
    # Get all isoforms covered by this primer
    primer_isoforms <- perfect_matches[primer_index == primer_idx, 
                                       unique(aligned_transcript)]
    
    data.table(
      rank = i,
      primer_index = primer_idx,
      primer_sequence = info$sequence,
      gene_id = info$gene_id,
      gene_name = info$gene_name,
      gene_strand = info$gene_strand,
      min_distance_to_end = info$min_distance,
      isoforms_covered = length(primer_isoforms),
      isoform_ids = paste(sort(primer_isoforms), collapse = ",")
    )
  })
  
  results_dt <- rbindlist(results)
  
  return(list(
    results = results_dt,
    covered_isoforms = covered_isoforms,
    all_isoforms = all_isoforms
  ))
}

#' Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript select_optimal_primer_set.R <input.tsv> <gene_id> [max_distance] [max_primers] [output.tsv]\n")
    cat("\nExample:\n")
    cat("  Rscript select_optimal_primer_set.R primer_alignment_summary.tsv ENSG00000104687 400 5 optimal_primers.tsv\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  gene_id <- args[2]
  max_distance <- if (length(args) >= 3) as.integer(args[3]) else 400
  max_primers <- if (length(args) >= 4) as.integer(args[4]) else 5
  output_file <- if (length(args) >= 5) args[5] else "optimal_primers.tsv"
  
  cat(sprintf("Parameters:\n"))
  cat(sprintf("  Input file: %s\n", input_file))
  cat(sprintf("  Gene ID: %s\n", gene_id))
  cat(sprintf("  Max distance: %d nt\n", max_distance))
  cat(sprintf("  Max primers: %d\n", max_primers))
  cat(sprintf("  Output file: %s\n", output_file))
  
  # Load data
  cat(sprintf("\nLoading alignment data from: %s\n", input_file))
  dt <- fread(input_file)
  cat(sprintf("  Loaded %d alignment records\n", nrow(dt)))
  cat(sprintf("  Unique genes: %d\n", length(unique(dt$gene_id))))
  
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("Processing gene: %s\n", gene_id))
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  
  # Filter by specificity
  specific <- filter_specific_primers(dt, gene_id)
  
  if (nrow(specific) == 0) {
    cat(sprintf("ERROR: No specific primers found for %s\n", gene_id))
    quit(status = 1)
  }
  
  # Filter by distance
  filtered <- filter_by_distance(specific, max_distance)
  
  if (nrow(filtered) == 0) {
    cat(sprintf("ERROR: No primers within distance threshold for %s\n", gene_id))
    quit(status = 1)
  }
  
  # Select optimal set
  result <- select_optimal_primers(filtered, max_primers)
  
  if (is.null(result)) {
    cat(sprintf("ERROR: No primers selected for %s\n", gene_id))
    quit(status = 1)
  }
  
  # Print summary
  coverage_pct <- (length(result$covered_isoforms) / length(result$all_isoforms)) * 100
  uncovered <- setdiff(result$all_isoforms, result$covered_isoforms)
  
  cat(sprintf("\nFINAL RESULTS for %s:\n", gene_id))
  cat(sprintf("  Selected %d primers\n", nrow(result$results)))
  cat(sprintf("  Coverage: %d/%d isoforms (%.1f%%)\n", 
              length(result$covered_isoforms), 
              length(result$all_isoforms), 
              coverage_pct))
  
  if (length(uncovered) > 0) {
    cat(sprintf("  Uncovered isoforms: %s\n", paste(sort(uncovered), collapse = ", ")))
  } else {
    cat("  All isoforms covered!\n")
  }
  
  # Write results
  fwrite(result$results, output_file, sep = "\t")
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("Results written to: %s\n", output_file))
}

# Run main if called from command line
if (!interactive()) {
  main()
}
