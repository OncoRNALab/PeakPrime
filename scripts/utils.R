# Filter genes with cutoff of 10 counts in both replicates
# A gene passes the filter if it has >= 10 counts in BOTH replicates 
# for at least ONE condition

filter_genes_by_replicates <- function(count_matrix, annotation, count_cutoff = 10) {
  # Ensure annotation is ordered to match count_matrix columns
  annotation <- annotation %>%
    filter(SampleID %in% colnames(count_matrix)) %>%
    arrange(match(SampleID, colnames(count_matrix)))
  
  # Get unique conditions
  conditions <- unique(annotation$Condition)
  
  # For each gene, check if it passes the filter in at least one condition
  genes_to_keep <- apply(count_matrix, 1, function(gene_counts) {
    # Check each condition
    for (condition in conditions) {
      # Get sample indices for this condition
      sample_indices <- which(annotation$Condition == condition)
      
      # Get counts for this condition
      condition_counts <- gene_counts[sample_indices]
      
      # Check if ALL replicates have >= cutoff counts
      if (all(condition_counts >= count_cutoff)) {
        return(TRUE)  # Gene passes in this condition
      }
    }
    return(FALSE)  # Gene doesn't pass in any condition
  })
  
  # Filter the count matrix
  filtered_matrix <- count_matrix[genes_to_keep, ]
  
  # Return results
  list(
    filtered_matrix = filtered_matrix,
    n_genes_kept = sum(genes_to_keep),
    n_genes_total = nrow(count_matrix),
    proportion_kept = sum(genes_to_keep) / nrow(count_matrix)
  )
}


#plot rlog count correlation

plot_count_corr <- function(plot_data){
  ggplot(plot_data, aes(x = Targeted, y = Random)) +
    geom_point(alpha = 0.3, size = 1, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray40", linewidth = 0.7) +
    facet_wrap(~ facet_label, scales = "free", ncol = 2) +
    labs(
      title = "Correlation Targeted vs Random Priming",
      x = "Targeted mean rlog expression",
      y = "Random mean rlog expression"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}