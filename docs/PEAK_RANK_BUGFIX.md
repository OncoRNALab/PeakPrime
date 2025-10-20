# Bug Fix: Peak Rank Column Row Mismatch

## Issue

When implementing the `peak_rank` column addition, a critical bug was introduced that caused pipeline failures when genes were filtered out due to exonic trimming.

### Error Message

```
Error in data.frame(gene = target_regions$gene_id, chr = seqnames(target_regions),  : 
  arguments imply differing number of rows: 46, 49, 1
Execution halted
```

### Root Cause

The error occurred because:

1. **Peak selection**: 49 genes selected (stored in `best_peaks` dataframe)
2. **Exonic trimming**: 3 genes failed forced exonic trimming and were removed from `target_regions`
3. **Output generation**: Tried to create dataframe using:
   - `target_regions` (46 rows) 
   - `best_peaks$actual_peak_rank` (49 rows) ← **Mismatch!**

The `target_regions` GRanges object was filtered, but the `best_peaks` dataframe was not kept in sync, causing a row count mismatch when creating the output dataframe.

## Example Case

```
Command: process_macs2_peaks.R --peak_rank 2 --force_exonic_trimming true --min_trimmed_length 30

Log output:
  Selected 49 best peaks (one per gene)
  Warning: Gene ENSG00000225663 trimmed region too short ( 19 bp <  30 bp)
  Removed 3 genes due to failed forced exonic trimming
  Extracting sequences...
  Error in data.frame(...)
```

### Affected Genes
- ENSG00000225663 (trimmed region: 19 bp < 30 bp minimum)
- 2 other genes (no exonic overlap or too short after trimming)

## Solution

Modified `bin/process_macs2_peaks.R` to keep `best_peaks` dataframe synchronized with `target_regions` filtering.

### Changes

#### 1. Forced Exonic Trimming Filter (line ~452)

**Before:**
```r
if(force_exonic_trimming) {
  failed_trimming <- target_regions$gene_id[!trimmed_to_exon & exonic_fractions == 0]
  if(length(failed_trimming) > 0) {
    keep_idx <- !(target_regions$gene_id %in% failed_trimming)
    genes_before_trimming <- target_regions$gene_id
    target_regions <- target_regions[keep_idx]
    # ❌ best_peaks not filtered!
    
    all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% failed_trimming] <- FALSE
    all_genes_qc$failure_reason[all_genes_qc$gene_id %in% failed_trimming] <- "Failed forced exonic trimming (no exonic overlap or trimmed region too short)"
    
    cat("Removed", length(failed_trimming), "genes due to failed forced exonic trimming\n")
  }
}
```

**After:**
```r
if(force_exonic_trimming) {
  failed_trimming <- target_regions$gene_id[!trimmed_to_exon & exonic_fractions == 0]
  if(length(failed_trimming) > 0) {
    keep_idx <- !(target_regions$gene_id %in% failed_trimming)
    genes_before_trimming <- target_regions$gene_id
    target_regions <- target_regions[keep_idx]
    
    # ✅ Also filter best_peaks to keep them in sync
    best_peaks <- best_peaks[!(best_peaks$gene_id %in% failed_trimming), ]
    
    all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% failed_trimming] <- FALSE
    all_genes_qc$failure_reason[all_genes_qc$gene_id %in% failed_trimming] <- "Failed forced exonic trimming (no exonic overlap or trimmed region too short)"
    
    cat("Removed", length(failed_trimming), "genes due to failed forced exonic trimming\n")
  }
}
```

#### 2. Exonic Fraction Filter (line ~466)

**Before:**
```r
if(!is.na(min_exonic_fraction) && !force_exonic_trimming) {
  keep_idx <- exonic_fractions >= min_exonic_fraction
  genes_before_filter <- target_regions$gene_id
  target_regions <- target_regions[keep_idx]
  genes_after_filter <- target_regions$gene_id
  # ❌ best_peaks not filtered!
  
  all_genes_qc$passes_exonic_filter <- FALSE
  all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% genes_after_filter] <- TRUE
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_after_filter] <- "Passed all filters"
  
  failed_exonic <- setdiff(genes_before_filter, genes_after_filter)
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% failed_exonic] <- paste0("Failed exonic fraction filter (<", min_exonic_fraction, ")")
  
  cat("Retained", length(target_regions), "regions after exonic fraction filtering\n")
}
```

**After:**
```r
if(!is.na(min_exonic_fraction) && !force_exonic_trimming) {
  keep_idx <- exonic_fractions >= min_exonic_fraction
  genes_before_filter <- target_regions$gene_id
  target_regions <- target_regions[keep_idx]
  genes_after_filter <- target_regions$gene_id
  
  # ✅ Also filter best_peaks to keep them in sync
  best_peaks <- best_peaks[best_peaks$gene_id %in% genes_after_filter, ]
  
  all_genes_qc$passes_exonic_filter <- FALSE
  all_genes_qc$passes_exonic_filter[all_genes_qc$gene_id %in% genes_after_filter] <- TRUE
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% genes_after_filter] <- "Passed all filters"
  
  failed_exonic <- setdiff(genes_before_filter, genes_after_filter)
  all_genes_qc$failure_reason[all_genes_qc$gene_id %in% failed_exonic] <- paste0("Failed exonic fraction filter (<", min_exonic_fraction, ")")
  
  cat("Retained", length(target_regions), "regions after exonic fraction filtering\n")
}
```

## Testing

### Before Fix
```bash
nextflow run main.nf \
  --peak_rank 2 \
  --force_exonic_trimming \
  --min_trimmed_length 30 \
  -resume

# Result: Error exit status (1)
# Error: arguments imply differing number of rows: 46, 49, 1
```

### After Fix
```bash
nextflow run main.nf \
  --peak_rank 2 \
  --force_exonic_trimming \
  --min_trimmed_length 30 \
  -resume

# Result: Success ✓
# Output: selected_peaks.tsv with 46 genes (3 removed due to failed trimming)
```

### Verification

Check that row counts match:

```bash
# Count rows in output (excluding header)
tail -n+2 selected_peaks.tsv | wc -l
# Should match: 46

# Verify peak_rank column exists and has valid values
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="peak_rank") print "Column",i":",$i} NR>1 {print $12}' selected_peaks.tsv | head
# Should show: Column 12: peak_rank
#              1
#              2
#              1  (for genes that fell back to rank 1)
```

## Impact

### Affected Scenarios

This bug affected pipelines when:
1. ✅ `--force_exonic_trimming` is enabled AND genes fail trimming
2. ✅ `--min_exonic_fraction` is set AND genes fail the threshold

### Not Affected

- ❌ No filtering (all selected peaks kept)
- ❌ All genes pass filters

## Prevention

The fix ensures that `best_peaks` is always kept synchronized with `target_regions` by:
1. Applying the same filtering conditions to both data structures
2. Filtering `best_peaks` immediately after filtering `target_regions`
3. Using `gene_id` as the common key for filtering

This maintains the invariant: **`length(target_regions) == nrow(best_peaks)`** throughout the script.

## Related Files

- `bin/process_macs2_peaks.R`: Core logic fix
- `docs/PEAK_RANK_COLUMN_IMPLEMENTATION.md`: Original feature documentation
- `docs/PEAK_RANK_BUGFIX.md`: This document

## Lesson Learned

When adding a new column that depends on synchronized data structures:
1. Identify all places where filtering occurs
2. Ensure ALL related data structures are filtered consistently
3. Test with scenarios that trigger filtering (not just the happy path)
4. Add assertions to verify invariants (e.g., row count matching)
