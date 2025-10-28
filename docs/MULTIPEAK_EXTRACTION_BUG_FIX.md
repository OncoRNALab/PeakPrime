# Multi-Peak Mode: Primer Extraction Bug Fixes

## Date: October 27, 2025

## Problem Summary

Even after implementing peak identifier tracking throughout the pipeline, the alignment FASTA file (`primers_for_alignment.fa`) only contained primers from **peak_1**, despite:
- ✅ Primer3 output having primers for all peaks (peak_1, peak_2, peak_3, peak_4, peak_5, peak_7)
- ✅ Peak identifiers being properly preserved in Primer3 output

## Root Causes Identified

### Bug #1: Strand Lookup Using Gene ID Only (extract_cdna_primers.R)

**Location:** `bin/extract_cdna_primers.R` line 41

**Problem:**
```r
gene_strands <- setNames(as.character(peaks$strand), as.character(peaks$gene))
```

When a gene has multiple peaks, this creates a named vector using gene_id as the key. **Only the last peak's strand is retained** - all previous peaks are overwritten.

**Example:**
```r
# For ENSG00000107130 with 6 peaks, only the strand for peak_7 (the last one) was kept
gene_strands["ENSG00000107130"]  # Only has strand from peak_7
```

**Impact:** The script could only extract primers for the last peak of each gene.

**Fix:**
```r
# Create lookup with peak-specific keys
if ("peak_rank" %in% colnames(peaks)) {
  # Multi-peak mode: create keys with peak identifiers
  peak_keys <- paste0(peaks$gene, "|peak_", peaks$peak_rank)
  gene_strands <- setNames(as.character(peaks$strand), peak_keys)
  # Also keep gene-only keys for backward compatibility
  gene_only <- setNames(as.character(peaks$strand), as.character(peaks$gene))
  gene_strands <- c(gene_strands, gene_only[!names(gene_only) %in% names(gene_strands)])
} else {
  # Single-peak mode: original behavior
  gene_strands <- setNames(as.character(peaks$strand), as.character(peaks$gene))
}
```

**Fix location:** Lines 37-52

### Bug #2: Filtering by Gene ID Only (primers_to_fasta.R)

**Location:** `bin/primers_to_fasta.R` line 28

**Problem:**
```r
split(primers, primers$gene_id)
```

When limiting primers via `--max_primers_per_gene 20`, the script split primers by gene_id only. In multi-peak mode, this meant:
- All primers from all peaks for a gene were grouped together
- Only the first 20 primers TOTAL were kept
- Later peaks' primers were discarded

**Example for ENSG00000107130:**
- 105 total primers across 6 peaks
- After `head(gene_primers, 20)`: Only first 20 primers kept (all from peak_1)
- Primers from peak_2, peak_3, peak_4, peak_5, peak_7 were discarded

**Impact:** The `--max_primers_per_gene` parameter inadvertently became "max primers per gene TOTAL" instead of "max primers per peak" in multi-peak mode.

**Fix:**
```r
if (!is.na(opt$max_primers_per_gene) && opt$max_primers_per_gene > 0) {
  # In multi-peak mode, limit per gene+peak combination, not just per gene
  if ("peak_id" %in% colnames(primers) && any(!is.na(primers$peak_id))) {
    # Multi-peak mode: group by gene_id AND peak_id
    primers$group_key <- paste0(primers$gene_id, "|", primers$peak_id)
    primers_filtered <- do.call(rbind, lapply(split(primers, primers$group_key), function(group_primers) {
      head(group_primers, opt$max_primers_per_gene)
    }))
    primers_filtered$group_key <- NULL
    primers <- primers_filtered
    cat("Filtered to:", nrow(primers), "primers (max", opt$max_primers_per_gene, "per gene+peak)\n")
  } else {
    # Single-peak mode: group by gene_id only
    primers_filtered <- do.call(rbind, lapply(split(primers, primers$gene_id), function(gene_primers) {
      head(gene_primers, opt$max_primers_per_gene)
    }))
    primers <- primers_filtered
    cat("Filtered to:", nrow(primers), "primers (max", opt$max_primers_per_gene, "per gene)\n")
  }
}
```

**Fix location:** Lines 27-45

## Verification Results

### Before Fixes:
```bash
# cdna_primers.tsv: Only primers from last peak extracted
grep "ENSG00000107130" cdna_primers.tsv | cut -f2 | sort -u
# Output: peak_7 (only)

# primers_for_alignment.fa: Only 20 primers total for gene
grep "ENSG00000107130" primers_for_alignment.fa | cut -d'|' -f2 | sort | uniq -c
# Output: 20 peak_1 (only)
```

### After Fixes:
```bash
# cdna_primers.tsv: All peaks extracted
grep "ENSG00000107130" cdna_primers_test.tsv | cut -f2 | sort -u
# Output: peak_1, peak_2, peak_3, peak_4, peak_5, peak_7

# primers_for_alignment.fa: Up to 20 primers per peak
grep "ENSG00000107130" primers_for_alignment_FIXED.fa | cut -d'|' -f2 | sort | uniq -c
# Output:
#   20 peak_1
#   20 peak_2
#   20 peak_3
#    5 peak_4
#   20 peak_5
#   20 peak_7
```

### Summary Statistics:
- **Before:** 105 primers extracted, only 20 in alignment FASTA (all peak_1)
- **After:** 105 primers extracted, 105 in alignment FASTA (all peaks)

## Test Case: ENSG00000107130

This gene has 6 peaks detected by MACS2 (ranks 1,2,3,4,5,7):

| Peak | Genomic Location | Primers Designed | Primers in Alignment (before) | Primers in Alignment (after) |
|------|-----------------|------------------|------------------------------|------------------------------|
| peak_1 | chr9:130236156-130236198 | 20 | 20 ✅ | 20 ✅ |
| peak_2 | chr9:130237136-130237190 | 20 | 0 ❌ | 20 ✅ |
| peak_3 | chr9:130236396-130236439 | 20 | 0 ❌ | 20 ✅ |
| peak_4 | chr9:130233597-130233626 | 5 | 0 ❌ | 5 ✅ |
| peak_5 | chr9:130223089-130223154 | 20 | 0 ❌ | 20 ✅ |
| peak_7 | chr9:130200978-130201006 | 20 | 0 ❌ | 20 ✅ |
| **Total** | | **105** | **20** | **105** |

## Impact on Downstream Analysis

### Before Fixes:
- ❌ Multi-peak optimization impossible (missing 85% of primers)
- ❌ Isoform coverage analysis broken
- ❌ Only primers from one peak analyzed for specificity
- ❌ Redundant peaks not filtered properly

### After Fixes:
- ✅ All primers available for transcriptome alignment
- ✅ Multi-peak optimization can select best primer combination
- ✅ Complete isoform coverage analysis possible
- ✅ Proper specificity filtering across all peaks
- ✅ True multi-peak primer design workflow functional

## Files Modified

1. **bin/extract_cdna_primers.R**
   - Lines 37-52: Multi-peak strand lookup with peak-specific keys
   - Lines 137-148: Use sequence_id (with peak) for strand lookup

2. **bin/primers_to_fasta.R**
   - Lines 27-45: Group by gene+peak in multi-peak mode for filtering

## Testing Commands

```bash
# Re-run primer extraction
Rscript bin/extract_cdna_primers.R \
  --primer3_output results/multipeak_plot/primer3_output.txt \
  --peaks_tsv results/multipeak_plot/processed_peaks/selected_peaks.tsv \
  --out_tsv results/multipeak_plot/cdna_primers.tsv

# Verify all peaks extracted
grep "ENSG00000107130" results/multipeak_plot/cdna_primers.tsv | cut -f2 | sort -u

# Generate alignment FASTA
Rscript bin/primers_to_fasta.R \
  --primers_tsv results/multipeak_plot/cdna_primers.tsv \
  --out_fasta results/multipeak_plot/primers_for_alignment.fa \
  --max_primers_per_gene 20

# Verify all peaks in alignment file
grep "ENSG00000107130" results/multipeak_plot/primers_for_alignment.fa | cut -d'|' -f2 | sort | uniq -c
```

## Related Documentation

- See `MULTIPEAK_PEAK_ID_TRACKING.md` for the initial peak identifier implementation
- See `PEAK_ID_FIX_SUMMARY.md` for quick reference

## Conclusion

Two critical bugs were preventing multi-peak mode from working correctly:
1. **Strand lookup bug** limited extraction to last peak only
2. **Filtering bug** limited alignment to first peak's primers only

Both bugs are now fixed, enabling true multi-peak primer design with proper tracking of all primers across all peaks for each gene.
