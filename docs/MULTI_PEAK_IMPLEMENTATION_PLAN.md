# Multi-Peak Primer Design Implementation Plan

## Executive Summary

**Goal:** Extend the pipeline to design primers for **all qualifying peaks per gene** instead of just one peak.

**Current Behavior:** `--peak_rank` parameter selects a single peak per gene (rank 1, 2, 3, etc.)

**New Behavior:** New `--select_all_peaks` mode that designs primers for all peaks meeting quality criteria.

**Impact Areas:**
1. ✅ Peak selection logic (`process_macs2_peaks.R`)
2. ✅ Primer3 execution (already handles multiple targets per gene)
3. ⚠️ Primer selection logic (`select_best_primer.py`)
4. ⚠️ Isoform optimization (`optimize_primer_isoforms.py`)
5. ⚠️ Plotting (`MakePlots_new.R`)
6. ⚠️ Output file structure

---

## Part 1: Detailed Analysis

### 1.1 Current Architecture

```
MACS2_CALLPEAK
    ↓
PROCESS_MACS2_PEAKS (Selects 1 peak/gene using --peak_rank)
    ↓ primer_targets.fa (1 sequence per gene)
    ↓ primer_targets.bed (1 region per gene)
    ↓
RUN_PRIMER3 (Designs primers for each target)
    ↓
PRIMERS_TO_FASTA (Extracts primer sequences)
    ↓
ALIGN_PRIMERS_TRANSCRIPTOME (Optional: aligns all primers)
    ↓
SELECT_BEST_PRIMERS (3-stage filtering → 1 primer per gene)
    ↓
OPTIMIZE_PRIMER_ISOFORMS (Optional: 1 optimal primer per gene)
    ↓
MAKEPLOTS_NEW (Plots 1 peak per gene)
```

### 1.2 Proposed Architecture (Multi-Peak Mode)

```
MACS2_CALLPEAK
    ↓
PROCESS_MACS2_PEAKS_MULTI (NEW: Selects ALL peaks/gene)
    ↓ primer_targets_multi.fa (N sequences per gene, e.g., ENSG001|peak1, ENSG001|peak2)
    ↓ primer_targets_multi.bed (N regions per gene)
    ↓
RUN_PRIMER3 (Designs primers for each target - no changes needed!)
    ↓
PRIMERS_TO_FASTA (Extracts primer sequences - works as-is)
    ↓
ALIGN_PRIMERS_TRANSCRIPTOME (Optional: aligns all primers - works as-is)
    ↓
SELECT_BEST_PRIMERS_MULTI (NEW: Filters primers per peak-gene combination)
    ↓ best_primers_per_peak.tsv (N rows per gene, one per qualifying peak)
    ↓
OPTIMIZE_PRIMERS_MULTIPEAK (NEW: Select best primer combination)
    ↓ optimized_primers_multi.tsv (Strategy: distance scoring + isoform coverage)
    ↓
MAKEPLOTS_MULTI (NEW: Plot all peaks per gene or selected subset)
```

---

## Part 2: Step-by-Step Implementation Plan

### Phase 1: Core Peak Selection Logic ⭐ PRIORITY ✅ COMPLETE

#### Step 1.1: Modify `process_macs2_peaks.R` ✅ DONE

**Status:** IMPLEMENTED AND TESTED (October 22, 2025)

**Location:** `bin/process_macs2_peaks.R`, lines 229-290

**Changes:**

1. **Add new parameter:**
   ```r
   make_option("--select_all_peaks", action="store_true", default=FALSE, 
               help="Select all peaks per gene instead of single peak (ignores --peak_rank)")
   ```

2. **Modify peak selection logic:**
   ```r
   if (opt$select_all_peaks) {
     # NEW: Select ALL peaks that meet criteria
     cat("Multi-peak mode: selecting ALL peaks per gene\n")
     
     selected_peaks <- ordered_peaks  # All ordered peaks
     
     # Add peak_rank column showing ranking (1=best, 2=second, etc.)
     selected_peaks <- selected_peaks %>%
       group_by(gene_id) %>%
       mutate(peak_rank = row_number()) %>%
       ungroup()
     
     cat("Selected", nrow(selected_peaks), "peaks for", 
         length(unique(selected_peaks$gene_id)), "genes\n")
     cat("  Peaks per gene: min=", 
         min(table(selected_peaks$gene_id)), 
         " max=", max(table(selected_peaks$gene_id)), 
         " median=", median(table(selected_peaks$gene_id)), "\n")
     
   } else {
     # EXISTING: Select single peak by rank
     # ... existing code ...
   }
   ```

3. **Update output file generation:**
   - Modify FASTA headers to include peak identifier: `>ENSG00000123456|peak_1`
   - Modify BED file to include peak rank in name column: `ENSG00000123456_peak_1`
   - Update TSV to include `peak_rank` column

4. **Update QC summary:**
   - Track how many peaks selected per gene
   - Store statistics for all peaks, not just "best" peak
   - New columns: `num_peaks_selected`, `peak_ranks_selected`

**Estimated Effort:** 4-6 hours

**Testing:**
```bash
# Test with sample data
Rscript bin/process_macs2_peaks.R \
  --peaks peaks.narrowPeak \
  --gtf genome.gtf \
  --genes genes.txt \
  --select_all_peaks \
  --out_fa test_multi.fa \
  --out_bed test_multi.bed \
  --out_peaks test_multi.tsv \
  --out_qc test_qc.tsv
  
# Verify: Each gene should have multiple entries
awk -F'\t' '{print $1}' test_multi.tsv | sort | uniq -c
```

---

### Phase 2: Primer Selection Adaptation ✅ NO CHANGES NEEDED

#### Step 2.1: Analysis of `select_best_primer.py` ✅ COMPLETE

**Status:** ANALYZED - NO MODIFICATIONS REQUIRED (October 22, 2025)

**Finding:** The script already works correctly with multi-peak format because:
- Groups by full `primer_id` which includes peak identifier
- Each gene-peak combination gets unique primer_ids
- No assumption of one primer per base gene
- Format: `ENSG00000197756_peak_1|1|F` is treated as independent from `ENSG00000197756_peak_2|1|F`

**Conclusion:** Automatically compatible with Phase 1 output.

**Estimated Effort:** 0 hours (analysis only)

---

### Phase 3: Isoform Optimization for Multi-Peak ✅ COMPLETE

#### Step 3.1: Create `optimize_primers_multipeak.py` ✅ DONE

**Status:** IMPLEMENTED (October 22, 2025)

**Location:** `bin/select_best_primer.py`

**Current Behavior:** Assumes 1 target region per gene, selects 1 best primer per gene

**New Behavior:** Process each gene-peak combination independently

**Changes:**

1. **Update gene identification logic:**
   ```python
   # Parse gene_id and peak_rank from primer_id
   # Format: ENSG00000123456|peak_1
   def parse_primer_id(primer_id):
       parts = primer_id.split('|')
       if len(parts) == 2:
           gene_id, peak_id = parts
           peak_rank = int(peak_id.split('_')[1])
           return gene_id, peak_rank
       else:
           # Backward compatible: single peak mode
           return primer_id, 1
   ```

2. **Modify filtering to work per gene-peak combination:**
   ```python
   # Group primers by (gene_id, peak_rank) instead of just gene_id
   primer_groups = primers_df.groupby(['gene_id', 'peak_rank'])
   
   # Apply 3-stage filtering per group
   for (gene_id, peak_rank), group in primer_groups:
       # Stage 1: Perfect alignment to target gene
       # Stage 2: 3' distance filtering
       # Stage 3: Unique gene mapping
       best_primer = select_best_for_peak(group)
       results.append(best_primer)
   ```

3. **Update output format:**
   - Add `peak_rank` column to output
   - Keep `gene_id` column for backwards compatibility
   - Add `gene_peak_id` column: `ENSG00000123456_peak_1`

**Estimated Effort:** 3-4 hours

---

### Phase 3: Isoform Optimization for Multi-Peak

#### Step 3.1: Create `optimize_primers_multipeak.py`

**Location:** `bin/optimize_primers_multipeak.py` (NEW FILE)

**Goal:** Select best primer combination across all peaks considering:
1. Distance to 3' end (shorter preferred)
2. Isoform coverage (maximize distinct isoforms)
3. Peak quality (prefer higher-ranked peaks)

**Algorithm:**

```python
"""
Multi-Peak Primer Optimization Algorithm

Input: best_primers_per_peak.tsv with columns:
  - gene_id
  - peak_rank
  - primer_id
  - distance_to_end_min
  - isoforms_targeted
  - peak_score

Strategy:
1. For each gene, we have N peaks with M primers each
2. Score each primer using weighted criteria:
   
   primer_score = (
       w1 * distance_score +     # Prefer closer to 3' end
       w2 * isoform_score +       # Prefer more isoforms
       w3 * peak_rank_score       # Prefer higher-ranked peaks
   )
   
   where:
     distance_score = 1 - (distance / max_distance)  # 0-1, higher is better
     isoform_score = isoforms_targeted / max_isoforms
     peak_rank_score = 1 / peak_rank  # Rank 1 = 1.0, rank 2 = 0.5, etc.
     
   Default weights: w1=0.5, w2=0.3, w3=0.2

3. Within each gene:
   - Score all primers across all peaks
   - Select top N primers (configurable, default=3)
   - Ensure no isoform overlap (similar to current optimization)

4. Across genes:
   - Resolve isoform conflicts (same as current optimization)
   - Prefer primers with higher scores
"""

class MultiPeakPrimerOptimizer:
    def __init__(self, distance_weight=0.5, isoform_weight=0.3, 
                 peak_rank_weight=0.2):
        self.w_dist = distance_weight
        self.w_iso = isoform_weight
        self.w_rank = peak_rank_weight
    
    def score_primer(self, primer_row, max_distance, max_isoforms):
        """Calculate composite score for a primer."""
        # Distance score (closer = better)
        dist_score = 1.0 - (primer_row['distance_to_end_min'] / max_distance)
        dist_score = max(0.0, min(1.0, dist_score))
        
        # Isoform score (more = better)
        iso_score = primer_row['isoforms_targeted'] / max_isoforms if max_isoforms > 0 else 0
        
        # Peak rank score (lower rank number = better)
        rank_score = 1.0 / primer_row['peak_rank']
        
        # Weighted combination
        total = (self.w_dist * dist_score + 
                 self.w_iso * iso_score + 
                 self.w_rank * rank_score)
        
        return total
    
    def optimize(self, primers_df, primers_per_gene=3):
        """
        Select best primers per gene considering multiple peaks.
        
        Returns:
        - DataFrame with top primers per gene
        - Maintains isoform uniqueness
        - Sorted by composite score
        """
        results = []
        
        # Calculate max values for normalization
        max_dist = primers_df['distance_to_end_min'].max()
        max_iso = primers_df['isoforms_targeted'].max()
        
        # Score all primers
        primers_df['composite_score'] = primers_df.apply(
            lambda row: self.score_primer(row, max_dist, max_iso),
            axis=1
        )
        
        # Select top primers per gene
        for gene_id, gene_group in primers_df.groupby('gene_id'):
            # Sort by score (descending)
            top_primers = gene_group.nlargest(primers_per_gene, 'composite_score')
            
            # Store with selection reason
            top_primers['selection_reason'] = [
                f"Rank {i+1}: score={row['composite_score']:.3f} "
                f"(dist={row['distance_to_end_min']}bp, "
                f"iso={row['isoforms_targeted']}, "
                f"peak={row['peak_rank']})"
                for i, (_, row) in enumerate(top_primers.iterrows())
            ]
            
            results.append(top_primers)
        
        optimized = pd.concat(results)
        
        # Resolve isoform conflicts (reuse existing logic)
        optimized = self.resolve_isoform_conflicts(optimized)
        
        return optimized
```

**Parameters:**
- `--primers_per_gene`: Number of primers to select per gene (default: 3)
- `--distance_weight`: Weight for distance score (default: 0.5)
- `--isoform_weight`: Weight for isoform coverage (default: 0.3)
- `--peak_rank_weight`: Weight for peak rank (default: 0.2)

**Estimated Effort:** 8-10 hours

---

### Phase 4: Plotting Modifications

#### Step 4.1: Extend `MakePlots_new.R`

**Location:** `bin/MakePlots_new.R`

**Current Behavior:** Plots single peak per gene

**New Behavior:** Plot all peaks for a gene (or selected subset)

**Changes:**

1. **Update peak reading logic:**
   ```r
   # Read peaks TSV
   peaks <- fread(peaks_tsv_file)
   
   # Filter to gene - may return multiple peaks
   gene_peaks <- peaks[peaks$gene_id == gene_id, ]
   
   if (nrow(gene_peaks) == 0) {
     stop("No peaks found for gene ", gene_id)
   }
   
   cat(sprintf("Found %d peak(s) for gene %s\n", nrow(gene_peaks), gene_id))
   ```

2. **Modify plotting to show multiple peaks:**
   ```r
   # Add peak regions as separate tracks
   for (i in 1:nrow(gene_peaks)) {
     peak_gr <- GRanges(
       seqnames = gene_peaks$peak_chr[i],
       ranges = IRanges(gene_peaks$peak_start[i], gene_peaks$peak_end[i])
     )
     
     # Add to plot with rank label
     peak_track <- geom_rect(
       aes(xmin = start(peak_gr), xmax = end(peak_gr), 
           ymin = -Inf, ymax = Inf),
       fill = peak_colors[i], alpha = 0.2
     )
     
     peak_label <- annotate(
       "text", x = (start(peak_gr) + end(peak_gr))/2, y = Inf,
       label = paste0("Peak ", gene_peaks$peak_rank[i]),
       vjust = 1.5, size = 3
     )
   }
   ```

3. **Update primer display:**
   ```r
   # Read primers - may now have multiple per gene
   primers <- fread(primer_bed_file)
   gene_primers <- primers[grepl(gene_id, primers$name), ]
   
   # Color-code primers by peak_rank
   gene_primers$peak_rank <- str_extract(gene_primers$name, "peak_(\\d+)")
   gene_primers$color <- peak_colors[gene_primers$peak_rank]
   
   # Plot primers with peak-specific colors
   primer_track <- geom_segment(
     data = gene_primers,
     aes(x = start, xend = end, y = primer_y, yend = primer_y, 
         color = peak_rank),
     linewidth = 2
   )
   ```

**Alternative Approach:** Generate separate plot per peak
```r
# Instead of combined plot, generate ENSG001_peak1.png, ENSG001_peak2.png, etc.
for (i in 1:nrow(gene_peaks)) {
  out_file <- sprintf("%s_peak%d.png", gene_id, gene_peaks$peak_rank[i])
  plot_single_peak(gene_peaks[i, ], out_file)
}
```

**Estimated Effort:** 6-8 hours

---

### Phase 5: Workflow Integration

#### Step 5.1: Create New Workflow Mode

**Option A: Extend Existing Workflow**
```groovy
// workflows/primer_design.nf

// Add parameter check
if (params.select_all_peaks) {
  cat("Running in multi-peak mode\n")
}

// Modify PROCESS_MACS2_PEAKS call
PROCESS_MACS2_PEAKS(
  MACS2_CALLPEAK.out.narrowpeak,
  genes_ch,
  gtf_ch,
  fasta_ch,
  params.select_all_peaks  // Pass new parameter
)

// Rest of workflow adapts automatically due to file format compatibility
```

**Option B: Create Separate Workflow** (RECOMMENDED)
```groovy
// workflows/primer_design_multipeak.nf

workflow primer_design_multipeak {
  take:
    bam_ch
    genes_ch
    gtf_ch
    fasta_ch
  
  main:
    // Use all existing modules with new parameters
    MACS2_CALLPEAK(bam_ch)
    
    PROCESS_MACS2_PEAKS_MULTI(
      MACS2_CALLPEAK.out.narrowpeak,
      genes_ch,
      gtf_ch,
      fasta_ch
    )
    
    // Primer3, alignment, etc. work unchanged
    RUN_PRIMER3(PROCESS_MACS2_PEAKS_MULTI.out.fasta)
    
    // ... rest of workflow
    
    // New optimization step
    if (params.optimize_multipeak) {
      OPTIMIZE_PRIMERS_MULTIPEAK(
        best_primers_per_peak,
        alignment_summary
      )
    }
  
  emit:
    primers_per_peak = best_primers_per_peak
    optimized_primers = optimized_primers_multi
    // ... other outputs
}
```

**Estimated Effort:** 4-6 hours

---

## Part 3: Configuration & Parameters

### New Parameters (add to `params.config`)

```groovy
params {
  // Multi-peak mode
  select_all_peaks = false          // Enable multi-peak mode
  max_peaks_per_gene = 0            // 0 = unlimited, N = limit to top N
  
  // Multi-peak optimization
  optimize_multipeak = false        // Enable multi-peak optimization
  primers_per_gene = 3              // Number of primers to select per gene
  distance_weight = 0.5             // Weight for distance scoring
  isoform_weight = 0.3              // Weight for isoform coverage
  peak_rank_weight = 0.2            // Weight for peak rank
  
  // Plotting
  plot_all_peaks = true             // Plot all peaks (true) or separate files (false)
  max_peaks_to_plot = 5             // Limit number of peaks per plot
}
```

### Command Line Examples

```bash
# Multi-peak mode: Design primers for all peaks
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --outdir results_multi

# With optimization: Select top 3 primers per gene
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --optimize_multipeak \
  --primers_per_gene 3 \
  --distance_weight 0.6 \
  --isoform_weight 0.3 \
  --peak_rank_weight 0.1

# Limit to top 5 peaks per gene
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --max_peaks_per_gene 5
```

---

## Part 4: Output File Structure

### Multi-Peak Mode Output

```
results_multi/
├── processed_peaks/
│   ├── primer_targets_multi.fa      # N sequences per gene
│   ├── primer_targets_multi.bed     # N regions per gene
│   ├── selected_peaks_multi.tsv     # All peaks with ranks
│   └── peaks_qc_summary_multi.tsv   # QC with peak counts
│
├── primer3_output/
│   └── ENSG*_peak*.txt             # One file per peak
│
├── primers/
│   ├── primers_all.fasta           # All primers for all peaks
│   ├── best_primers_per_peak.tsv   # Best primer per peak
│   └── optimized_primers_multi.tsv # Top N primers per gene
│
├── alignment/
│   ├── primer_alignment_summary.tsv
│   └── primer_alignment_report.tsv
│
└── plots/
    ├── plot_ENSG001_all_peaks.png  # Combined plot (if plot_all_peaks=true)
    ├── plot_ENSG001_peak1.png      # Individual plots (if plot_all_peaks=false)
    ├── plot_ENSG001_peak2.png
    └── plot_ENSG001_peak3.png
```

### File Format Examples

**primer_targets_multi.fa:**
```
>ENSG00000123456|peak_1 chr1:1000-1300(+) rank=1 score=150
ATCGATCGATCG...
>ENSG00000123456|peak_2 chr1:5000-5200(+) rank=2 score=120
GCTAGCTAGCTA...
>ENSG00000234567|peak_1 chr2:2000-2400(-) rank=1 score=180
TGCATGCATGCA...
```

**best_primers_per_peak.tsv:**
```
gene_id           peak_rank  primer_id                 distance_to_end  isoforms_targeted  composite_score
ENSG00000123456   1          ENSG00000123456|peak_1    150              5                  0.85
ENSG00000123456   2          ENSG00000123456|peak_2    200              3                  0.72
ENSG00000234567   1          ENSG00000234567|peak_1    100              8                  0.92
```

**optimized_primers_multi.tsv:**
```
gene_id           selected_primers                                          selection_reason
ENSG00000123456   ENSG00000123456|peak_1,ENSG00000123456|peak_2            Rank 1: score=0.85 (dist=150bp, iso=5, peak=1); Rank 2: score=0.72 (dist=200bp, iso=3, peak=2)
ENSG00000234567   ENSG00000234567|peak_1                                    Rank 1: score=0.92 (dist=100bp, iso=8, peak=1)
```

---

## Part 5: Testing Strategy

### Test Cases

#### Test 1: Single Peak per Gene (Baseline)
```bash
# Should behave identically to current workflow
nextflow run main.nf --bam test.bam --genes genes.txt --peak_rank 1
nextflow run main.nf --bam test.bam --genes genes.txt --select_all_peaks --max_peaks_per_gene 1
# Compare outputs - should be identical
```

#### Test 2: Multiple Peaks per Gene
```bash
# Gene with 3 peaks
nextflow run main.nf --bam test.bam --genes genes_multipeak.txt --select_all_peaks

# Verify:
# 1. primer_targets.fa has 3 entries for gene
# 2. Primer3 ran 3 times
# 3. 3 sets of primers generated
```

#### Test 3: Optimization Scoring
```bash
# Test scoring system with known data
python bin/optimize_primers_multipeak.py \
  --best_primers test_primers_per_peak.tsv \
  --primers_per_gene 2 \
  --distance_weight 1.0 --isoform_weight 0.0  # Only distance
  
# Verify: Selected primers have shortest distances
```

#### Test 4: Plotting
```bash
# Generate plots for multi-peak gene
nextflow run main.nf --bam test.bam --genes single_gene.txt \
  --select_all_peaks --makeplots

# Verify:
# - Plot shows all peaks
# - Primers color-coded by peak
# - Legend indicates peak ranks
```

---

## Part 6: Implementation Timeline

### Phase 1: Core Peak Selection (Week 1)
- **Days 1-2:** Modify `process_macs2_peaks.R` for multi-peak selection
- **Day 3:** Test peak selection with sample data
- **Day 4:** Update output file formats
- **Day 5:** Integration testing

### Phase 2: Primer Selection (Week 2)
- **Days 1-2:** Modify `select_best_primer.py` for per-peak filtering
- **Day 3:** Test with multi-peak data
- **Days 4-5:** Update documentation and error handling

### Phase 3: Optimization (Week 3)
- **Days 1-3:** Implement `optimize_primers_multipeak.py`
- **Day 4:** Test scoring algorithm with various weights
- **Day 5:** Integration with existing isoform optimization

### Phase 4: Plotting & Workflow (Week 4)
- **Days 1-2:** Modify `MakePlots_new.R` for multi-peak display
- **Day 3:** Create/update workflow
- **Days 4-5:** End-to-end testing and documentation

---

## Part 7: Risk Assessment & Mitigation

### Risk 1: Primer3 Performance
**Issue:** Running Primer3 for many peaks may be slow
**Mitigation:** 
- Add `max_peaks_per_gene` parameter
- Consider parallelizing Primer3 calls
- Profile and optimize if needed

### Risk 2: Output File Size
**Issue:** Multi-peak mode generates N× more files
**Mitigation:**
- Implement file compression
- Add summary files
- Clear documentation on disk usage

### Risk 3: Plotting Complexity
**Issue:** Too many peaks make plots unreadable
**Mitigation:**
- Add `max_peaks_to_plot` parameter
- Option for separate plots per peak
- Smart layout algorithm

### Risk 4: Backward Compatibility
**Issue:** Changes may break existing workflows
**Mitigation:**
- Default to single-peak mode (`select_all_peaks=false`)
- Maintain existing parameter names
- Comprehensive testing

---

## Part 8: Documentation Requirements

### User Documentation
1. **Feature Guide:** `docs/features/MULTI_PEAK_MODE.md`
   - When to use multi-peak mode
   - Parameter descriptions
   - Example commands
   - Output interpretation

2. **Update Existing Docs:**
   - `docs/features/PEAK_RANKING.md` - Add multi-peak section
   - `docs/TECHNICAL_APPENDIX.md` - Update file formats
   - `docs/pipeline_steps.md` - Add multi-peak workflow

### Developer Documentation
1. **Implementation Notes:** `docs/MULTI_PEAK_IMPLEMENTATION.md`
   - Architecture decisions
   - Algorithm details
   - Testing procedures

---

## Part 9: Success Criteria

### Functional Requirements
- ✅ Select all qualifying peaks per gene
- ✅ Design primers for each peak
- ✅ Intelligently select best primer combination
- ✅ Generate informative plots
- ✅ Maintain backward compatibility

### Performance Requirements
- ✅ Handle genes with 10+ peaks
- ✅ Complete 38-gene analysis in <2 hours (multi-peak mode)
- ✅ Plotting parallelization still works

### Quality Requirements
- ✅ Distance-based scoring gives sensible results
- ✅ Isoform optimization works with multiple primers per gene
- ✅ Plots are readable and informative
- ✅ All tests pass

---

## Part 10: Summary & Recommendations

### Recommended Approach

**Option: Incremental Implementation** ✅ RECOMMENDED

1. **Start with:** Peak selection logic (Phase 1)
   - Lowest risk, high value
   - Can test independently
   - Provides foundation for other changes

2. **Then add:** Primer selection adaptation (Phase 2)
   - Natural next step
   - Reuses existing filtering logic

3. **Then add:** Optimization scoring (Phase 3)
   - Most complex part
   - Requires phases 1-2 complete
   - Can be optional initially

4. **Finally add:** Plotting enhancements (Phase 4)
   - Nice-to-have, not critical
   - Can use simple approach initially

### Key Design Decisions

1. **Use separate workflow?** 
   - **NO** - Extend existing workflow with parameter flag
   - Easier to maintain, better code reuse

2. **Optimize all primers or subset?**
   - **SUBSET** - Use scoring to select top N
   - More practical, avoids overwhelming users

3. **Combined or separate plots?**
   - **CONFIGURABLE** - Support both modes
   - Combined for overview, separate for detail

4. **Backward compatibility?**
   - **YES** - Default to single-peak mode
   - Existing workflows unaffected

### Estimated Total Effort

- **Development:** 3-4 weeks (full-time)
- **Testing:** 1 week
- **Documentation:** 3-4 days
- **Total:** ~1 month

### Next Steps

1. **Review this plan** with team/users
2. **Create GitHub issues** for each phase
3. **Set up test data** with multi-peak genes
4. **Begin Phase 1** implementation
5. **Iterate** based on testing results

---

*Plan created: October 21, 2025*
*Status: Ready for implementation*
