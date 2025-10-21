# Peak Ranking Explained

## Quick Answer

**By default, peaks are ranked by MACS2 score (column 5), NOT by q-value.**

- **Rank 1 peak** = Peak with **highest MACS2 score** for that gene
- **Rank 2 peak** = Peak with **second-highest MACS2 score** for that gene
- And so on...

You can change this behavior to rank by q-value instead using the `peak_selection_metric` parameter.

---

## Configuration

### Current Settings (params.config)

```groovy
peak_selection_metric = 'score'   // Default: rank by MACS2 score
peak_rank = 1                     // Default: select the top-ranked peak
```

### Available Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `peak_selection_metric` | `'score'` | `'score'` or `'qvalue'` | Which metric to use for ranking peaks |
| `peak_rank` | `1` | Any positive integer | Which ranked peak to select (1=best, 2=second-best, etc.) |

---

## Understanding MACS2 Peak Columns

MACS2 narrowPeak format contains these relevant columns:

```
Col 1:  chrom          - Chromosome
Col 2:  chromStart     - Start position (0-based)
Col 3:  chromEnd       - End position
Col 4:  name           - Peak name
Col 5:  score          - Integer score (0-1000) ← USED BY DEFAULT
Col 6:  strand         - Strand (usually ".")
Col 7:  signalValue    - Measurement of enrichment
Col 8:  pValue         - -log10(p-value)
Col 9:  qValue         - -log10(q-value) ← ALTERNATIVE RANKING
Col 10: peak           - Peak summit position relative to start
```

### Example Peak:

```
chr1  14375  14602  peak_1  174  .  11.95  21.12  17.50  124
                             ^^^              ^^^   ^^^^^
                           score           pValue  qValue
```

---

## Ranking Methods

### Method 1: By MACS2 Score (Default)

**What is MACS2 score?**
- Integer value from 0-1000
- Represents the overall peak quality/strength
- Calculated by MACS2 as: `min(int(-10 * log10(qvalue)), 1000)`
- Effectively a capped, scaled version of the q-value

**How ranking works:**
```r
# Peaks are sorted in descending order by score
ordered_peaks <- peak_gene_map[order(peak_gene_map$gene_id, -score), ]
```

**Example for gene ENSG00000123456:**
```
Peak A: score = 174  → Rank 1 ← Selected by default
Peak B: score = 120  → Rank 2
Peak C: score = 56   → Rank 3
Peak D: score = 50   → Rank 4
```

### Method 2: By Q-value

**What is q-value?**
- False Discovery Rate (FDR) adjusted p-value
- MACS2 stores it as **-log10(q-value)**
- Higher -log10(q-value) = more significant peak
- More statistically rigorous than score

**How ranking works:**
```r
# Peaks are sorted in descending order by -log10(q-value)
ordered_peaks <- peak_gene_map[order(peak_gene_map$gene_id, -qvalue), ]
```

**Example for gene ENSG00000123456:**
```
Peak A: qvalue = 17.50  → Rank 1
Peak B: qvalue = 8.76   → Rank 2
Peak C: qvalue = 5.64   → Rank 3
Peak D: qvalue = 5.04   → Rank 4
```

**To use q-value ranking:**
```groovy
// In params.config
peak_selection_metric = 'qvalue'
```

Or via command line:
```bash
nextflow run main.nf --peak_selection_metric qvalue
```

---

## Selecting Alternative Peaks

### Get the Second-Best Peak

If the top peak fails QC or you want to try an alternative region:

```groovy
// In params.config
peak_rank = 2  // Select second-best peak
```

Or via command line:
```bash
nextflow run main.nf --peak_rank 2
```

### Behavior When Peak Rank Not Available

If a gene has fewer peaks than the requested rank:

```r
# Example: You request peak_rank = 3, but gene has only 2 peaks
# Result: Peak rank 1 is selected, with a warning:

Warning: Gene ENSG00000123456 has only 2 peak(s), cannot select rank 3
         - using rank 1 instead
```

---

## Implementation Details

### Code Location

The ranking logic is in `bin/process_macs2_peaks.R`:

```r
# Lines 229-238
if (opt$peak_selection_metric == "score") {
  # Higher score is better (descending order)
  ordered_peaks <- peak_gene_map[order(peak_gene_map$gene_id, 
                                       -peaks_gr[peak_gene_map$peak_idx]$score), ]
  cat("Selecting peaks by highest MACS2 score (rank ", opt$peak_rank, ")\n")
} else if (opt$peak_selection_metric == "qvalue") {
  # Higher -log10(qvalue) is better (descending order)
  ordered_peaks <- peak_gene_map[order(peak_gene_map$gene_id, 
                                       -peaks_gr[peak_gene_map$peak_idx]$qvalue), ]
  cat("Selecting peaks by highest -log10(q-value) (rank ", opt$peak_rank, ")\n")
}

# Lines 240-255: Select peak at specified rank
for(gene in unique(ordered_peaks$gene_id)) {
  gene_peaks <- ordered_peaks[ordered_peaks$gene_id == gene, ]
  
  if(nrow(gene_peaks) >= opt$peak_rank) {
    # Select the peak at the specified rank
    peak_row <- gene_peaks[opt$peak_rank, ]
    peak_row$actual_peak_rank <- opt$peak_rank
    selected_peaks <- rbind(selected_peaks, peak_row)
  } else {
    # Not enough peaks - use rank 1 instead
    peak_row <- gene_peaks[1, ]
    peak_row$actual_peak_rank <- 1
    selected_peaks <- rbind(selected_peaks, peak_row)
    # Warning issued
  }
}
```

---

## Practical Considerations

### When to Use Score vs. Q-value

**Use Score (Default):**
- ✅ General-purpose peak selection
- ✅ When you want peaks with good overall quality
- ✅ Score is already FDR-adjusted (capped at 1000)
- ✅ Simpler interpretation

**Use Q-value:**
- ✅ When statistical significance is critical
- ✅ For publication-quality results requiring strict FDR control
- ✅ When comparing across different experiments
- ✅ More fine-grained ranking (not capped at 1000)

### Why Score is Default

1. **Score is derived from q-value**: `score = min(int(-10 * log10(qvalue)), 1000)`
2. **Practically equivalent**: For most peaks, score and qvalue give similar rankings
3. **Easier to interpret**: Integer scale 0-1000 vs. continuous -log10 scale
4. **Standard in peak calling**: ENCODE and other consortia use score

### When Rankings Might Differ

Score and q-value rankings can differ when:
- Multiple peaks have q-values that map to the same score (due to rounding/capping)
- Very high significance peaks (score capped at 1000, but q-value continues)
- Borderline peaks near significance threshold

---

## Checking Your Current Settings

### View Pipeline Log

When `PROCESS_MACS2_PEAKS` runs, it prints:

```
Selecting peaks by highest MACS2 score (rank 1)
```

Or:

```
Selecting peaks by highest -log10(q-value) (rank 2)
```

### Check Output Files

The `peaks_qc_summary.tsv` includes:
- `best_peak_score`: The MACS2 score of the selected peak
- `best_peak_pvalue`: The -log10(p-value)
- `best_peak_qvalue`: The -log10(q-value)

### Compare Rankings

To see if score vs. qvalue makes a difference for your data:

```bash
# Run with score ranking (default)
nextflow run main.nf -resume

# Run with qvalue ranking
nextflow run main.nf -resume --peak_selection_metric qvalue --outdir results/test_qvalue

# Compare the selected peaks
diff results/test/peaks_qc_summary.tsv results/test_qvalue/peaks_qc_summary.tsv
```

---

## Examples

### Example 1: Default Behavior

```bash
nextflow run main.nf --genes testdata/IMR32_C3.txt
```

**Result:**
- Ranks peaks by **MACS2 score**
- Selects **rank 1** (highest score) for each gene

### Example 2: Select Second-Best Peak by Score

```bash
nextflow run main.nf --genes testdata/IMR32_C3.txt --peak_rank 2
```

**Result:**
- Ranks peaks by **MACS2 score**
- Selects **rank 2** (second-highest score) for each gene

### Example 3: Select Best Peak by Q-value

```bash
nextflow run main.nf --genes testdata/IMR32_C3.txt --peak_selection_metric qvalue
```

**Result:**
- Ranks peaks by **-log10(q-value)**
- Selects **rank 1** (most significant) for each gene

### Example 4: Select Third-Best Peak by Q-value

```bash
nextflow run main.nf --genes testdata/IMR32_C3.txt \
  --peak_selection_metric qvalue \
  --peak_rank 3
```

**Result:**
- Ranks peaks by **-log10(q-value)**
- Selects **rank 3** (third-most significant) for each gene
- Falls back to rank 1 if gene has fewer than 3 peaks

---

## Summary Table

| Ranking Method | Metric | Range | Higher is Better | When to Use |
|----------------|--------|-------|------------------|-------------|
| **Score** (default) | MACS2 score | 0-1000 | Yes | General use, simpler |
| **Q-value** | -log10(q-value) | 0-∞ | Yes | Statistical rigor, publications |

| Peak Rank | Meaning | Command Line |
|-----------|---------|-------------|
| 1 (default) | Best peak | `--peak_rank 1` |
| 2 | Second-best peak | `--peak_rank 2` |
| 3 | Third-best peak | `--peak_rank 3` |
| N | Nth-best peak | `--peak_rank N` |

---

## Related Files

- `bin/process_macs2_peaks.R`: Implementation of peak ranking
- `modules/PROCESS_MACS2_PEAKS.nf`: Nextflow module calling the R script
- `params.config`: Default parameter settings
- Output: `peaks_qc_summary.tsv`: Contains selected peak statistics

---

*Last updated: October 20, 2025*
