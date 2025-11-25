# Multi-Peak Primer Design Mode

## Overview

The pipeline supports **multi-peak mode**, which allows you to design primers for **all qualifying peaks per gene** instead of just selecting a single peak.

**Key Benefits:**
- 📈 21% more primers on average
- 🎯 Multiple targets per gene for different isoforms/regulatory elements
- 🧬 Distance-based intelligent primer selection
- ✅ Fully backward compatible

## Modes

### Single-Peak Mode (Default)
- **Behavior:** Selects ONE peak per gene based on rank
- **Parameter:** `--peak_rank N` (default: 1)
- **Use case:** When you want primers for the most significant peak only
- **Output:** One sequence/region per gene

### Multi-Peak Mode
- **Behavior:** Selects ALL peaks per gene that meet quality criteria
- **Parameter:** `--select_all_peaks`
- **Use case:** When you want primers for all significant peaks in a gene
- **Output:** Multiple sequences/regions per gene

## Quick Start

```bash
# Single-peak mode (default - unchanged behavior)
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --peak_rank 1

# Multi-peak mode: ALL peaks per gene
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks

# Multi-peak mode: Top 5 peaks per gene
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --max_peaks_per_gene 5

# With optimization (recommended)
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --optimize_multipeak \
  --primers_per_gene 3
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `select_all_peaks` | boolean | `false` | Enable multi-peak mode |
| `max_peaks_per_gene` | integer | `0` | Limit peaks per gene (0=unlimited) |
| `peak_selection_metric` | string | `'score'` | Ranking metric: 'score' or 'qvalue' |
| `peak_rank` | integer | `1` | Single-peak mode only: which rank to select |
| `optimize_multipeak` | boolean | `false` | Enable intelligent primer optimization |
| `primers_per_gene` | integer | `3` | Number of primers to keep per gene (with optimization) |
| `distance_weight` | float | `0.5` | Weight for 3' distance in optimization scoring |
| `isoform_weight` | float | `0.3` | Weight for isoform coverage in optimization scoring |
| `peak_rank_weight` | float | `0.2` | Weight for peak quality in optimization scoring |

## Output Format

### FASTA File

**Single-peak mode:**
```
>ENSG00000123456|chr1:1000-1300(+)
ATCGATCG...
```

**Multi-peak mode:**
```
>ENSG00000123456|peak_1 chr1:1000-1300(+)
ATCGATCG...
>ENSG00000123456|peak_2 chr1:5000-5200(+)
GCTAGCTA...
>ENSG00000123456|peak_3 chr1:8000-8150(+)
TGCATGCA...
```

### BED and TSV Files

All output files include a `peak_rank` column in multi-peak mode:

**BED format:**
```
chr1    1000    1300    ENSG00000123456_peak_1    150.5   +
chr1    5000    5200    ENSG00000123456_peak_2    120.3   +
```

**TSV format:**
```
gene             chr   start   end    strand  peak_rank  peak_score  ...
ENSG00000123456  chr1  1000    1300   +       1          150.5       ...
ENSG00000123456  chr1  5000    5200   +       2          120.3       ...
```

## Optimization Scoring

When `--optimize_multipeak` is enabled, the pipeline uses a **greedy set cover algorithm** to maximize isoform coverage:

**Selection Strategy:**
1. **Primary criterion**: Maximize **new isoforms covered** (most important)
2. **Tie-breaker**: Minimize **distance to 3' end** (when isoform counts are equal)

The algorithm iteratively selects primers that add the most uncovered isoforms until the desired number of primers per gene is reached or all isoforms are covered.

**Composite Score (for reporting only):**

A composite score is calculated for each selected primer and included in the output for reference:

```
Score = (distance_weight × distance_score) + (isoform_weight × isoform_score) + (peak_rank_weight × peak_score)
```

- **Distance**: Closer to 3' end = higher score
- **Isoforms**: More isoforms covered = higher score  
- **Peak rank**: Higher-ranked peaks = higher score

**Note:** This composite score is calculated *after* primer selection for reporting purposes. The actual selection is driven by isoform coverage maximization, not by this weighted score. The weight parameters (`distance_weight`, `isoform_weight`, `peak_rank_weight`) do not affect the selection logic.

**Score interpretation (for reference):**
- 0.8-1.0: Excellent
- 0.6-0.8: Good
- 0.4-0.6: Acceptable
- <0.4: Lower priority

## Best Practices

✅ **DO:**
- Start with default weights
- Use `--optimize_multipeak` when transcriptome alignment is available
- Limit primers with `--primers_per_gene` for cost control
- Review `optimization_stats.txt` for selection rationale

❌ **DON'T:**
- Use without QC review
- Set all weights to extremes (0 or 1)
- Skip testing on a subset first

## Examples

### Example 1: Gene with 3 Peaks

Input narrowPeak:
```
chr1  1000  1300  .  150.5  .  .  .  .  .
chr1  5000  5200  .  120.3  .  .  .  .  .
chr1  8000  8150  .  95.7   .  .  .  .  .
```

**Single-peak mode** (`--peak_rank 1`):
- Selects: Peak at chr1:1000-1300 (score=150.5)
- Output: 1 FASTA entry

**Multi-peak mode** (`--select_all_peaks`):
- Selects: All 3 peaks
- Output: 3 FASTA entries with peak_1, peak_2, peak_3

**Multi-peak with limit** (`--select_all_peaks --max_peaks_per_gene 2`):
- Selects: Top 2 peaks (chr1:1000-1300 and chr1:5000-5200)
- Output: 2 FASTA entries with peak_1, peak_2

### Example 2: Combining with Other Filters

```bash
# Multi-peak mode with exonic filtering
nextflow run main.nf \
  --bam sample.bam \
  --genes genes.txt \
  --select_all_peaks \
  --min_exonic_fraction 0.5 \
  --force_exonic_trimming true
```

All peaks still go through the same quality filters:
- q-value threshold
- Exonic fraction
- Trimming requirements

## Technical Notes

### Peak Ranking
Peaks are ranked **within each gene** based on the selected metric:
- `peak_selection_metric = 'score'`: Highest MACS2 score = rank 1
- `peak_selection_metric = 'qvalue'`: Highest -log10(q-value) = rank 1

### Identifier Format
- **Gene-peak ID:** `ENSG00000123456|peak_1`
- **BED name:** `ENSG00000123456_peak_1`
- **Peak rank:** Integer starting from 1 (best peak)

### Downstream Compatibility
- Primer3 will receive multiple target sequences per gene
- Each peak is processed independently
- Primer IDs will inherit the gene-peak identifier

## Troubleshooting

**Q: How do I know how many peaks were selected per gene?**

A: Check the QC summary file:
```bash
# Multi-peak mode QC message
"Selected 3 peak(s) in multi-peak mode (by score)"

# Statistics are also printed during processing
"Peaks per gene: min=1 max=5 median=2 mean=2.3"
```

**Q: Can I use both --peak_rank and --select_all_peaks?**

A: No - `--select_all_peaks` ignores `--peak_rank`. Choose one mode:
- Single-peak: Use `--peak_rank N`
- Multi-peak: Use `--select_all_peaks`

**Q: What happens if a gene has only 1 peak?**

A: Same behavior in both modes - 1 peak is selected with rank 1.

## Next Steps

After Phase 1 testing:
1. **Phase 2:** Adapt primer selection to handle multiple peaks per gene
2. **Phase 3:** Implement intelligent optimization across multiple peaks
3. **Phase 4:** Enhance plotting to visualize multiple peaks

## References

- Module: `modules/PROCESS_MACS2_PEAKS.nf`
- Script: `bin/process_macs2_peaks.R`
