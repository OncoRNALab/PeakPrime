# Multi-Peak Primer Design - User Guide# Multi-Peak Primer Design Mode



## Overview## Overview



Multi-peak mode designs primers for **ALL qualifying peaks per gene** instead of just one, enabling comprehensive experimental coverage.The pipeline now supports **multi-peak mode**, which allows you to design primers for **all qualifying peaks per gene** instead of just selecting a single peak.



**Key Benefits:****Status:** Phase 1 Complete ✅

- 📈 21% more primers on average

- 🎯 Multiple targets per gene for different isoforms/regulatory elements## Modes

- 🧬 Distance-based intelligent primer selection

- ✅ Fully backward compatible### Single-Peak Mode (Default)

- **Behavior:** Selects ONE peak per gene based on rank

## Quick Start- **Parameter:** `--peak_rank N` (default: 1)

- **Use case:** When you want primers for the most significant peak only

```bash- **Output:** One sequence/region per gene

# Basic multi-peak mode

nextflow run main.nf \### Multi-Peak Mode (New)

  --bam sample.bam \- **Behavior:** Selects ALL peaks per gene that meet quality criteria

  --genes genes.txt \- **Parameter:** `--select_all_peaks`

  --select_all_peaks- **Use case:** When you want primers for all significant peaks in a gene

- **Output:** Multiple sequences/regions per gene

# With optimization (recommended)

nextflow run main.nf \## Usage

  --bam sample.bam \

  --genes genes.txt \### Command Line

  --select_all_peaks \

  --optimize_multipeak \```bash

  --primers_per_gene 3# Single-peak mode (default - unchanged behavior)

```nextflow run main.nf \

  --bam sample.bam \

## Parameters  --genes genes.txt \

  --peak_rank 1

| Parameter | Default | Description |

|-----------|---------|-------------|# Multi-peak mode: ALL peaks per gene

| `--select_all_peaks` | `false` | Enable multi-peak mode |nextflow run main.nf \

| `--max_peaks_per_gene` | `0` | Limit peaks (0=unlimited) |  --bam sample.bam \

| `--optimize_multipeak` | `false` | Enable smart selection |  --genes genes.txt \

| `--primers_per_gene` | `3` | Primers to keep per gene |  --select_all_peaks

| `--distance_weight` | `0.5` | 3' distance importance |

| `--isoform_weight` | `0.3` | Isoform coverage importance |# Multi-peak mode: Top 5 peaks per gene

| `--peak_rank_weight` | `0.2` | Peak quality importance |nextflow run main.nf \

  --bam sample.bam \

## Output Format  --genes genes.txt \

  --select_all_peaks \

### Peak Naming  --max_peaks_per_gene 5

``````

ENSG00000197756_peak_1  ← Rank 1 (best)

ENSG00000197756_peak_2  ← Rank 2### Configuration File

ENSG00000197756_peak_3  ← Rank 3

``````groovy

params {

### Optimized Primers  // Enable multi-peak mode

```  select_all_peaks = true

base_gene_id     gene_peak_id            score   reason  

ENSG00000197756  ENSG00000197756_peak_1  0.850   Rank 1: dist=150bp, iso=5, peak=1  // Optional: limit number of peaks per gene

ENSG00000197756  ENSG00000197756_peak_2  0.720   Rank 2: dist=200bp, iso=3, peak=2  max_peaks_per_gene = 3  // 0 = unlimited (default)

```  

  // Peak ranking still applies

## Examples  peak_selection_metric = 'score'  // or 'qvalue'

}

### Example 1: Research Grade (All Peaks)```

```bash

nextflow run main.nf --select_all_peaks --optimize_multipeak --primers_per_gene 5## Parameters

```

Best for: Comprehensive studies, multiple isoform analysis| Parameter | Type | Default | Description |

|-----------|------|---------|-------------|

### Example 2: Production (Top 2 Per Gene)| `select_all_peaks` | boolean | `false` | Enable multi-peak mode |

```bash| `max_peaks_per_gene` | integer | `0` | Limit peaks per gene (0=unlimited) |

nextflow run main.nf --select_all_peaks --optimize_multipeak --primers_per_gene 2| `peak_selection_metric` | string | `'score'` | Ranking metric: 'score' or 'qvalue' |

```| `peak_rank` | integer | `1` | Single-peak mode only: which rank to select |

Best for: Cost-effective coverage, validation experiments

## Output Format

### Example 3: Distance-Prioritized  

```bash### FASTA File

nextflow run main.nf --select_all_peaks --optimize_multipeak \

  --distance_weight 0.7 --isoform_weight 0.2 --peak_rank_weight 0.1**Single-peak mode:**

``````

Best for: 3' RNA-seq protocols where distance is critical>ENSG00000123456|chr1:1000-1300(+)

ATCGATCG...

## Optimization Scoring```



Primers are scored using weighted criteria:**Multi-peak mode:**

```

```>ENSG00000123456|peak_1 chr1:1000-1300(+)

Score = 0.5×(distance) + 0.3×(isoforms) + 0.2×(peak_rank)ATCGATCG...

```>ENSG00000123456|peak_2 chr1:5000-5200(+)

GCTAGCTA...

- **Distance**: Closer to 3' end = higher score>ENSG00000123456|peak_3 chr1:8000-8150(+)

- **Isoforms**: More isoforms covered = higher score  TGCATGCA...

- **Peak rank**: Higher-ranked peaks = higher score```



**Score interpretation:**### BED File

- 0.8-1.0: Excellent

- 0.6-0.8: Good**Single-peak mode:**

- 0.4-0.6: Acceptable```

- <0.4: Lower prioritychr1    1000    1300    ENSG00000123456_peak    150.5   +

```

## Best Practices

**Multi-peak mode:**

✅ **DO:**```

- Start with default weightschr1    1000    1300    ENSG00000123456_peak_1    150.5   +

- Use `--optimize_multipeak` when transcriptome alignment is availablechr1    5000    5200    ENSG00000123456_peak_2    120.3   +

- Limit primers with `--primers_per_gene` for cost controlchr1    8000    8150    ENSG00000123456_peak_3    95.7    +

- Review `optimization_stats.txt` for selection rationale```



❌ **DON'T:**### TSV File

- Use without QC review

- Set all weights to extremes (0 or 1)All output files now include a `peak_rank` column:

- Skip testing on subset first

```

## Troubleshootinggene             chr   start   end    strand  peak_rank  peak_score  ...

ENSG00000123456  chr1  1000    1300   +       1          150.5       ...

**Problem**: No multi-peak genes  ENSG00000123456  chr1  5000    5200   +       2          120.3       ...

**Solution**: Lower `--qvalue_threshold`, check peak fileENSG00000123456  chr1  8000    8150   +       3          95.7        ...

```

**Problem**: Too many primers  

**Solution**: Use `--primers_per_gene 2` or `--max_peaks_per_gene 2`## Implementation Status



**Problem**: Optimization fails  ### ✅ Phase 1: Complete

**Solution**: Check Python environment, pandas installed- [x] Peak selection logic modified

- [x] Multi-peak parameter support

## Technical Details- [x] Output file format updated (FASTA, BED, TSV)

- [x] Peak rank tracking

- **Backward compatible**: Default behavior unchanged- [x] QC summary updated

- **Format**: gene_id_peak_N naming convention- [x] Module integration

- **Pipeline stage**: After SELECT_BEST_PRIMERS, before results- [x] Configuration parameters

- **Dependencies**: Python 3, pandas

### 🔄 Phase 2: In Progress

## Testing- [ ] Adapt primer selection for (gene, peak) combinations

- [ ] Update `select_best_primer.py`

Comprehensive tests in:

- `test_phase1_multipeak.sh` - Peak selection### ⏳ Phase 3: Planned

- `test_phase4_workflow.sh` - End-to-end integration- [ ] Multi-peak primer optimization

- [ ] Distance-based scoring

Results:- [ ] Isoform coverage optimization

- ✅ 197 peaks for 163 genes (+21%)

- ✅ Max 3 peaks per gene### ⏳ Phase 4: Planned

- ✅ Distance scoring validated- [ ] Plotting for multiple peaks

- [ ] Combined or separate visualizations

---

## Backward Compatibility

**Status**: ✅ Production Ready  

**Version**: 1.0  ✅ **Fully backward compatible** - existing workflows work unchanged:

**Date**: October 22, 2025- Default behavior: `select_all_peaks = false`

- Single-peak mode is the default

For detailed implementation, see `docs/MULTI_PEAK_IMPLEMENTATION_PLAN.md`- All existing parameters still work

- Output format unchanged in single-peak mode

## Testing

Run the Phase 1 test script:

```bash
bash test_multipeak_phase1.sh
```

This will:
1. Test single-peak mode (baseline)
2. Test multi-peak mode (all peaks)
3. Test multi-peak mode with limits
4. Verify output formats
5. Show statistics

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

- Implementation Plan: `docs/MULTI_PEAK_IMPLEMENTATION_PLAN.md`
- Test Script: `test_multipeak_phase1.sh`
- Module: `modules/PROCESS_MACS2_PEAKS.nf`
- Script: `bin/process_macs2_peaks.R`
