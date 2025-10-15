# Distance-Based Filtering for 3' RNA-seq Protocols

## Overview

PeakPrime now includes **3-stage filtering** optimized for 3' RNA-seq library preparation protocols (QuantSeq, Lexogen, etc.), where only the 3' end of transcripts is amplified and sequenced.

## Biological Rationale

In 3' RNA-seq protocols:
- **Only the 3' terminal region** of transcripts is captured during library preparation
- Off-target primer alignments **far from the 3' end** are never amplified or sequenced
- These distant alignments are **biologically irrelevant** and should not disqualify primers

**Example:**
```
Primer for gene ACTB:
✓ Perfect match to ACTB (200bp from 3' end)    → Will be amplified
✗ Perfect match to ACTG1 (5000bp from 3' end)  → Won't be amplified (ignored)

Old behavior: REJECT (multi-gene hit)
New behavior: ACCEPT (only relevant hit considered)
```

## 3-Stage Filtering Strategy

### Stage 1: Perfect Alignment Quality
- **Filter:** `mismatches == 0`
- **Purpose:** Keep only perfect matches to transcriptome
- **Removes:** Low-quality alignments with mismatches

### Stage 2: Distance to 3' End
- **Filter:** `distance_to_end <= threshold`
- **Purpose:** Keep only alignments near the 3' end (protocol-relevant)
- **Removes:** Alignments far from 3' end (not amplified)
- **Parameter:** `--distance_threshold` (default: 1000bp)

### Stage 3: Unique Gene Mapping
- **Filter:** Alignments map to only one gene
- **Purpose:** Ensure gene-specific amplification
- **Removes:** Multi-gene hits **within the relevant region**
- **Key:** Only considers alignments that passed Stage 2

## Usage

### Command Line Parameter

```bash
nextflow run main.nf \
  --bam input.bam \
  --gtf annotations.gtf \
  --genes gene_list.txt \
  --distance_threshold 1000 \
  ... other parameters ...
```

### Parameter Values

| Threshold | Description | Use Case |
|-----------|-------------|----------|
| **500-1000bp** | Conservative (default) | Immediate 3' UTR region only |
| **1500-2000bp** | Moderate | Typical 3' bias range for most protocols |
| **3000-5000bp** | Permissive | Terminal exons, long 3' UTRs |
| **>5000bp** | Very permissive | May include some non-3' regions |

**Default:** 1000bp (captures immediate 3' UTR and terminal coding sequence)

### Configuration File

Edit `params.config`:
```groovy
params {
  distance_threshold = 1000  // Change to desired value
}
```

## Output Format

The `best_primers.tsv` file now includes:

| Column | Description |
|--------|-------------|
| `aligned_gene_name` | Gene this primer maps to (with 0 mismatches near 3' end) |
| `zero_mismatch_alignments` | Number of perfect alignments near 3' end |
| **`distance_to_end_min`** | **Closest distance to 3' end (bp)** |
| `num_alignments` | Total alignments (including far hits) |

**Example output:**
```
gene_id              primer_index  aligned_gene_name  zero_mismatch_alignments  distance_to_end_min  num_alignments
ENSG00000123456      5             ACTB               3                         289                  8
```

Interpretation:
- Primer has **3 perfect alignments** to ACTB transcripts near the 3' end
- Closest alignment is **289bp** from the 3' end
- Total of **8 alignments** (including 5 that are >1000bp from 3' end, now ignored)

## Filtering Statistics

The pipeline prints detailed statistics:

```
=== 3-STAGE FILTERING FOR 3' RNA-SEQ ===
Total alignments in summary: 2514
STAGE 1 (mismatches=0): 2344 alignments retained (93.2%)
STAGE 2 (distance≤1000.0bp from 3' end): 1384 alignments retained (59.0%)
  Filtered out: 960 alignments too far from 3' end
STAGE 3 (unique gene mapping): 315 primers retained
  Filtered out: 20 primers with multi-gene hits near 3' end

=== FINAL RESULT: 315 primers passed all filters ===
```

### Threshold Comparison Example

Real data (FirstOrder_auto_rank2):

| Threshold | Primers Passing | Improvement vs 1000bp |
|-----------|----------------|----------------------|
| 1000bp    | 315 primers    | Baseline |
| 2000bp    | 335 primers    | +20 primers (+6.3%) |
| 5000bp    | ~350+ primers  | +35+ primers (~11%) |

**Recommendation:** Start with 1000bp, increase if you need more gene coverage.

## Strand-Aware Distance Calculation

The pipeline correctly handles strand orientation:

- **Positive strand genes (`+`):**
  - `distance_to_end = transcript_length - primer_start`
  - Measures from primer to 3' end of transcript

- **Negative strand genes (`-`):**
  - `distance_to_end = primer_start`
  - Measures from primer to 3' end (transcript start in genomic coordinates)

This ensures the filter works correctly regardless of gene strand.

## When to Use This Feature

### ✅ **Use distance filtering when:**
- Using **3' RNA-seq** protocols (QuantSeq, Lexogen, 3'READS, etc.)
- You want to **rescue primers** with irrelevant off-target hits
- You need **better gene coverage** in your primer set

### ❌ **Don't use (or use large threshold) when:**
- Using **full-length RNA-seq** or whole-transcript amplification
- Using **random priming** instead of oligo-dT
- You require **ultra-strict** specificity regardless of protocol

## Technical Details

### Handling Missing Values
- Alignments with `distance_to_end = "NA"` (unknown transcript length) are treated as `infinity`
- These alignments are **excluded** (conservative approach)

### Multi-Gene Hits
Primers are **rejected** if they have perfect alignments to multiple genes **within the distance threshold**:

```
Primer X alignments:
  - Gene A: 250bp from 3' end ✓ (within threshold)
  - Gene B: 300bp from 3' end ✓ (within threshold)
  → REJECTED: Multi-gene hit near 3' end

Primer Y alignments:
  - Gene A: 250bp from 3' end ✓ (within threshold)
  - Gene C: 4000bp from 3' end ✗ (beyond threshold, ignored)
  → ACCEPTED: Only one gene within relevant region
```

## Validation

The feature has been tested with real data showing:
- ✅ Correct filtering logic (Stage 1 → 2 → 3)
- ✅ Proper distance calculations (all values ≤ threshold)
- ✅ Strand-aware handling
- ✅ Meaningful statistics output
- ✅ 6-11% increase in primers passing with moderate threshold increase

## See Also

- `bin/select_best_primer.py` - Implementation
- `modules/SELECT_BEST_PRIMERS.nf` - Nextflow process
- `params.config` - Default parameter values
