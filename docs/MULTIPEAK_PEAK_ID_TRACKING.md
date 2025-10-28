# Multi-Peak Mode: Peak Identifier Tracking Implementation

## Overview

This document describes the implementation of peak identifier tracking throughout the primer design pipeline in multi-peak mode. This ensures that primers can be traced back to the specific genomic peaks they were designed from.

## Problem Statement

**Original Issue:**
In multi-peak mode, when a gene has multiple MACS2 peaks (e.g., ENSG00000107130 with 7 peaks), the pipeline:
1. ✅ Correctly generated `selected_peaks.fa` with peak identifiers: `>ENSG00000107130|peak_1`, `>ENSG00000107130|peak_2`, etc.
2. ❌ Lost peak identifiers during Primer3 input conversion, using only `SEQUENCE_ID=ENSG00000107130`
3. ❌ Could not link primers back to specific peaks in downstream analysis

**Impact on Analysis:**
- **Primer-Peak Linkage Lost**: Cannot determine which primers correspond to which genomic peaks
- **Isoform Optimization Broken**: Multi-peak optimization cannot properly assign primers to peaks for coverage analysis
- **Quality Control Impaired**: Cannot trace failed primers back to problematic peaks
- **Alignment Ambiguity**: Transcriptome alignment results cannot be mapped to genomic coordinates
- **Plotting Failures**: Cannot visualize primer locations relative to correct peak regions

## Solution Implementation

### 1. FASTA Header Format (`process_macs2_peaks.R`)

**Location:** Lines 676-692

**Multi-peak mode FASTA headers:**
```
>ENSG00000107130|peak_1 chr9:130236156-130236198(+)
>ENSG00000107130|peak_2 chr9:130237136-130237190(+)
```

**Format:** `gene_id|peak_N genomic_coordinates`

**Code:**
```r
names(sequences) <- paste0(
  target_regions$gene_id, "|peak_", best_peaks$actual_peak_rank, " ",
  seqnames(target_regions), ":",
  start(target_regions), "-",
  end(target_regions), "(",
  strand(target_regions), ")"
)
```

### 2. Primer3 Input Generation (`fasta_to_primer3.py`)

**Location:** Lines 43-68

**MODIFIED:** Preserves peak identifiers in `SEQUENCE_ID` field

**Before:**
```
SEQUENCE_ID=ENSG00000107130  # Peak identifier lost!
```

**After:**
```
SEQUENCE_ID=ENSG00000107130|peak_1  # Peak identifier preserved
SEQUENCE_ID=ENSG00000107130|peak_2
```

**Implementation:**
```python
if '|' in header:
    parts = header.split('|')
    gene_id = parts[0]
    
    # Check if this is multi-peak format (has peak_N identifier)
    if len(parts) > 1 and parts[1].strip().split()[0].startswith('peak_'):
        peak_info = parts[1].strip().split()[0]  # Get "peak_1"
        sequence_id = f"{gene_id}|{peak_info}"
    else:
        sequence_id = gene_id
```

### 3. Primer Extraction (`extract_cdna_primers.R`)

**Location:** Lines 103-126

**MODIFIED:** Added `peak_id` column to output table

**Output columns:**
- `gene_id`: Gene identifier (e.g., "ENSG00000107130")
- `peak_id`: Peak identifier (e.g., "peak_1") or NA for single-peak mode
- `primer_index`: Primer number from Primer3
- `primer_type`: LEFT/RIGHT
- `primer_sequence`: Primer sequence
- `gene_strand`: +/-
- `rationale`: Why this primer type was selected

**Example output (`cdna_primers.tsv`):**
```
gene_id              peak_id    primer_index    primer_type    primer_sequence
ENSG00000107130      peak_1     0               LEFT           CAGTGGTGAGTAGCTGCTTT
ENSG00000107130      peak_1     1               LEFT           CAGTGGTGAGTAGCTGCTTTT
ENSG00000107130      peak_2     0               LEFT           TTGGGGTGCCGTCCTGTC
```

### 4. Alignment FASTA Headers (`primers_to_fasta.R`)

**Location:** Lines 36-52

**MODIFIED:** Includes peak_id in FASTA headers for transcriptome alignment

**Multi-peak format:**
```
>ENSG00000107130|peak_1|idx0|LEFT|strand+
>ENSG00000107130|peak_1|idx1|LEFT|strand+
>ENSG00000107130|peak_2|idx0|LEFT|strand+
```

**Single-peak format (unchanged):**
```
>ENSG00000107130|idx0|LEFT|strand+
```

**Implementation:**
```r
if ("peak_id" %in% colnames(primers) && any(!is.na(primers$peak_id))) {
  # Multi-peak mode
  headers <- paste0(
    primers$gene_id, "|",
    primers$peak_id, "|",
    "idx", primers$primer_index, "|",
    primers$primer_type, "|",
    "strand", primers$gene_strand
  )
} else {
  # Single-peak mode
  headers <- paste0(
    primers$gene_id, "|",
    "idx", primers$primer_index, "|",
    primers$primer_type, "|",
    "strand", primers$gene_strand
  )
}
```

### 5. Alignment Analysis (`analyze_primer_alignments.py`)

**Location:** Lines 92-106

**MODIFIED:** Constructs `primer_id` with peak identifier when present

**Multi-peak primer_id:**
```
ENSG00000107130|peak_1|idx0|LEFT|strand+
```

**Implementation:**
```python
if 'peak_id' in primers.columns and primers['peak_id'].notna().any():
    # Multi-peak mode: include peak_id in primer_id
    primers['primer_id'] = (primers['gene_id'].astype(str) + "|" +
                           primers['peak_id'].astype(str) + "|" +
                           "idx" + primers['primer_index'].astype(str) + "|" +
                           primers['primer_type'].astype(str) + "|" +
                           "strand" + primers['gene_strand'].astype(str))
```

### 6. Best Primer Selection (`select_best_primer.py`)

**Location:** Lines 43-77

**MODIFIED:** Handles peak_id column in primer_id construction

**Effect:** Ensures multi-gene filtering and distance filtering work correctly per peak

## Data Flow Summary

```
┌─────────────────────────────────────────────────────────────────────┐
│ 1. MACS2 Peak Detection                                              │
│    Gene ENSG00000107130 has 7 peaks at different genomic locations  │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 2. Peak Processing (process_macs2_peaks.R)                          │
│    selected_peaks.fa:                                                │
│      >ENSG00000107130|peak_1 chr9:130236156-130236198(+)           │
│      >ENSG00000107130|peak_2 chr9:130237136-130237190(+)           │
│    selected_peaks.tsv:                                               │
│      gene=ENSG00000107130 peak_rank=1 chr=9 start=130236156         │
│      gene=ENSG00000107130 peak_rank=2 chr=9 start=130237136         │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 3. Primer3 Input (fasta_to_primer3.py) **CRITICAL STEP**            │
│    primer3_input.txt:                                                │
│      SEQUENCE_ID=ENSG00000107130|peak_1  ← Peak ID preserved!       │
│      SEQUENCE_TEMPLATE=TTTCGGAGGGGGTTGGTGGGGAGGTC...                │
│      =                                                                │
│      SEQUENCE_ID=ENSG00000107130|peak_2                             │
│      SEQUENCE_TEMPLATE=TTGGGGTGCCGTCCTGTCTGAACCTG...                │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 4. Primer3 Design                                                    │
│    primer3_output.txt:                                               │
│      SEQUENCE_ID=ENSG00000107130|peak_1  ← Peak ID in output        │
│      PRIMER_LEFT_0_SEQUENCE=CAGTGGTGAGTAGCTGCTTT                    │
│      PRIMER_LEFT_1_SEQUENCE=CAGTGGTGAGTAGCTGCTTTT                   │
│      =                                                                │
│      SEQUENCE_ID=ENSG00000107130|peak_2                             │
│      PRIMER_LEFT_0_SEQUENCE=TTGGGGTGCCGTCCTGTC                      │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 5. Primer Extraction (extract_cdna_primers.R)                       │
│    cdna_primers.tsv:                                                 │
│      gene_id          peak_id  primer_index  primer_sequence        │
│      ENSG00000107130  peak_1   0             CAGTGGTGAGTAGCTGCTTT   │
│      ENSG00000107130  peak_1   1             CAGTGGTGAGTAGCTGCTTTT  │
│      ENSG00000107130  peak_2   0             TTGGGGTGCCGTCCTGTC     │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 6. Alignment FASTA (primers_to_fasta.R)                             │
│    primers_for_alignment.fa:                                         │
│      >ENSG00000107130|peak_1|idx0|LEFT|strand+                      │
│      >ENSG00000107130|peak_1|idx1|LEFT|strand+                      │
│      >ENSG00000107130|peak_2|idx0|LEFT|strand+                      │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 7. Transcriptome Alignment (bowtie2)                                │
│    primers_alignment.bam: Headers preserved in read names           │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 8. Alignment Analysis (analyze_primer_alignments.py)                │
│    primer_alignment_summary.tsv:                                     │
│      gene_id          peak_id  primer_index  aligned_gene_name      │
│      ENSG00000107130  peak_1   0             FOXE1                  │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 9. Best Primer Selection (select_best_primer.py)                    │
│    best_primers.tsv: Filtered primers with peak tracking            │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│ 10. Multi-Peak Optimization (optimize_primers_multipeak.py)         │
│     optimized_primers_multipeak.tsv: Peak-aware isoform selection   │
└─────────────────────────────────────────────────────────────────────┘
```

## Testing

### Test Case: ENSG00000107130

**Input:**
- 7 MACS2 peaks detected
- Multi-peak mode: `--select_all_peaks --optimize_multipeak`

**Expected Output:**

1. **selected_peaks.fa** should have 7 sequences with headers:
   ```
   >ENSG00000107130|peak_1 chr9:...
   >ENSG00000107130|peak_2 chr9:...
   ...
   ```

2. **primer3_input.txt** should have 7 blocks with:
   ```
   SEQUENCE_ID=ENSG00000107130|peak_1
   SEQUENCE_ID=ENSG00000107130|peak_2
   ...
   ```

3. **cdna_primers.tsv** should have `peak_id` column:
   ```
   gene_id              peak_id    primer_index
   ENSG00000107130      peak_1     0
   ENSG00000107130      peak_1     1
   ENSG00000107130      peak_2     0
   ```

4. **primers_for_alignment.fa** headers should include peak_id:
   ```
   >ENSG00000107130|peak_1|idx0|LEFT|strand+
   ```

5. **best_primers.tsv** should retain peak_id information

### Verification Commands

```bash
# Check peak identifiers in FASTA
grep ">" results/multipeak_plot/processed_peaks/selected_peaks.fa | head

# Check Primer3 input preserves peak IDs
grep "SEQUENCE_ID=" results/multipeak_plot/primer3_input.txt | head

# Check cdna_primers.tsv has peak_id column
head results/multipeak_plot/cdna_primers.tsv

# Check alignment FASTA headers
grep ">" results/multipeak_plot/primers_for_alignment.fa | head

# Check best primers retain peak tracking
head results/multipeak_plot/best_primers.tsv
```

## Backward Compatibility

**Single-Peak Mode:**
All scripts check for presence of `peak_id` column and handle both modes:
- If `peak_id` present and has non-NA values → Multi-peak mode
- Otherwise → Single-peak mode (original behavior)

**Single-peak FASTA headers (unchanged):**
```
>ENSG00000123456 chr1:1000-1300(+)
```

**Single-peak SEQUENCE_ID (unchanged):**
```
SEQUENCE_ID=ENSG00000123456
```

## Files Modified

1. **bin/fasta_to_primer3.py** - Lines 43-68: Preserve peak_N in SEQUENCE_ID
2. **bin/extract_cdna_primers.R** - Lines 103-126: Add peak_id column parsing
3. **bin/primers_to_fasta.R** - Lines 36-52: Include peak_id in alignment headers
4. **bin/analyze_primer_alignments.py** - Lines 92-106: Handle peak_id in primer_id
5. **bin/select_best_primer.py** - Lines 43-77: Support peak_id in filtering

## Benefits

1. ✅ **Complete Traceability**: Every primer can be traced to its genomic peak
2. ✅ **Accurate Isoform Analysis**: Optimization scripts can correctly assign primers to peaks
3. ✅ **Better Quality Control**: Failed primers can be linked to problematic peaks
4. ✅ **Proper Visualization**: Plots can show primers at correct genomic locations
5. ✅ **Downstream Analysis**: External tools can map primers back to genome
6. ✅ **Backward Compatible**: Single-peak mode behavior unchanged

## Date Implemented

October 27, 2025
