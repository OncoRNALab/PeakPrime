# Isoform Optimization Integration

## Overview
The isoform optimization module has been successfully integrated into both primer design workflows. It runs automatically after primer selection to maximize distinct isoform coverage.

## Integration Points

### Peak-Based Workflow (`primer_design.nf`)
- **Module Call**: After `SELECT_BEST_PRIMERS`
- **Input Files**:
  - `best_primers.tsv` from `SELECT_BEST_PRIMERS`
  - `primer_alignment_summary.tsv` from `ANALYZE_PRIMER_ALIGNMENTS`
- **Output File**: `best_primers_optimal.tsv`
- **Location**: Line ~110 in workflow file

### Distance-Based Workflow (`distance_primer_design.nf`)
- **Module Call**: After `SELECT_BEST_PRIMERS`
- **Input Files**:
  - `best_primers.tsv` from `SELECT_BEST_PRIMERS`
  - `primer_alignment_summary.tsv` from `ANALYZE_PRIMER_ALIGNMENTS`
- **Output File**: `best_primers_optimal.tsv`
- **Location**: Line ~108 in workflow file
- **Note**: Only runs when `--transcriptome_index` is provided

## Workflow Changes

### Files Modified
1. **workflows/primer_design.nf**
   - Added include statement for `OPTIMIZE_PRIMER_ISOFORMS`
   - Added module call after `SELECT_BEST_PRIMERS`
   - Added `optimized_primers` to workflow emit section

2. **workflows/distance_primer_design.nf**
   - Added include statement for `OPTIMIZE_PRIMER_ISOFORMS`
   - Added module call after `SELECT_BEST_PRIMERS`
   - Added `optimized_primers` to workflow emit section
   - Set `optimized_primers = Channel.empty()` when alignment is skipped

3. **modules/OPTIMIZE_PRIMER_ISOFORMS.nf** (created)
   - New module that wraps `bin/optimize_primer_isoforms.py`
   - Uses conda environment: Python 3.10, pandas 2.0
   - Published to `params.outdir`

## Pipeline Flow

### With Transcriptome Alignment (both workflows)
```
SELECT_BEST_PRIMERS
    ↓ best_primers.tsv
    ↓
OPTIMIZE_PRIMER_ISOFORMS ← primer_alignment_summary.tsv
    ↓
best_primers_optimal.tsv (published)
```

### Without Transcriptome Alignment (distance workflow only)
```
Best primers not available
    ↓
optimized_primers = Channel.empty()
```

## Output Files

### best_primers.tsv
- Original best primer selection
- Contains all primers that passed quality filters
- May include primers targeting same isoforms

### best_primers_optimal.tsv
- Isoform-optimized primer selection
- Maximizes distinct isoform coverage
- Eliminates redundant primers targeting same isoforms
- Subset of primers from `best_primers.tsv`

## Usage

No additional parameters required. The optimization runs automatically when:
- Running peak-based workflow with transcriptome alignment
- Running distance-based workflow with `--transcriptome_index` provided

The `--distance_threshold` parameter (if set) is automatically passed to the optimization module.

## Testing

To test the integration:

```bash
# Peak-based workflow
nextflow run main.nf \
  --bam sample.bam \
  --gtf genes.gtf \
  --genome genome.fa \
  --transcriptome_fasta transcripts.fa \
  --transcriptome_index transcripts

# Distance-based workflow
nextflow run main.nf \
  --distance_mode \
  --genes genes.txt \
  --transcriptome_fasta transcripts.fa \
  --transcriptome_index transcripts
```

Check for `best_primers_optimal.tsv` in the output directory.

## Technical Details

### Algorithm
1. **Phase 1**: Select best primer per gene (maximum isoform coverage)
2. **Phase 2**: Resolve inter-gene conflicts (remove redundant primers)
3. **Iterative Refinement**: Continue until no primers share isoforms

### Key Design Decisions
- Uses composite keys `(gene_id, primer_index)` for unique primer identification
- Matches primers using `gene_id` from both input files (not gene symbols)
- Handles both Ensembl gene IDs and gene symbols automatically
- Preserves all columns from original `best_primers.tsv`

### Dependencies
- Python 3.10+
- pandas 2.0+
- Configured via conda environment in module

## Related Documentation
- [Isoform Optimization Plan](ISOFORM_OPTIMIZATION_PLAN.md)
- [Isoform Optimization Usage](ISOFORM_OPTIMIZATION_USAGE.md)
- [Isoform Optimization Summary](ISOFORM_OPTIMIZATION_SUMMARY.md)
