# Isoform-Optimized Primer Selection

## Overview

The `optimize_primer_isoforms.py` script optimizes primer selection to maximize the number of distinct transcript isoforms targeted while ensuring no isoform is targeted by multiple primers.

## Problem Statement

### Current Situation

The `best_primers.tsv` file contains primers that meet quality requirements:
- ✓ Perfect alignment to target gene (0 mismatches)
- ✓ Within distance threshold from 3' end
- ✓ Unique gene mapping (no off-target hits near 3' end)

However, multiple primers may be available for each gene, and:
1. **Within-gene redundancy**: Different primers for the same gene may target overlapping sets of isoforms
2. **Cross-gene redundancy**: Primers for different genes may target the same isoforms (if genes share isoforms)

### Goal

Select **one optimal primer per gene** that:
- Maximizes total distinct isoforms covered
- Eliminates isoform redundancy between primers
- Maintains one-to-one gene-to-primer relationship

## Algorithm

### Phase 1: Within-Gene Optimization

For each gene with multiple primers:
1. Count distinct isoforms targeted by each primer
2. Select primer with maximum isoform coverage
3. Tie-breaking (in order):
   - Minimum distance to 3' end
   - Fewer total alignments (more specific)
   - Lexicographically first primer_id

### Phase 2: Cross-Gene Deduplication

Iteratively resolve conflicts where multiple primers target the same isoform:

```
while conflicts exist:
    for each conflicting isoform:
        identify all primers targeting it
        score primers by:
            - total isoforms targeted (higher is better)
            - number of alternative primers available (fewer is higher priority)
            - distance to 3' end (closer is better)
        keep primer with best score
        remove others
        for genes that lost their primer:
            try to select alternative primer without conflicts
    update conflict list
```

### Result

- One primer per gene
- No isoform targeted by multiple primers
- Maximum distinct isoform coverage

## Usage

### Basic Usage

```bash
python bin/optimize_primer_isoforms.py \
  --best_primers results/best_primers.tsv \
  --alignment_summary results/primer_alignment_summary.tsv \
  --output results/best_primers_optimal.tsv
```

### With Custom Distance Threshold

```bash
python bin/optimize_primer_isoforms.py \
  --best_primers results/best_primers.tsv \
  --alignment_summary results/primer_alignment_summary.tsv \
  --output results/best_primers_optimal.tsv \
  --distance_threshold 1500
```

### Parameters

- `--best_primers`: Input file with validated primers (required)
- `--alignment_summary`: Input file with transcript alignment details (required)
- `--output`: Output file for optimized primers (default: `best_primers_optimal.tsv`)
- `--distance_threshold`: Maximum distance from 3' end in bp (default: 1000)

## Input Files

### best_primers.tsv

Required columns:
- `gene_id`: Target gene identifier
- `primer_id`: Unique primer identifier
- `distance_to_end_min`: Distance to 3' end (optional, for tie-breaking)
- `zero_mismatch_alignments`: Number of perfect alignments (optional, for tie-breaking)

### primer_alignment_summary.tsv

Required columns:
- `primer_id`: Links to best_primers.tsv
- `aligned_transcript`: Transcript/isoform ID
- `aligned_gene_name`: Gene that transcript belongs to
- `mismatches`: Number of mismatches (must be 0 for inclusion)
- `distance_to_end`: Distance from 3' end (optional, filtered by threshold)

## Output File

### best_primers_optimal.tsv

Contains all columns from `best_primers.tsv` plus:

- `isoforms_targeted`: Number of distinct isoforms this primer targets
- `target_transcripts`: Comma-separated list of transcript IDs

### Example Output

```
gene_id          primer_id              isoforms_targeted  target_transcripts
ENSG00000001234  ENSG00000001234|idx0   3                  ENST00000001,ENST00000002,ENST00000003
ENSG00000005678  ENSG00000005678|idx1   2                  ENST00000010,ENST00000011
```

## Example Scenarios

### Scenario 1: Multiple Primers Per Gene

**Input:**
- Gene A: Primer A1 (targets isoforms 1, 2, 3), Primer A2 (targets isoforms 2, 3)
- Gene B: Primer B1 (targets isoforms 10, 11)

**Result:**
- Select A1 (3 isoforms > 2 isoforms)
- Select B1 (only option)
- Total: 5 distinct isoforms

### Scenario 2: Cross-Gene Isoform Conflict

**Input:**
- Gene A: Primer A1 (targets isoforms 1, 2, 3)
- Gene B: Primer B1 (targets isoforms 3, 4, 5)  ← Conflict on isoform 3!

**Resolution:**
- Both primers target 3 isoforms (tie)
- Check alternatives: both genes have only 1 primer (tie)
- Compare distance: A1=100bp, B1=150bp
- Keep A1 (closer to 3' end)
- Gene B has no alternative → removed from output

**Result:**
- Gene A: Primer A1 (isoforms 1, 2, 3)
- Gene B: No primer
- Total: 3 distinct isoforms (no conflicts)

### Scenario 3: Complex Multi-Way Conflict

**Input:**
- Gene A: Primer A1 (isoforms 1, 2), Primer A2 (isoform 3)
- Gene B: Primer B1 (isoforms 2, 4), Primer B2 (isoform 5)
- Gene C: Primer C1 (isoforms 2, 6)  ← All conflict on isoform 2!

**Initial Selection:**
- Gene A: A1 (2 isoforms)
- Gene B: B1 (2 isoforms)
- Gene C: C1 (2 isoforms)

**Conflict Resolution (isoform 2):**
- A1, B1, C1 all target isoform 2
- All have 2 isoforms (tie)
- Check alternatives: A has A2, B has B2, C has only C1
- C has fewer alternatives → C1 wins
- Remove A1 and B1

**Re-selection:**
- Gene A: Try A2 (isoform 3) → No conflict! ✓
- Gene B: Try B2 (isoform 5) → No conflict! ✓

**Final Result:**
- Gene A: A2 (isoform 3)
- Gene B: B2 (isoform 5)
- Gene C: C1 (isoforms 2, 6)
- Total: 4 distinct isoforms (3, 5, 2, 6)

## Performance

### Expected Performance

- **Input size**: 100-1000 genes, 1-10 primers per gene
- **Runtime**: < 1 minute for typical datasets
- **Memory**: < 100 MB

### Complexity

- **Time**: O(G × P × I × C) where:
  - G = number of genes
  - P = primers per gene
  - I = isoforms per primer
  - C = conflict resolution iterations
- **Space**: O(G × P + T) where T = total isoforms

## Validation

### Output Verification

```bash
# Count primers per gene (should be 1)
awk -F'\t' 'NR>1 {print $1}' best_primers_optimal.tsv | sort | uniq -c | sort -rn

# Check for duplicate isoforms (should be empty)
awk -F'\t' 'NR>1 {print $NF}' best_primers_optimal.tsv | \
  tr ',' '\n' | sort | uniq -d

# Compare coverage
echo "Before optimization:"
wc -l < best_primers.tsv
echo "After optimization:"
wc -l < best_primers_optimal.tsv
```

### Quality Checks

1. **One primer per gene**: `awk 'NR>1{print $1}' output.tsv | sort | uniq -d` should be empty
2. **No isoform conflicts**: All transcript IDs should be unique across all primers
3. **Coverage comparison**: Total distinct isoforms should be maximized

## Integration with Pipeline

### Option 1: Add to Nextflow Pipeline

Create a new module:

```groovy
// modules/OPTIMIZE_PRIMER_ISOFORMS.nf
process OPTIMIZE_PRIMER_ISOFORMS {
  tag 'optimize_isoforms'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.10 pandas=2.0"

  input:
  path best_primers
  path alignment_summary

  output:
  path 'best_primers_optimal.tsv'

  script:
  """
  optimize_primer_isoforms.py \\
    --best_primers ${best_primers} \\
    --alignment_summary ${alignment_summary} \\
    --output best_primers_optimal.tsv \\
    --distance_threshold ${params.distance_threshold ?: 1000}
  """
}
```

Add to workflow:

```groovy
// After SELECT_BEST_PRIMERS
optimized_primers = OPTIMIZE_PRIMER_ISOFORMS(
  best_primers,
  primer_alignment_summary
)
```

### Option 2: Post-Processing Script

Run after pipeline completes:

```bash
# After pipeline finishes
python bin/optimize_primer_isoforms.py \
  --best_primers results/best_primers.tsv \
  --alignment_summary results/primer_alignment_summary.tsv \
  --output results/best_primers_optimal.tsv
```

## Limitations and Edge Cases

### Limitations

1. **Requires transcript-level alignment data**: Cannot optimize without `primer_alignment_summary.tsv`
2. **Greedy algorithm**: May not find global optimum (but provides good solution quickly)
3. **Gene priority**: Genes with fewer alternative primers get priority in conflicts

### Edge Cases

| Case | Behavior |
|------|----------|
| No alignment summary | Exit with error |
| No isoform data | Copy input to output (no optimization) |
| Single primer per gene | Automatic selection (no within-gene optimization) |
| Gene loses all primers | Excluded from output (reported in summary) |
| All primers conflict | Iterative removal until no conflicts remain |

## Troubleshooting

### Issue: "No isoform data available"

**Cause**: Alignment summary missing or empty

**Solution**: Verify `primer_alignment_summary.tsv` exists and contains transcript IDs

### Issue: "Many genes removed"

**Cause**: High cross-gene isoform overlap

**Solutions**:
- Increase `--distance_threshold` to be more permissive
- Check if genes actually share transcripts (expected for gene families)
- Review input primer quality

### Issue: "Could not resolve conflicts"

**Cause**: Circular dependencies or insufficient alternatives

**Solution**: Check iteration logs, may need manual review of conflicting primers

## References

- Related scripts:
  - `bin/select_best_primer.py`: Initial 3-stage filtering
  - `bin/analyze_primer_alignments.py`: Generate alignment summary
  
- Documentation:
  - `docs/ISOFORM_OPTIMIZATION_PLAN.md`: Detailed algorithm design
