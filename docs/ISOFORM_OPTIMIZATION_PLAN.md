# Isoform-Optimized Primer Selection Plan

## Overview

Create an optimized primer selection algorithm that maximizes the number of distinct isoforms (transcripts) targeted across all genes, avoiding redundant isoform coverage.

## Current State

- **Input**: `best_primers.tsv` - primers that meet quality requirements (distance, mismatches, uniqueness)
- **Limitation**: May have multiple primers per gene targeting overlapping isoforms
- **Goal**: Select one optimal primer per gene that maximizes total unique isoform coverage

## Algorithm Design

### Step-by-Step Plan

#### Phase 1: Data Collection and Preparation

1. **Load Input Data**
   - Read `best_primers.tsv` (validated primers)
   - Read `primer_alignment_summary.tsv` (detailed alignment data with transcript IDs)
   
2. **Extract Isoform Information**
   - For each primer, collect all aligned transcripts (isoforms)
   - Filter for:
     - Perfect matches (mismatches = 0)
     - Within distance threshold (distance_to_end ≤ threshold)
     - Aligned to target gene (aligned_gene_name matches gene_id)
   
3. **Build Primer-Isoform Mapping**
   - Create data structure: `{primer_id: set(transcript_ids)}`
   - Count isoforms per primer: `{primer_id: isoform_count}`
   - Group primers by gene: `{gene_id: [primer_ids]}`

#### Phase 2: Within-Gene Optimization

4. **Select Best Primer Per Gene**
   - For each gene with multiple primers:
     - Rank primers by number of distinct isoforms targeted (descending)
     - Select primer with maximum isoform coverage
     - Tie-breaking criteria (in order):
       1. Minimum distance to 3' end (closer is better)
       2. Fewer total alignments (more specific)
       3. Lexicographically first primer_id (deterministic)

#### Phase 3: Cross-Gene Deduplication

5. **Build Global Isoform Assignment**
   - Create initial set of selected primers (one per gene)
   - Track which isoforms are targeted by each selected primer
   
6. **Detect Isoform Conflicts**
   - Identify isoforms targeted by multiple primers
   - For each conflicting isoform, note which primers target it
   
7. **Iterative Conflict Resolution**
   - While conflicts exist:
     - For each conflicting isoform:
       - Find all primers targeting it
       - Keep primer with highest total isoform count
       - For tied primers, use tie-breaking criteria:
         1. Gene with fewer alternative primers available
         2. Minimum distance to 3' end
         3. Lexicographically first primer_id
     - Remove losing primers from selection
     - For genes that lost their primer:
       - Select next-best available primer
       - Re-check for conflicts
     - Update conflict list
   
8. **Verify Solution**
   - Ensure each selected primer targets only unique isoforms
   - Count total distinct isoforms covered
   - Verify one primer per gene (or zero if no valid primer)

#### Phase 4: Output Generation

9. **Create Output File**
   - Write `best_primers_optimal.tsv` with columns:
     - All original columns from `best_primers.tsv`
     - `isoforms_targeted`: Number of distinct isoforms
     - `target_transcripts`: Comma-separated list of transcript IDs
     - `optimization_rank`: Selection order/priority
   
10. **Generate Summary Report**
    - Total genes processed
    - Total primers in input
    - Total primers in output (one per gene)
    - Total distinct isoforms covered
    - Genes with no primers (if any)
    - Conflict resolution statistics

## Data Structures

```python
# Core data structures
primer_isoforms: Dict[str, Set[str]]          # primer_id -> {transcript_ids}
gene_primers: Dict[str, List[str]]             # gene_id -> [primer_ids]
selected_primers: Dict[str, str]               # gene_id -> primer_id
isoform_ownership: Dict[str, Set[str]]         # transcript_id -> {primer_ids}
```

## Edge Cases

1. **No isoform data available**
   - Fallback to current best_primers.tsv (no optimization)
   
2. **Primer aligns to multiple genes**
   - Already filtered out by best_primers selection
   
3. **Gene has only one primer**
   - Automatically selected (no within-gene optimization needed)
   
4. **Gene has no primers after conflict resolution**
   - Report in summary, exclude from output
   
5. **Isoform has no length/distance information**
   - Skip distance-based filtering, use only perfect matches
   
6. **Tie in isoform count and all criteria**
   - Use deterministic lexicographic ordering of primer_id

## Performance Considerations

- **Expected input size**: 100-1000 genes × 1-10 primers/gene
- **Complexity**: O(G × P × I) where G=genes, P=primers/gene, I=isoforms/primer
- **Optimization**: Use sets for isoform operations (fast intersection/difference)

## Testing Strategy

1. **Unit tests**:
   - Single gene with multiple primers
   - Multiple genes with overlapping isoforms
   - Edge cases (no data, single primer, etc.)

2. **Integration test**:
   - Run on actual pipeline output
   - Verify output format matches input
   - Confirm isoform uniqueness

3. **Validation**:
   - Compare total isoform coverage: optimized vs. naive selection
   - Verify no isoform is targeted by multiple primers

## Success Metrics

- **Coverage**: Maximized number of distinct isoforms
- **Specificity**: Each isoform targeted by at most one primer
- **Completeness**: One primer per gene (where possible)
- **Determinism**: Same input always produces same output
