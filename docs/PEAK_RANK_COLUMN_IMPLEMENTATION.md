# Peak Rank Column Implementation

## Summary

Added a `peak_rank` column to `selected_peaks.tsv` output to track the actual rank of the selected peak for each gene. This helps identify cases where the same peak is selected for different `peak_rank` parameter values.

## Problem

When running the pipeline with different `peak_rank` values (e.g., `peak_rank=1` vs `peak_rank=2`), some genes show the same peak selected for both ranks. This occurs when:

1. MACS2 calls only one peak for a gene
2. The requested rank (e.g., rank 2) exceeds the number of available peaks
3. The script falls back to selecting the best available peak (rank 1)

**Example**:
- Gene: `ENSG00000247596`
- MACS2 parameters: `shift=0, extsize=150`
- Only 1 peak called by MACS2
- Both `peak_rank=1` and `peak_rank=2` select the same peak

Without a rank indicator in the output, it's difficult to identify which genes had insufficient peaks for the requested rank.

## Solution

### Code Changes

Modified `bin/process_macs2_peaks.R`:

1. **Peak selection loop** (lines 241-262):
   - Added `actual_peak_rank` field to track the actual rank of the selected peak
   - When sufficient peaks exist: `actual_peak_rank = requested_rank`
   - When insufficient peaks exist: `actual_peak_rank = 1` (best available)
   - Updated warning message to indicate fallback to rank 1

2. **Output dataframe** (lines 546-580):
   - Added `peak_rank` column to `peaks_df`
   - Populated from `best_peaks$actual_peak_rank`
   - Included in empty dataframe template for consistency

### Documentation Updates

1. **docs/MANUSCRIPT.md**:
   - Added `peak_rank` to column descriptions
   - Added note explaining fallback behavior
   - Updated example output

2. **docs/TECHNICAL_APPENDIX.md**:
   - Added `peak_rank` column to format specification
   - Added explanatory note

## Behavior

### Before Changes

```r
# When peak_rank=2 but only 1 peak exists:
Warning: Gene ENSG00000247596 has only 1 peak(s), cannot select rank 2
# Gene is excluded from output (genes_without_rank)
```

### After Changes

```r
# When peak_rank=2 but only 1 peak exists:
Warning: Gene ENSG00000247596 has only 1 peak(s), cannot select rank 2 - using rank 1 instead
# Gene is included in output with peak_rank=1
```

### Output Format

**Old** `selected_peaks.tsv`:
```
gene	chr	start	end	strand	peak_score	peak_pvalue	peak_qvalue	exonic_fraction	trimmed_to_exon	strategy
ENSG00000067191	chr1	1000000	1000120	+	45.2	1.5e-05	2.3e-04	0.95	FALSE	macs2_peak_boundaries_by_score
```

**New** `selected_peaks.tsv`:
```
gene	chr	start	end	strand	summit_pos	peak_score	peak_pvalue	peak_qvalue	exonic_fraction	trimmed_to_exon	peak_rank	strategy
ENSG00000067191	chr1	1000000	1000120	+	1000060	45.2	1.5e-05	2.3e-04	0.95	FALSE	1	macs2_peak_boundaries_by_score
ENSG00000247596	chr2	2000000	2000150	-	2000075	32.1	3.2e-04	1.1e-03	0.87	FALSE	1	macs2_peak_boundaries_by_score
```

If you run with `peak_rank=2`, the second gene will still have `peak_rank=1` in the output, making it clear that the requested rank was not available.

## Downstream Impact Assessment

Analyzed all scripts that read `selected_peaks.tsv`:

### ✅ No Breaking Changes

All downstream scripts use named column access (not positional), so adding a column is safe:

1. **bin/extract_cdna_primers.R**:
   - Accesses: `peaks$gene`, `peaks$strand`
   - Impact: ✅ None

2. **bin/filter_genes_with_peaks.R**:
   - Accesses: `peaks$gene`
   - Impact: ✅ None

3. **bin/MakePlots_new.R**:
   - Accesses: `peaks_df$gene`, filters by gene ID
   - Impact: ✅ None

4. **bin/summarize_pipeline_results.py**:
   - Accesses: `peaks_df['gene']`
   - Impact: ✅ None

5. **preprocess_for_standalone.R**:
   - Reads entire table without specific column access
   - Impact: ✅ None (extra column ignored)

## Testing

To verify the fix works correctly:

1. Run pipeline with `peak_rank=1`:
   ```bash
   nextflow run main.nf --peak_rank 1 --outdir results_rank1
   ```

2. Run pipeline with `peak_rank=2`:
   ```bash
   nextflow run main.nf --peak_rank 2 --outdir results_rank2
   ```

3. Compare outputs:
   ```bash
   # Check for genes with peak_rank=1 in the rank2 output
   awk '$13 == 1' results_rank2/processed_peaks/selected_peaks.tsv
   ```

4. Verify MACS2 narrowPeak files:
   ```bash
   # For each gene with peak_rank=1 in rank2 output,
   # confirm only 1 peak was called
   grep ENSG00000247596 results_rank2/macs2_peaks/*_peaks.narrowPeak | wc -l
   ```

## Example Usage

Identify genes where rank 2 wasn't available:

```bash
# Extract genes that got rank 1 when rank 2 was requested
awk -F'\t' 'NR==1 || $13==1' results_rank2/processed_peaks/selected_peaks.tsv > genes_with_insufficient_peaks.tsv
```

Compare peak selections between different ranks:

```bash
# Join outputs to compare
join -t $'\t' -1 1 -2 1 \
  <(tail -n+2 results_rank1/processed_peaks/selected_peaks.tsv | sort -k1,1) \
  <(tail -n+2 results_rank2/processed_peaks/selected_peaks.tsv | sort -k1,1) \
  | awk -F'\t' '$13==$26' > same_peak_selected.tsv
```

## Benefits

1. **Transparency**: Users can immediately see which peaks were actually selected
2. **Quality Control**: Easy to identify genes with limited peak calls
3. **Debugging**: Helps diagnose why different rank parameters produce similar results
4. **Reproducibility**: Output files now contain complete information about selection decisions

## Files Modified

- `bin/process_macs2_peaks.R`: Core logic changes
- `docs/MANUSCRIPT.md`: Documentation update
- `docs/TECHNICAL_APPENDIX.md`: Format specification update

## Backward Compatibility

⚠️ **Minor breaking change**: Output format has one additional column

- **Impact**: Low - all downstream scripts use named column access
- **Migration**: None required for existing pipelines
- **Old data**: Can still be processed (scripts will ignore missing column)
