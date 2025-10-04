# Optimal Primer Set Selection for Maximum Isoform Coverage

## Overview

This tool selects the optimal combination of primers that maximizes the number of amplified isoforms for each target gene. It implements a greedy set-cover algorithm with quality filters.

## Algorithm

The selection process follows these steps:

### 1. Quality Filtering

**Specificity Filter:**
- Keep only primers with **perfect matches** (0 mismatches) to the target gene
- Remove primers that have perfect matches to other genes (off-targets)
- Primers matching multiple isoforms of the **same gene** are allowed

**Distance Filter:**
- Keep only primers within a threshold distance from the transcript 3' end
- Default: 400 nucleotides
- Rationale: Primers too far from the end may not amplify efficiently in RT-PCR

### 2. Greedy Set Cover Algorithm

For each gene, iteratively select primers that maximize isoform coverage:

```
1. Start with:
   - Selected primers = empty set
   - Covered isoforms = empty set

2. While (iterations < max_primers) AND (uncovered isoforms exist):
   a. For each remaining primer:
      - Calculate how many NEW isoforms it would cover
   
   b. Select the primer that covers the most NEW isoforms
   
   c. Add primer to selected set
   
   d. Add its isoforms to covered set
   
   e. Remove primer from candidates

3. Return selected primer set
```

### 3. Output

For each gene, the tool reports:
- Selected primers (ranked by selection order)
- Primer sequences and metadata
- Number of isoforms covered by each primer
- Total coverage (% of all isoforms)
- Uncovered isoforms (if any)

## Why Greedy Set Cover?

This is a classic **set cover problem** where:
- **Universe** = all isoforms of the target gene
- **Sets** = isoforms covered by each primer
- **Goal** = find minimum primers covering maximum isoforms

The greedy algorithm:
- ✅ Efficient: O(n²) time complexity
- ✅ Near-optimal: Guaranteed within ln(n) of optimal solution
- ✅ Practical: Works well for typical gene/isoform counts (5-20 isoforms)
- ✅ Interpretable: Selection order shows primer importance

## Installation

No installation required. Scripts use standard Python/R libraries:

**Python:** pandas (included in most Python distributions)
**R:** data.table (install via `install.packages("data.table")`)

## Usage

### Python Version

```bash
# Single gene
python bin/select_optimal_primer_set.py \
  -i results/my_run/primer_alignment_summary.tsv \
  -g ENSG00000104687 \
  -o optimal_primers_GSR.tsv

# Custom parameters
python bin/select_optimal_primer_set.py \
  -i results/my_run/primer_alignment_summary.tsv \
  -g ENSG00000104687 \
  -d 300 \          # Max distance: 300nt
  -m 3 \            # Max 3 primers
  -o optimal_primers_GSR.tsv

# All genes in file
python bin/select_optimal_primer_set.py \
  -i results/my_run/primer_alignment_summary.tsv \
  --all-genes \
  -o optimal_primers_all.tsv
```

### R Version

```bash
# Single gene
Rscript bin/select_optimal_primer_set.R \
  results/my_run/primer_alignment_summary.tsv \
  ENSG00000104687 \
  400 \             # Max distance
  5 \               # Max primers
  optimal_primers_GSR.tsv

# Defaults (400nt, 5 primers)
Rscript bin/select_optimal_primer_set.R \
  results/my_run/primer_alignment_summary.tsv \
  ENSG00000104687
```

## Output Format

TSV file with columns:

| Column | Description |
|--------|-------------|
| `rank` | Selection order (1 = best) |
| `primer_index` | Primer index from original design |
| `primer_sequence` | DNA sequence |
| `gene_id` | Target gene Ensembl ID |
| `gene_name` | Target gene symbol |
| `gene_strand` | Gene strand (+/-) |
| `min_distance_to_end` | Closest distance to any isoform end |
| `isoforms_covered` | Number of isoforms covered |
| `isoform_ids` | Comma-separated transcript IDs |

## Example Output

```
Processing gene: ENSG00000104687
===============================================================================

Filtering primers for gene: ENSG00000104687
  Total alignments for ENSG00000104687: 156
  Unique primers: 12
  Specific primers (no off-target): 8

Filtering by distance to end (max=400nt)
  Primers before: 8, after: 5, removed: 3

Selecting optimal primer combination (max=5 primers)
  Total candidate primers: 5
  Total isoforms to cover: 3

  Iteration 1: Selected primer 3
    - Covers 2 new isoforms: ENST00000221130, ENST00000643653
    - Total coverage: 2/3 (66.7%)
    - Distance to end: 378nt

  Iteration 2: Selected primer 7
    - Covers 1 new isoforms: ENST00000644249
    - Total coverage: 3/3 (100.0%)
    - Distance to end: 95nt

FINAL RESULTS for ENSG00000104687:
  Selected 2 primers
  Coverage: 3/3 isoforms (100.0%)
  All isoforms covered!
```

## Parameters

### Distance Threshold (`-d`, `--max-distance`)

**Default:** 400 nucleotides

**Rationale:**
- cDNA synthesis from mRNA starts at 3' poly(A) tail
- Primers closer to 3' end have higher amplification efficiency
- 400nt is a conservative threshold for most RT-PCR applications

**Adjust if:**
- Your protocol uses longer cDNA synthesis: increase (e.g., 600-800nt)
- You want only highly-expressed 3' UTRs: decrease (e.g., 200-300nt)
- Transcripts are very short: increase to avoid filtering all primers

### Maximum Primers (`-m`, `--max-primers`)

**Default:** 5 primers

**Rationale:**
- Most genes have 5-10 major isoforms
- Diminishing returns after 3-5 primers (often 80%+ coverage)
- Lab validation cost increases linearly with primer count

**Adjust if:**
- Gene has many isoforms (>10): increase to 7-10
- Validation budget is limited: decrease to 2-3
- You want exhaustive coverage: increase to 10+

## Interpreting Results

### Coverage Percentage

- **100%**: All isoforms covered (ideal)
- **80-99%**: Excellent coverage; uncovered isoforms may be rare
- **60-79%**: Good coverage; check if uncovered isoforms are important
- **<60%**: Poor coverage; consider:
  - Relaxing distance threshold
  - Checking if primers have off-target issues
  - Designing additional primers manually

### Primer Rank

Rank 1 = most important primer (covers most isoforms)
- For validation, test primers in rank order
- If resources limited, use top 2-3 primers
- Lower-ranked primers often cover rare isoforms

### Distance to End

- **<200nt**: Excellent for RT-qPCR
- **200-400nt**: Good for most applications
- **>400nt**: May have reduced efficiency; test empirically

## Troubleshooting

### "No specific primers found"

**Cause:** All primers have off-target perfect matches

**Solutions:**
1. Check primer design parameters (increase specificity in primer3)
2. Review alignment summary for common off-targets
3. Consider using a more stringent transcriptome reference

### "No primers within distance threshold"

**Cause:** All primers are too far from transcript end

**Solutions:**
1. Increase distance threshold (e.g., `-d 600`)
2. Check if peak calling/primer design is targeting wrong regions
3. Verify transcript annotations (some may have wrong 3' ends)

### "Coverage < 50%"

**Possible causes:**
- Gene has many divergent isoforms with different 3' UTRs
- Some isoforms lack the target peak region
- Distance threshold is too stringent

**Solutions:**
1. Increase max_primers to 7-10
2. Relax distance threshold
3. Design primers for multiple peak regions
4. Accept that some isoforms may not be amplifiable with common primers

## Integration with Pipeline

This tool is designed to work with the PeakPrime pipeline output:

```bash
# 1. Run PeakPrime pipeline with alignment
nextflow run main.nf --bam sample.bam --genes targets.txt \
  --transcriptome_index transcriptome --transcriptome_fasta transcripts.fa

# 2. Select optimal primers
python bin/select_optimal_primer_set.py \
  -i results/my_run/primer_alignment_summary.tsv \
  --all-genes -o results/my_run/optimal_primers.tsv

# 3. View results in Shiny app
Rscript app_standalone.R results/my_run
```

## Algorithm Complexity

- **Time:** O(n × m) where n = primers, m = isoforms
  - Typical: <1 second per gene
  - Large genes (>50 primers, >20 isoforms): 1-5 seconds

- **Space:** O(n × m) for primer-isoform mapping
  - Typical: <1 MB per gene

## References

**Set Cover Problem:**
- Greedy approximation: Chvatal, 1979
- ln(n) approximation bound: Johnson, 1974

**Primer Design for Isoforms:**
- RT-qPCR primer design: Bustin et al., 2009
- Isoform-specific amplification: Leparc et al., 2009

## Citation

If you use this tool, please cite the PeakPrime pipeline paper (in preparation).

---

**Author:** PeakPrime Development Team  
**Version:** 1.0  
**Last Updated:** October 2025
