# Exonic Filter Failure Analysis: ENSG00000108960

## Executive Summary

**Gene:** ENSG00000108960  
**Chromosome:** chr17 (minus strand)  
**Gene Coordinates:** 55392622-55421924 (29,303 bp)  
**Peak Coordinates:** 55419614-55419911 (298 bp)  
**Exonic Overlap:** 0 bp  
**Exonic Fraction:** 0.00%  
**Failure Reason:** Failed forced exonic trimming (no exonic overlap or trimmed region too short)  

**Conclusion:** ❌ The peak falls entirely in an **intronic region** and has zero overlap with any exon.

---

## Detailed Analysis

### 1. Peak Location

```
Peak: chr17:55419614-55419911
Length: 298 bp
```

### 2. Exon Coordinates

Retrieved 10 unique exons from GTF file for ENSG00000108960:

| Exon # | Coordinates | Distance from Peak |
|--------|-------------|-------------------|
| 1 | 55392622-55394534 | 25,080 bp upstream |
| 2 | 55392643-55394534 | 25,080 bp upstream |
| 3 | 55401469-55401538 | 18,076 bp upstream |
| 4 | 55403767-55403868 | 15,746 bp upstream |
| 5 | 55404455-55404547 | 15,067 bp upstream |
| 6 | 55407746-55407820 | 11,794 bp upstream |
| 7 | 55411257-55411417 | 8,197 bp upstream |
| 8 | 55414151-55414232 | 5,382 bp upstream |
| **9** | **55420092-55420533** | **181 bp downstream** ⭐ |
| 10 | 55421670-55421924 | 1,759 bp downstream |

⭐ **Closest exon**: Exon 9 (55420092-55420533) starts 181 bp **after** the peak ends

### 3. Overlap Calculation

```python
Peak:     chr17:55419614--------55419911  (298 bp)
                                    ↓ 181 bp gap
Exon 9:   chr17:              55420092--------55420533  (442 bp)
```

**Result:**
- Total exonic overlap: **0 bp**
- Exonic fraction: **0.0000** (0.00%)
- Peak falls entirely in **intron** between exon 8 and exon 9

---

## Why This Fails the Filter

### Exonic Filter Requirements

The pipeline applies a **forced exonic trimming** filter that requires:
1. Peak must overlap at least one exon
2. Trimmed exonic region must be long enough for primer design

### Why ENSG00000108960 Fails

1. **No exonic overlap**: Peak is entirely intronic (gap of 181 bp to nearest exon)
2. **Cannot trim to exons**: No exonic sequence available in peak region
3. **Primer design impossible**: Primers require exonic targets for proper amplification

---

## Biological Interpretation

### What Does an Intronic Peak Mean?

An intronic peak can indicate several biological scenarios:

#### 1. **Immature/Unspliced RNA** ✅ Most likely
- Pre-mRNA contains introns before splicing
- Peak may represent nascent transcription
- Common in nuclear RNA or RNA-seq with incomplete splicing

#### 2. **Retained Intron**
- Some isoforms may retain this intron
- Check for alternative splicing events
- May be tissue-specific or condition-specific

#### 3. **Non-coding RNA**
- Peak may represent a long non-coding RNA (lncRNA)
- Could be an intronic regulatory element
- Check GENCODE/Ensembl for annotated ncRNAs

#### 4. **Technical Artifacts**
- Genomic DNA contamination
- Mapping artifacts
- Peak calling false positive

#### 5. **Pre-mRNA Stabilization**
- Intron may be retained in stable pre-mRNA
- Possible regulatory mechanism
- Check RNA stability data

### Gene-Specific Context

**ENSG00000108960** - This gene should be looked up to understand:
- What protein does it encode?
- Is it known to have retained introns?
- Are there annotated regulatory elements in this region?
- What is the expression pattern in your tissue/condition?

---

## Recommendations

### Immediate Actions

1. **Visualize the peak** in IGV or UCSC Genome Browser:
   ```bash
   # Open IGV and load:
   # - Your BigWig file
   # - MACS2 peak file
   # - GTF annotation
   # Navigate to: chr17:55419614-55419911
   ```

2. **Check if other genes have similar issues**:
   ```bash
   # Count how many genes failed exonic filter
   grep "Failed forced exonic trimming" results/test/peaks_qc_summary.tsv | wc -l
   
   # List all failed genes
   grep "Failed forced exonic trimming" results/test/peaks_qc_summary.tsv | cut -f1
   ```

3. **Review peak calling parameters**:
   - Are peaks too sensitive (low q-value threshold)?
   - Consider adjusting MACS2 parameters for stricter peak calling
   - Check if `--broad` mode is appropriate for your data

### Long-term Considerations

1. **Accept intronic peaks** if biologically relevant:
   - Some genes naturally have intronic expression
   - May need to adjust primer design strategy
   - Consider designing primers in flanking exons

2. **Filter out intronic peaks** if not relevant:
   - Keep current exonic filter
   - Focus on genes with clear exonic peaks
   - Reduces false positives in primer design

3. **Investigate retained introns**:
   - Use tools like LeafCutter or rMATS
   - Check for condition-specific splicing
   - May reveal interesting biology

---

## Command History

### Commands Used for Analysis

```bash
# 1. Extract exon coordinates from GTF
grep -w "ENSG00000108960" /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf | \
  grep -w "exon" | \
  awk '{print $1 "\t" $4 "\t" $5}' | \
  sort -k2,2n -k3,3n | \
  uniq

# 2. Calculate overlap with Python
python3 << 'EOF'
peak_start = 55419614
peak_end = 55419911
exons = [(55420092, 55420533)]  # Closest exon
overlap_start = max(peak_start, exons[0][0])
overlap_end = min(peak_end, exons[0][1])
overlap = max(0, overlap_end - overlap_start + 1)
print(f"Overlap: {overlap} bp")
EOF
```

---

## Related Genes

Check these other genes that failed the same filter:

```bash
# From peaks_qc_summary.tsv
ENSG00000177700  # Failed forced exonic trimming
ENSG00000112787  # Failed forced exonic trimming  
ENSG00000071894  # Failed forced exonic trimming
```

Investigate if they also have intronic peaks or different failure modes.

---

## Conclusion

**ENSG00000108960 failed the exonic filter correctly** - the peak has zero exonic overlap because it falls entirely in an intron. The filter is working as designed to exclude peaks that cannot be used for primer design on mature mRNA.

**Next steps:**
1. Decide if intronic peaks are biologically relevant for your study
2. Review other failed genes to understand the scope of intronic peaks
3. Consider adjusting peak calling or filtering parameters if needed
4. Visualize peaks in genome browser to confirm biological validity

---

*Analysis performed: 2025-01-XX*  
*Pipeline version: Nextflow DSL2*  
*GTF version: Ensembl 109*
