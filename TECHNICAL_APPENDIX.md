# PeakPrime Technical Appendix

## Appendix A: File Format Specifications

### A.1 Input File Formats

#### Gene List Format
```
# gene_list.txt - One Ensembl gene ID per line
ENSG00000067191
ENSG00000089009  
ENSG00000108468
ENSG00000111249
ENSG00000115414
```

#### Primer3 Settings Format
```
# primer3_settings.txt - Key-value pairs for Primer3 configuration
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=0
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_MIN_SIZE=16
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=27
PRIMER_LEFT_MIN_SIZE=18
PRIMER_LEFT_OPT_SIZE=20
PRIMER_LEFT_MAX_SIZE=22
PRIMER_LEFT_MIN_TM=57.0
PRIMER_LEFT_OPT_TM=60.0
PRIMER_LEFT_MAX_TM=63.0
PRIMER_LEFT_MIN_GC=35.0
PRIMER_LEFT_MAX_GC=65.0
PRIMER_RIGHT_MIN_SIZE=16
PRIMER_RIGHT_OPT_SIZE=16
PRIMER_RIGHT_MAX_SIZE=18
PRIMER_RIGHT_MIN_TM=50.0
PRIMER_RIGHT_OPT_TM=54.0
PRIMER_RIGHT_MAX_TM=58.0
PRIMER_PRODUCT_SIZE_RANGE=40-200
PRIMER_NUM_RETURN=20
PRIMER_EXPLAIN_FLAG=1
=
```

### A.2 Output File Examples

#### Selected Peaks TSV Format
```
gene	chr	start	end	strand	peak_score	peak_pvalue	peak_qvalue	exonic_fraction	selection_method
ENSG00000067191	chr1	12345678	12345798	+	45.2	1.5e-05	2.3e-04	0.95	best_peak
ENSG00000089009	chr2	23456789	23456909	-	32.1	3.2e-04	1.1e-03	0.87	best_peak
ENSG00000108468	chr3	34567890	34568010	+	28.9	5.7e-04	1.8e-03	0.92	best_peak
```

#### cDNA Primers TSV Format
```
gene	primer_sequence	primer_type	start_pos	end_pos	tm	gc_percent	product_size
ENSG00000067191	ATCGATCGATCGATCGATCG	LEFT	1	20	60.2	50.0	120
ENSG00000067191	CGATCGATCGATCGATCGAT	RIGHT	101	120	59.8	50.0	120
ENSG00000089009	GCTAGCTAGCTAGCTAGCTA	LEFT	1	20	59.5	50.0	95
ENSG00000089009	TAGCTAGCTAGCTAGCTAGC	RIGHT	76	95	60.1	50.0	95
```

## Appendix B: Parameter Optimization Guidelines

### B.1 MACS2 Parameter Tuning

#### P-value Threshold Selection
- **High coverage samples**: Use stricter thresholds (0.01-0.001)
- **Low coverage samples**: Use relaxed thresholds (0.05-0.1)
- **Exploratory analysis**: Start with 0.05, adjust based on results

#### Peak Score Filtering
- **Default**: No minimum score (0)
- **High stringency**: Minimum score 10-20
- **Quality-focused**: Minimum score 5-10

### B.2 Quality Control Parameter Guidelines

#### Exonic Fraction Requirements
- **Exon-specific primers**: 0.8-0.9 minimum exonic fraction
- **General transcripts**: 0.5-0.7 minimum exonic fraction
- **Exploratory**: No minimum (null)

#### Coverage Smoothing
- **Default**: 31 bp window (odd number required)
- **High resolution**: 11-21 bp window
- **Low noise**: 41-61 bp window

### B.3 Primer3 Optimization

#### Product Size Ranges by Application
- **qPCR**: 75-150 bp (optimal for efficiency)
- **Standard PCR**: 100-300 bp
- **Sequencing primers**: 150-250 bp
- **Cloning**: 200-500 bp

#### Temperature Optimization
- **Standard conditions**: 57-63°C range
- **High stringency**: 60-65°C range
- **Multiplex PCR**: Narrow range (59-61°C)

## Appendix C: Computational Resource Guidelines

### C.1 Memory Requirements by Dataset Size

#### BAM File Size Guidelines
- **< 1 GB**: 4 GB RAM sufficient
- **1-5 GB**: 8 GB RAM recommended
- **5-20 GB**: 16 GB RAM recommended
- **> 20 GB**: 32+ GB RAM required

#### Gene List Size Impact
- **< 50 genes**: Minimal impact on resources
- **50-500 genes**: Linear scaling with memory
- **> 500 genes**: Consider batch processing

### C.2 Runtime Estimates

#### Local Execution (single core)
- **Peak calling**: 5-30 minutes depending on BAM size
- **Peak processing**: 2-10 minutes depending on gene count
- **Primer design**: 1-5 minutes per 100 genes
- **Transcriptome QC**: 10-60 minutes depending on primer count
- **Visualization**: 1-2 minutes per gene

#### HPC Execution (parallel)
- **Overall runtime**: 10-60 minutes for typical datasets
- **Scalability**: Linear with number of available cores
- **Bottlenecks**: Peak calling and transcriptome alignment

### C.3 Storage Requirements

#### Temporary Files
- **MACS2 outputs**: 100MB - 1GB depending on peak count
- **Intermediate sequences**: 10-100MB depending on target count
- **BigWig files**: 10-50% of original BAM size

#### Final Outputs
- **Essential results**: 10-100MB for most analyses
- **With visualization**: +50-200MB for plot files
- **With QC alignments**: +100MB-1GB for alignment files

## Appendix D: Integration with External Tools

### D.1 Genome Browser Integration

#### BigWig File Usage
```bash
# Upload to UCSC Genome Browser
# Use generated *.bam.bw files as custom tracks

# IGV (Integrative Genomics Viewer)
igv.sh -g hg38 sample.bam.bw

# Command line viewing with pyBigWig
import pyBigWig
bw = pyBigWig.open("sample.bam.bw")
values = bw.values("chr1", 1000000, 1001000)
```

#### BED File Integration
```bash
# Load primer targets in genome browsers
# Use selected_peaks.bed as custom annotation track

# Convert to other formats
bedToBigBed selected_peaks.bed hg38.chrom.sizes selected_peaks.bb
```

### D.2 Downstream Analysis Integration

#### Primer Database Integration
```python
# Example: Import primers into database
import pandas as pd
primers = pd.read_csv("cdna_primers.tsv", sep="\t")
# Insert into primer management system
```

#### Experimental Validation Pipeline
```bash
# Generate primer order sheets
awk -F'\t' 'NR>1{print $1"\t"$2"\t"$6"°C\t"$7"%GC"}' cdna_primers.tsv > primer_orders.txt

# Quality prediction scores
# Integrate with primer analysis tools like Primer-BLAST
```

### D.3 Workflow Management Integration

#### Snakemake Integration
```python
# Snakemake rule example
rule peakprime_design:
    input:
        bam="samples/{sample}.bam",
        gtf="reference/annotations.gtf",
        genes="gene_lists/{geneset}.txt"
    output:
        primers="results/{sample}_{geneset}/cdna_primers.tsv"
    shell:
        "nextflow run PeakPrime/main.nf --bam {input.bam} --gtf {input.gtf} "
        "--genes {input.genes} --outdir results/{wildcards.sample}_{wildcards.geneset}"
```

#### CWL Integration
```yaml
# Common Workflow Language specification
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [nextflow, run, main.nf]
inputs:
  bam_file:
    type: File
    inputBinding:
      prefix: --bam
  gtf_file:
    type: File
    inputBinding:
      prefix: --gtf
  gene_list:
    type: File
    inputBinding:
      prefix: --genes
```

## Appendix E: Quality Control Metrics Interpretation

### E.1 Peak Quality Assessment

#### Peak Score Interpretation
- **Score > 50**: Excellent peak quality, high confidence
- **Score 20-50**: Good peak quality, acceptable for most applications
- **Score 10-20**: Moderate peak quality, may require validation
- **Score < 10**: Low peak quality, consider alternative targets

#### P-value Guidelines
- **P < 0.001**: Very high confidence peaks
- **P 0.001-0.01**: High confidence peaks
- **P 0.01-0.05**: Moderate confidence peaks
- **P > 0.05**: Low confidence, may be false positives

### E.2 Primer Quality Metrics

#### Melting Temperature Assessment
- **ΔTm < 2°C**: Excellent primer pair balance
- **ΔTm 2-5°C**: Acceptable for most applications
- **ΔTm > 5°C**: May cause PCR efficiency issues

#### GC Content Guidelines
- **40-60%**: Optimal range for most applications
- **35-40% or 60-65%**: Acceptable with careful optimization
- **< 35% or > 65%**: May require special PCR conditions

### E.3 Transcriptome QC Classification

#### Quality Categories
- **PERFECT**: Single perfect match, no off-targets
- **GOOD**: Strong target match, minimal cross-reactivity
- **MODERATE**: Acceptable match with some off-targets
- **POOR**: Multiple strong matches or weak target match
- **FAIL**: No acceptable target matches

#### Recommended Actions
- **PERFECT/GOOD**: Use directly for experiments
- **MODERATE**: Consider additional validation
- **POOR/FAIL**: Redesign with different parameters

## Appendix F: Troubleshooting Decision Tree

### F.1 No Peaks Detected

```
No MACS2 peaks found
├── Check BAM file
│   ├── Indexed properly? → Fix indexing
│   ├── Contains target regions? → Verify alignment
│   └── Sufficient coverage? → Check sequencing depth
├── Adjust MACS2 parameters
│   ├── Increase p-value threshold → Try 0.1 or 0.2
│   ├── Remove peak score filter → Set to 0
│   └── Check MACS2 logs → Review error messages
└── Verify gene list
    ├── Correct gene IDs? → Match with GTF
    ├── Genes in BAM regions? → Check coordinate overlap
    └── Alternative gene names? → Use Ensembl IDs
```

### F.2 Poor Primer Quality

```
Primer3 fails or poor primers
├── Check target sequences
│   ├── Too short? → Increase padding
│   ├── Repetitive? → Filter low complexity
│   └── GC extreme? → Adjust GC parameters
├── Adjust Primer3 settings
│   ├── Relax temperature range → ±3°C from optimal
│   ├── Increase size range → Allow broader lengths
│   └── Reduce stringency → Lower specificity requirements
└── Review peak selection
    ├── Poor peak quality? → Increase score threshold
    ├── Low exonic overlap? → Adjust exonic fraction
    └── Edge effects? → Increase padding
```

### F.3 Performance Issues

```
Slow execution or failures
├── Memory issues
│   ├── Increase allocation → Use --max_memory
│   ├── Use HPC profile → Switch to PBS
│   └── Process in batches → Split gene lists
├── CPU bottlenecks
│   ├── Parallel execution → Use multiple cores
│   ├── Optimize I/O → Use fast storage
│   └── Check dependencies → Update software
└── Disk space
    ├── Clean temp files → Use --cleanup
    ├── Optimize outputs → Remove unnecessary files
    └── Monitor usage → Use df and du commands
```

---

*Technical Appendix for PeakPrime version 1.0*  
*Last updated: [Current Date]*