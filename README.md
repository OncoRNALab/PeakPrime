# PeakPrime — RNA-seq primer design pipeline


PeakPrime quickly finds peak-covered, exon-aware windows in RNA-seq data and designs strand-correct cDNA primers with robust QC.

Why PeakPrime?

- Peak-centered first: fast and biologically sensible default
- Sliding-window rescue: can search for better windows when the peak window fails QC
- Strand-aware primer selection for cDNA assays
- Optional transcriptome alignment QC to detect cross-reactivity

## Features

- **Isoform‑agnostic Coverage analysis**: Extracts high-coverage regions from RNA-seq BAM files using BigWig conversion
- **Advanced QC**: Coverage trimming, gap detection, sliding-window optimization
- **Strand-specific primer selection**: Ensures primers match mRNA orientation for cDNA amplification
- **Peak detection**: Identifies optimal target regions with configurable parameters
- **Primer design**: Uses Primer3 for robust primer pair generation
- **Transcriptome alignment QC**: Bowtie2-based specificity checks and Python analysis

## Requirements

- Nextflow (≥22.04.0)
- Conda/Mamba for environment management
- Required input files:
  - BAM file from RNA-seq alignment
  - GTF annotation file
  - Gene list (one gene per line)
  - Primer3 settings file
  - Optional: Pre-built Bowtie2 transcriptome index

## Quick Start

## Example Usage

### Basic Usage

```bash
nextflow run main.nf \
  --bam sample.bam \
  --gtf annotations.gtf \
  --genes gene_list.txt \
  --primer3_settings primer3_settings.txt \
  --outdir results/
```

### With Advanced Quality Control

```bash
nextflow run main.nf \
  --bam sample.bam \
  --gtf annotations.gtf \
  --genes gene_list.txt \
  --primer3_settings primer3_settings.txt \
  --sliding_window \
  --min_window_mean_pct 20 \
  --max_gap 50 \
  --trim_low_coverage_pct 10 \
  --trim_to_exon \
  --min_exonic_fraction 0.8 \
  --outdir results/
```

### With Transcriptome QC and Gene Name Mapping

```bash
# Using pre-built Bowtie2 index with gene name mapping (recommended)
nextflow run main.nf \
  --bam sample.bam \
  --gtf annotations.gtf \
  --genes gene_list.txt \
  --primer3_settings primer3_settings.txt \
  --transcriptome_index /path/to/transcriptome_index_prefix \
  --transcriptome_fasta /path/to/transcriptome.fasta \
  --max_primers_per_gene 5 \
  --outdir results/
```

### Complete Example with All Features

```bash
nextflow run main.nf \
  --bam sample_rnaseq.bam \
  --gtf gencode.v44.annotation.gtf \
  --genes target_genes.txt \
  --primer3_settings config/primer3_settings.txt \
  --fasta GRCh38.primary_assembly.fa \
  --pad 75 \
  --smooth_k 51 \
  --sliding_window \
  --min_window_mean_pct 15 \
  --max_gap 30 \
  --trim_low_coverage_pct 5 \
  --trim_to_exon \
  --min_exonic_fraction 0.75 \
  --transcriptome_index indexes/gencode_transcriptome \
  --transcriptome_fasta gencode.v44.transcripts.fa \
  --max_primers_per_gene 3 \
  --outdir results_comprehensive/
```

### Setting Up Transcriptome Index

To use transcriptome alignment QC, first build a Bowtie2 index from your transcriptome FASTA:

```bash
# Build Bowtie2 index (one-time setup)
bowtie2-build transcriptome.fasta transcriptome_index

# Then use the prefix in your pipeline
nextflow run main.nf \
  --transcriptome_index transcriptome_index \
  [other parameters...]
```

## Parameters

### Required Parameters

- `--bam`: Path to RNA-seq BAM file
- `--gtf`: Path to GTF annotation file  
- `--genes`: Path to gene list file (one gene per line)
- `--primer3_settings`: Path to Primer3 settings file
- `--outdir`: Output directory for results

### Optional Parameters

#### Coverage Analysis
- `--fasta`: Reference genome FASTA file for sequence extraction
- `--pad`: Padding around gene regions in bp (default: 60)
- `--smooth_k`: Smoothing kernel size for coverage (default: 31, must be odd)
- `--genome_package`: R/Bioconductor genome package name (default: 'BSgenome.Hsapiens.UCSC.hg38')

#### Peak Detection and Window Selection
- `--sliding_window`: Enable sliding window analysis for better peak selection
- `--min_window_mean`: Minimum mean coverage required across the final window (absolute value)
- `--min_window_mean_pct`: Minimum coverage as percentage of gene's peak coverage (0-100, overrides `--min_window_mean`)
- `--max_gap`: Maximum allowed longest zero-coverage run within the final window
- `--search_slop`: Extra bases around exon span for coverage import (default: 1000)

#### Window Refinement
- `--trim_to_exon`: Trim final window to boundaries of exon containing the peak (prevents intronic spillover)
- `--trim_low_coverage_pct`: Trim window ends with coverage below X% of window peak (0-100, e.g., 10 for 10%)
- `--min_exonic_fraction`: Minimum required exonic fraction of window (0-1, for quality control)

#### Transcriptome Alignment QC
- `--transcriptome_index`: Path to pre-built Bowtie2 transcriptome index prefix
- `--transcriptome_fasta`: Transcriptome FASTA file for gene name mapping (optional but recommended)
- `--max_primers_per_gene`: Maximum primers per gene for alignment (default: 3)

## Output Files

### Core Pipeline Outputs
- `primer_targets.fa`: Target sequences for primer design in FASTA format
- `primer_targets.bed`: Genomic coordinates of target regions in BED format
- `peaks.tsv`: Detailed peak information with coverage statistics and quality metrics
- `qc_coverage_summary.tsv`: Comprehensive quality control summary for each gene
- `primer3_input.txt`: Formatted input file for Primer3
- `primer3_output.txt`: Raw Primer3 results with all primer candidates
- `cdna_primers.tsv`: Final cDNA-appropriate primer pairs (strand-specific selection)

### Transcriptome Alignment QC (if enabled)
- `primers_for_alignment.fa`: Primer sequences in FASTA format for alignment
- `primers_alignment.bam`: Bowtie2 alignment results of primers against transcriptome
- `primers_alignment.bam.bai`: BAM index file
- `alignment_stats.txt`: Bowtie2 alignment statistics
- `primer_alignment_report.tsv`: Summary statistics and quality classification per primer
- `primer_alignment_summary.tsv`: Detailed per-alignment summary with gene names (MAPQ=255 only)

### Output File Details

#### `peaks.tsv`
Contains information about selected primer target regions:
- `gene`: Gene identifier
- `chr`: Chromosome
- `start`/`end`: Genomic coordinates of target window
- `strand`: Gene strand orientation
- `exonic_fraction`: Fraction of window that overlaps exons (0-1)
- `trimmed_to_exon`: Whether window was trimmed to exon boundaries
- `fail_exonic_fraction`: Whether window failed exonic fraction threshold

#### `qc_coverage_summary.tsv`  
Comprehensive quality control metrics for each gene:
- `gene`: Gene identifier
- `total_exonic_bases`: Total exonic length for the gene
- `max_cov`/`mean_cov`/`median_cov`: Coverage statistics across all exons
- `peak_pos`: Position of maximum coverage after smoothing
- `window_start`/`window_end`: Final target window coordinates
- `window_mean`/`window_median`/`window_min`/`window_max`: Window coverage statistics
- `window_zeros`: Number of zero-coverage positions in window
- `longest_zero_run`: Length of longest consecutive zero-coverage stretch
- `pass_min_mean`/`pass_max_gap`: Quality control pass/fail flags
- `strategy`: Window selection strategy used (e.g., "peak_centered", "sliding_best", "peak_centered+trimmed")
- `exonic_fraction`: Fraction of final window overlapping exons
- `fail_exonic_fraction`: Whether window failed exonic fraction requirements
- `trimmed_to_exon`: Whether window was trimmed to exon boundaries

#### `cdna_primers.tsv`
Final strand-appropriate primer pairs for cDNA amplification:
- `gene_id`: Target gene identifier
- `primer_index`: Primer pair index (from Primer3)
- `primer_type`: LEFT or RIGHT primer
- `primer_sequence`: Primer sequence
- `primer_tm`: Melting temperature
- `primer_gc`: GC content percentage
- `gene_strand`: Strand orientation of target gene

#### `primer_alignment_report.tsv` (if transcriptome QC enabled)
Summary statistics per primer:
- `gene_id`: Target gene identifier
- `primer_index`: Primer pair index
- `primer_type`: LEFT or RIGHT primer
- `primer_sequence`: Primer sequence
- `gene_strand`: Target gene strand
- `num_alignments`: Total number of transcriptome alignments
- `alignment_quality`: Quality classification (PERFECT/GOOD/MODERATE/POOR/FAIL)
- `best_mapq`: Highest MAPQ score among all alignments

#### `primer_alignment_summary.tsv` (if transcriptome QC enabled)
Detailed alignment information for primers that align in forward direction only (5' to 3' orientation):
- `gene_id`: Target gene identifier
- `primer_index`: Primer pair index
- `primer_type`: LEFT or RIGHT primer
- `primer_sequence`: Primer sequence (reported in 5' to 3' orientation)
- `gene_strand`: Target gene strand
- `aligned_transcript`: Transcript ID from alignment
- `aligned_gene_name`: Gene name extracted from transcriptome FASTA headers
- `alignment_start`/`alignment_end`: Alignment coordinates on transcript
- `alignment_length`: Length of alignment
- `alignment_strand`: Always "forward" (reverse alignments are filtered out)
- `alignment_mapq`: Mapping quality score (255 = perfect/unique alignment)
- `mismatches`: Number of mismatches in alignment
- `distance_to_end`: Distance from primer start to transcript end (strand-specific)
- `transcript_length`: Total length of aligned transcript

## Quality Control Features

### Coverage Analysis and Peak Selection
The pipeline implements sophisticated coverage analysis to identify optimal primer target regions:

- **Smoothed Coverage**: Uses running mean smoothing to reduce noise in coverage data
- **Peak Detection**: Identifies the highest coverage position within each gene's exonic regions
- **Window Selection**: Creates windows around peaks with configurable padding

### Advanced Window Quality Control

#### Sliding Window Analysis
When `--sliding_window` is enabled, the pipeline can search for better target windows if the initial peak-centered window fails quality thresholds:
- Searches across the entire exonic span of the gene
- Finds the window with maximum mean coverage that meets quality criteria
- Respects gap and coverage thresholds during search

#### Coverage-Based Trimming
Multiple trimming strategies ensure high-quality target regions:

- **Low Coverage Trimming** (`--trim_low_coverage_pct`): Removes window ends with coverage below a specified percentage of the window's peak coverage
- **Exonic Boundary Trimming** (`--trim_to_exon`): Trims windows to the boundaries of the exon containing the coverage peak
- **Quality Thresholds**: Configurable minimum mean coverage requirements (absolute or percentage-based)

#### Gap Analysis
- **Zero-Coverage Detection**: Identifies stretches of zero coverage within target windows  
- **Gap Filtering** (`--max_gap`): Rejects or searches for alternative windows if gaps exceed thresholds
- **Quality Classification**: Flags windows that fail quality criteria

### Strand-Specific Primer Selection
For cDNA amplification, the pipeline selects biologically appropriate primers based on correct orientation:

- **Positive-strand genes**: Uses RIGHT primers (reverse complement of template matches mRNA)
- **Negative-strand genes**: Uses LEFT primers (direct match of template matches mRNA) 
- **Biological Logic**: Ensures primers match the mRNA sequence orientation for proper cDNA amplification
- **Template Handling**: Genomic sequences are automatically extracted in proper strand orientation

### Transcriptome Alignment QC
Fast Python-based quality control through transcriptome alignment:

#### Alignment Strategy
- **Comprehensive Mapping**: Uses Bowtie2 with `-a` flag to report all valid alignments
- **Forward Strand Filtering**: Only reports primers that align in forward direction (5' to 3')
- **Specificity Assessment**: Evaluates primer specificity against the entire transcriptome
- **Cross-reactivity Detection**: Identifies primers that may amplify unintended targets

#### Quality Classification System
Primers are classified based on transcriptome alignment patterns:
- **PERFECT (1 alignment)**: Unique alignment, highest specificity
- **GOOD (2-5 alignments)**: Limited cross-reactivity, generally acceptable
- **MODERATE (6-20 alignments)**: Some cross-reactivity, use with caution
- **POOR (>20 alignments)**: High cross-reactivity, likely problematic
- **FAIL (0 alignments)**: No transcriptome match, may indicate design issues

#### Gene Name Resolution
- **FASTA Header Parsing**: Extracts gene names from transcriptome FASTA headers
- **Cross-reference Mapping**: Maps transcript IDs to actual gene symbols
- **Specificity Context**: Shows which genes primers actually target vs. intended targets

## Example Files

### Gene List Format
```
GENE1
GENE2
GENE3
```

### Primer3 Settings Example
```
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_NUM_RETURN=5
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=20.0
PRIMER_MAX_GC=80.0
PRIMER_MAX_POLY_X=4
PRIMER_INTERNAL_MAX_POLY_X=4
PRIMER_SALT_MONOVALENT=50.0
PRIMER_DNA_CONC=50.0
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_SELF_ANY=12
PRIMER_MAX_SELF_END=8
PRIMER_PAIR_MAX_COMPL_ANY=12
PRIMER_PAIR_MAX_COMPL_END=8
=
```

## Troubleshooting

### Common Issues

1. **Missing BAM index**: The pipeline automatically creates BAM indices if missing
2. **Memory issues**: Adjust Nextflow executor settings for large BAM files
3. **Gene not found**: Ensure gene names in the gene list match those in the GTF file
4. **No primers generated**: Check coverage depth and adjust quality control parameters

### Performance Tips

1. **Use pre-built transcriptome indices**: Much faster than building from FASTA each time
2. **Limit primer count**: Use `--max_primers_per_gene` to reduce alignment time
3. **Adjust coverage parameters**: Fine-tune quality thresholds based on your data depth

## Citation

If you use this pipeline in your research, please cite the relevant tools:
- Nextflow
- Primer3
- Bowtie2 (if using transcriptome alignment QC)
- R/Bioconductor packages used in the analysis
