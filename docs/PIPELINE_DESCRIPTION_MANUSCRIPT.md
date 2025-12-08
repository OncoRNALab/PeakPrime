# Pipeline Description

## 1. Overview of the Pipeline

PeakPrime is a Nextflow-based computational pipeline that automates the design of strand-specific cDNA primers for targeted enrichment in 3′-end RNA-seq experiments. The pipeline integrates statistical peak detection on RNA-seq coverage profiles with automated primer design, quality control, and specificity screening. PeakPrime operates in a modular fashion: it identifies high-confidence transcript 3′ regions using MACS2 peak calling, constrains candidate windows to exonic sequences, designs primer pairs with Primer3, and validates specificity via transcriptome alignment using Bowtie2. The workflow is portable, reproducible, and scales to thousands of genes, producing publication-ready visualizations and comprehensive quality-control reports.

The pipeline implements two distinct operational modes. In **peak-based mode** (default), PeakPrime leverages empirical RNA-seq coverage data to identify statistically significant peaks in transcript 3′ regions, ensuring primers target genuine transcriptional hotspots. In **distance-based mode**, primers are designed at a fixed distance from annotated 3′ ends without requiring coverage data, enabling primer design from sequence information alone. Both modes support optional isoform-aware optimization to maximize transcriptome coverage with minimal primer sets.

## 2. Operating Modes

PeakPrime supports three primary operating configurations:

### Peak-Based Single-Peak Mode (Default)
Identifies the single highest-confidence peak per gene from RNA-seq coverage profiles and designs primers targeting that region. This mode is optimal for standard targeted enrichment where one primer set per gene is sufficient. The pipeline selects peaks using configurable metrics (peak score or q-value) and ranks peaks per gene to allow selection of the best (rank=1), second-best (rank=2), or any ranked peak.

### Peak-Based Multi-Peak Mode
Extracts multiple significant peaks per gene (configured via `--select_all_peaks` and `--max_peaks_per_gene`) and designs distinct primer sets for each. This mode is designed for genes with complex 3′ UTR structures, alternative polyadenylation sites, or scenarios requiring broader isoform coverage. When combined with `--optimize_multipeak`, the pipeline applies a distance-weighted scoring algorithm to select optimal primer combinations across peaks.

### Distance-Based Mode
Designs primers at a user-specified distance from annotated 3′ transcript ends without RNA-seq input. This mode is useful for preliminary primer design, cross-platform primer sets, or when empirical coverage data are unavailable. Primers are designed on a fixed-length template extracted from transcript sequences, either provided as FASTA or fetched from Ensembl.

## 3. Inputs and Configuration

Table A presents the complete input specification.

**Table A: Inputs**

| Name | Type/Format | Required? | Validation | Used in Steps | Notes |
|------|-------------|-----------|------------|---------------|-------|
| `bam` | BAM (indexed) | Required (peak-based modes) | Indexed (`.bam.bai` must exist), coordinate-sorted, splice-aware alignment recommended | Peak calling (MACS2), coverage generation (megadepth) | Must contain reads aligned to same genome build as GTF |
| `gtf` | GTF | Required | Valid GTF format (GFF3 not supported), gene_id and transcript_id attributes present | Peak processing, exon overlap calculation, plotting | Must match genome build of BAM |
| `genes` | Text file | Required | One Ensembl gene ID per line (format: `ENSG[0-9]+`), no header | Peak selection, filtering, plotting | Gene IDs must exist in GTF |
| `genome_package` | R/Bioconductor BSgenome package name | Optional | Valid BSgenome package name installed in R environment | Sequence extraction for primer design | Default: `BSgenome.Hsapiens.UCSC.hg38` |
| `fasta` | FASTA | Optional | Valid FASTA format, indexed (`.fai`) | Alternative to BSgenome for sequence extraction | Use when BSgenome package unavailable |
| `transcriptome_index` | Bowtie2 index prefix | Optional | All index files (`*.bt2`) must exist | Primer specificity checking via transcriptome alignment | Prefix only (e.g., `transcriptome_idx`) |
| `transcript_mapping` | TSV | Optional | Tab-delimited: transcript_id, gene_id, gene_name | Annotation of alignment results | Default: `resources/transcript_gene_mapping.GRCh38.109.tsv` |
| `transcript_fasta` | FASTA | Optional (distance mode) | Valid transcript FASTA | Distance-based mode sequence source | Alternative to Ensembl fetch |
| `primer3_settings` | Text | Required | Primer3-compatible parameter file | Primer3 configuration | Default: `config/primer3_settings.txt` |

### Configuration Parameters

All parameters are specified via command-line flags (`--parameter value`) or in `params.config`. Key parameters include:

- **MACS2 peak calling**: `--macs2_qvalue_threshold` (default: 0.05), `--macs2_min_peak_score` (default: 0), `--macs2_extsize` (default: null, auto-detect), `--macs2_shift` (default: null, auto-detect)
- **Peak selection**: `--peak_selection_metric` (default: 'score'), `--peak_rank` (default: 1)
- **Multi-peak mode**: `--select_all_peaks` (default: false), `--max_peaks_per_gene` (default: 0, unlimited)
- **Exonic trimming**: `--force_exonic_trimming` (default: true), `--min_trimmed_length` (default: 30 bp), `--min_exonic_fraction` (default: null)
- **Distance-based mode**: `--distance_mode` (default: false), `--template_length` (required for distance mode)
- **Optimization**: `--optimize_multipeak` (default: false), `--distance_weight` (default: 0.5), `--isoform_weight` (default: 0.3), `--peak_rank_weight` (default: 0.2)
- **Transcriptome QC**: `--max_primers_per_gene` (default: 20), `--distance_threshold` (default: 400 bp)

## 4. Preprocessing

Before peak calling, the pipeline performs minimal preprocessing:

1. **BAM Validation and Indexing**: Input BAM files are checked for coordinate-sorted order and indexed using `samtools index` (v1.16) if `.bam.bai` is absent. This step ensures compatibility with downstream tools requiring random access (MACS2, megadepth).

2. **GTF Filtering** (plotting mode only): When `--makeplots` is enabled, the GTF annotation is pre-filtered to retain only features (genes, transcripts, exons, UTRs) matching the input gene list. This optimization reduces memory consumption and accelerates feature extraction during visualization. Filtering is performed using R/Bioconductor `rtracklayer` (version not specified).

3. **Gene List Validation**: The input gene list is parsed to verify that all identifiers match the pattern `ENSG[0-9]+` and exist in the provided GTF. Genes absent from the annotation trigger warnings but do not halt execution.

Assumption: Input BAM files represent spliced alignments generated by splice-aware aligners (e.g., STAR, HISAT2). The pipeline does not perform read alignment or quality filtering; these steps must be completed upstream.

## 5. Step-by-Step Processing Stages

### Stage 1: Peak Calling with MACS2

**Purpose**: Identify statistically significant coverage peaks in 3′ regions of target genes from RNA-seq BAM files.

**Algorithm**: MACS2 (Model-based Analysis of ChIP-Seq) is applied in single-sample mode to call peaks. MACS2 models fragment size distribution (unless `--nomodel` is specified via `--macs2_extsize`), calculates local enrichment statistics, and identifies regions exceeding a q-value threshold.

**Parameters**:
- `--gsize hs` (human genome size)
- `--qvalue` (default: 0.05, configurable via `--macs2_qvalue_threshold`)
- `--bdg` (generate bedGraph for coverage tracks)
- `--keep-dup auto` (automatic duplicate handling)
- `--nomodel --extsize <value> --shift <value>` (only when `--macs2_extsize` is specified)

**Key Data Structures**: Produces narrowPeak format files (BED6+4) with columns: chr, start, end, peak_name, score, strand, fold_enrichment, -log10(p-value), -log10(q-value), relative_summit_position.

**Complexity**: O(N) where N = number of reads. Runtime scales linearly with BAM file size; typical execution time is 5–20 minutes for 10–50 million reads on a single core.

**Failure Modes**: 
- Insufficient coverage depth (< 1M reads) may prevent fragment size modeling; use `--macs2_extsize 50 --macs2_shift -25` for 3′ RNA-seq data.
- Genes with no reads produce zero peaks; logged but non-fatal.

**Module**: `MACS2_CALLPEAK`

### Stage 2: Coverage Track Generation

**Purpose**: Generate normalized BigWig coverage files for visualization and QC.

**Algorithm**: Uses `megadepth` (version not specified) to compute per-base coverage from BAM files and output in BigWig format. Coverage is reported as raw counts (unnormalized).

**Parameters**: None (uses defaults).

**Complexity**: O(N), where N = genome positions with coverage. Typically completes in 2–5 minutes for human RNA-seq.

**Module**: `MEGADEPTH_BW`

### Stage 3: Peak Processing and Selection

**Purpose**: Filter MACS2 peaks by quality metrics, select best peak(s) per gene, constrain to exonic regions, and extract target sequences.

**Algorithm**: R script (`bin/process_macs2_peaks.R`) performs the following operations:
1. Import narrowPeak file and GTF annotation
2. Compute overlap between peaks and exonic features (exons + UTRs) per gene
3. Calculate exonic fraction: (peak ∩ exons) / peak_length
4. Rank peaks per gene by selected metric (score or q-value)
5. Select peak(s) per gene: rank 1 (best), rank 2 (second-best), or all peaks (multi-peak mode)
6. Apply quality filters: `q-value ≤ macs2_qvalue_threshold`, `score ≥ macs2_min_peak_score`, `exonic_fraction ≥ min_exonic_fraction` (if specified)
7. Trim peak boundaries to exonic regions if `--force_exonic_trimming` is enabled
8. Extract genomic sequences for surviving peaks using BSgenome or FASTA

**Key Parameters**:
- `--peak_selection_metric`: 'score' (default) or 'qvalue'
- `--peak_rank`: integer ≥ 1 (default: 1)
- `--min_exonic_fraction`: 0–1 (default: null, no filter)
- `--force_exonic_trimming`: boolean (default: true)
- `--min_trimmed_length`: minimum bp after trimming (default: 30)

**Data Structures**:
- Input: GRanges objects (peaks, exons)
- Output: data.frame with columns: gene_id, chr, start, end, strand, peak_rank, peak_score, qvalue, pvalue, exonic_fraction, trimmed_start, trimmed_end

**Complexity**: O(P × E), where P = number of peaks, E = number of exonic features. Typically < 1 minute for 1000 genes.

**Failure Modes**:
- Peaks entirely outside exons are discarded if `--min_exonic_fraction` > 0
- Trimmed regions shorter than `--min_trimmed_length` are rejected (logged in QC summary)

**Outputs**:
- `selected_peaks.tsv`: Peak coordinates and QC metrics
- `primer_targets.bed`: BED6 format target regions
- `primer_targets.fa`: FASTA sequences for Primer3 input
- `peaks_qc_summary.tsv`: Comprehensive QC table

**Module**: `PROCESS_MACS2_PEAKS`

### Stage 4: Primer3 Input Formatting

**Purpose**: Convert FASTA sequences to Primer3 input format with sequence IDs and target region annotations.

**Algorithm**: Python script (`bin/fasta_to_primer3.py`) parses FASTA headers to extract gene IDs and peak identifiers, formats each sequence as a Primer3 SEQUENCE_TEMPLATE, and appends configuration parameters from `primer3_settings.txt`.

**Peak Identifier Preservation**: In multi-peak mode, peak identifiers (e.g., `ENSG00000123456|peak_1`) are preserved in SEQUENCE_ID to enable downstream tracking of primers to specific peaks.

**Complexity**: O(G), where G = number of genes. Completes in seconds.

**Module**: `MAKE_PRIMER3_INPUT`

### Stage 5: Primer Design with Primer3

**Purpose**: Generate candidate primer pairs for each target sequence.

**Algorithm**: Executes `primer3_core` (version not specified) with settings from `config/primer3_settings.txt`. Primer3 evaluates thermodynamic properties (melting temperature, GC content, secondary structure) and designs optimal forward/reverse primer pairs.

**Default Settings** (from `config/primer3_settings.txt`):
- Product size: 75–150 bp
- Primer length: 18–25 bp (optimum: 20 bp)
- Melting temperature: 57–63°C (optimum: 60°C)
- GC content: 20–80%
- Max poly-X: 4 nucleotides
- Max self-complementarity: 12 (any), 8 (3′ end)
- Max primer-pair complementarity: 12 (any), 8 (3′ end)
- Number of primer pairs returned: 5

**Complexity**: O(G × L²), where G = genes, L = average template length. Typical runtime: 1–5 minutes for 1000 genes.

**Failure Modes**: Templates shorter than minimum product size (75 bp default) fail silently; reported in downstream summaries.

**Module**: `RUN_PRIMER3`

### Stage 6: cDNA Primer Extraction

**Purpose**: Extract primers appropriate for cDNA amplification from Primer3 output, ensuring strand compatibility.

**Algorithm**: R script (`bin/extract_cdna_primers.R`) parses Primer3 output and performs strand-aware selection:
1. For genes on the **plus strand**: select Primer3 LEFT primers (bind cDNA minus strand)
2. For genes on the **minus strand**: select Primer3 RIGHT primers (bind cDNA plus strand, reverse-complemented from genomic coordinates)
3. Preserve peak identifiers in multi-peak mode using composite keys (gene_id|peak_id)
4. Compute primer coordinates relative to genomic positions
5. Output TSV with columns: gene_id, peak_id (multi-peak mode), primer_index, sequence, tm, gc_percent, start, end, strand, length

**Key Logic**: The script uses a named vector lookup where keys are `gene_id|peak_id` (multi-peak) or `gene_id` (single-peak) to avoid overwriting strand information when multiple peaks exist per gene.

**Complexity**: O(P), where P = number of primer pairs. Completes in seconds.

**Failure Modes**: Genes with no Primer3 results produce empty primer lists (logged).

**Module**: `EXTRACT_CDNA_PRIMERS`

### Stage 7: Primer FASTA Conversion

**Purpose**: Convert primer TSV to FASTA format for transcriptome alignment, limiting to top N primers per gene/peak to manage computational cost.

**Algorithm**: R script (`bin/primers_to_fasta.R`) filters primers and formats headers:
1. In multi-peak mode: group primers by `gene_id|peak_id` and retain top `--max_primers_per_gene` per group
2. In single-peak mode: group by `gene_id` only
3. Generate FASTA headers: `>gene_id|peak_id|primer_index` (multi-peak) or `>gene_id|primer_index` (single-peak)

**Bug Fix (Implemented)**: Prior implementation grouped by `gene_id` only, limiting to 20 primers total in multi-peak mode. Current implementation uses composite keys to enforce per-peak limits.

**Complexity**: O(P), where P = primers. Completes in seconds.

**Module**: `PRIMERS_TO_FASTA`

### Stage 8: Transcriptome Specificity Screening (Optional)

**Purpose**: Align primers to transcriptome to detect off-target binding and cross-reactivity.

**Algorithm**: Bowtie2 (v2.x, version not specified) aligns primer sequences to a transcriptome index in local alignment mode with high sensitivity:
- `bowtie2 --local -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 -a`
- Reports all alignments (`-a`) to identify multi-mapping primers
- Local mode permits partial alignments (relevant for primers at exon junctions)

**Complexity**: O(P × T), where P = primers, T = transcriptome size. Runtime: 5–15 minutes for 1000 genes on 4 cores.

**Assumption**: Transcriptome index is built from the same annotation version (e.g., Gencode v38, Ensembl 109) as the input GTF to ensure transcript ID consistency.

**Module**: `ALIGN_PRIMERS_TRANSCRIPTOME`

### Stage 9: Alignment Analysis

**Purpose**: Summarize primer alignments, classify specificity, and calculate distance to 3′ ends.

**Algorithm**: Python script (`bin/analyze_primer_alignments.py`) parses alignment BAM:
1. For each primer, count alignments by mismatch category (0, 1, 2+ mismatches)
2. Calculate distance from alignment start to transcript 3′ end using GTF coordinates
3. Filter alignments > `--distance_threshold` bp from 3′ end (default: 400 bp)
4. Classify primer specificity:
   - **PERFECT**: 0 mismatches, single target gene
   - **GOOD**: 0 mismatches, multiple transcripts of same gene
   - **MODERATE**: 1 mismatch, or 0 mismatches to off-target gene
   - **POOR**: 2+ mismatches required for any alignment
   - **FAIL**: No valid alignments

**Key Data Structures**: Output TSV with columns: gene_id, primer_index, transcript_id, mismatches, distance_to_3prime, specificity_class

**Complexity**: O(A), where A = number of alignments. Typically < 1 minute for 1000 genes.

**Failure Modes**: Transcripts absent from GTF are labeled "unknown_transcript" but alignments are retained.

**Outputs**:
- `primer_alignment_report.tsv`: Per-primer specificity classification
- `primer_alignment_summary.tsv`: Detailed per-alignment records

**Module**: `ANALYZE_PRIMER_ALIGNMENTS`

### Stage 10: Best Primer Selection

**Purpose**: Select optimal primer per gene based on specificity and thermodynamic properties.

**Algorithm**: Python script (`bin/select_best_primer.py`) applies hierarchical filtering:
1. Filter primers by specificity class: PERFECT > GOOD > MODERATE (exclude POOR/FAIL)
2. Within each class, filter by Tm (58–62°C) and GC% (40–60%)
3. Rank by specificity class, then by ascending off-target count (1-mismatch hits)
4. Select top-ranked primer per gene

**Complexity**: O(P log P), where P = primers per gene. Sorting dominates; completes in seconds.

**Module**: `SELECT_BEST_PRIMERS`

### Stage 11: Isoform-Aware Optimization (Optional)

**Purpose**: Optimize primer selection to maximize isoform coverage with minimal primer sets.

**Algorithm**: Two distinct optimization strategies depending on operating mode:

#### Single-Peak Optimization (`optimize_primer_isoforms.py`)
Applies greedy set cover algorithm:
1. Filter alignment summary to primers in `best_primers.tsv` using composite keys `(gene_id, primer_index)` to prevent cross-contamination
2. Build primer-to-isoform mapping (distance ≤ 1000 bp from 3′ end, 0 mismatches)
3. Iteratively select primers covering the most uncovered isoforms until all reachable isoforms are covered or no primers remain
4. Report: primer_id, isoforms_targeted (count), target_transcripts (comma-separated list)

**Bug Fix (Implemented)**: Prior implementation filtered by `primer_index` alone, causing cross-contamination from other genes with the same primer index. Current implementation uses composite keys to ensure gene-specific filtering.

#### Multi-Peak Optimization (`optimize_primers_multipeak.py`)
Applies distance-weighted scoring:
1. For each gene, collect all primers from all peaks
2. Calculate composite score: `score = distance_weight × (1 - norm_distance) + isoform_weight × norm_isoforms + peak_rank_weight × (1 - norm_peak_rank)`
3. Select top `--primers_per_gene` primers by score
4. Report: primer_id, peak_id, distance_to_3prime, isoforms_targeted, peak_rank, composite_score

**Complexity**: O(P × I), where P = primers, I = isoforms per gene. Typical runtime: < 30 seconds for 1000 genes.

**Modules**: `OPTIMIZE_PRIMER_ISOFORMS`, `OPTIMIZE_PRIMERS_MULTIPEAK`

### Stage 12: Visualization (Optional)

**Purpose**: Generate publication-ready coverage plots with gene structure, peaks, and primer locations.

**Algorithm**: R script (`bin/MakePlots_new.R`) using ggplot2:
1. Import BigWig coverage for gene region ± flanking bases
2. Overlay gene structure (exons, introns, UTRs) from GTF
3. Annotate selected peak window and all MACS2 peaks within gene
4. Mark primer locations as directional arrows
5. Display QC metrics (peak score, q-value, exonic fraction) in subtitle

**Complexity**: O(G × L), where G = genes, L = average gene length. Parallelized (4 processes); 5–10 minutes for 1000 genes.

**Module**: `MAKEPLOTS_NEW`

## 6. Intermediate Artifacts and Data Flow

The pipeline produces a sequence of intermediate files enabling modular execution and debugging:

1. **BAM → MACS2**: `*_peaks.narrowPeak`, `*_summits.bed`, `*_peaks.xls`
2. **MACS2 → Peak Processing**: `selected_peaks.tsv`, `primer_targets.bed`, `primer_targets.fa`, `peaks_qc_summary.tsv`
3. **Peak Processing → Primer3**: `primer3_input.txt`
4. **Primer3 → Extraction**: `primer3_output.txt` → `cdna_primers.tsv`
5. **Extraction → Alignment**: `primers_for_alignment.fa`
6. **Alignment → Analysis**: `primers_alignment.bam` → `primer_alignment_report.tsv`, `primer_alignment_summary.tsv`
7. **Analysis → Selection**: `best_primers.tsv`
8. **Selection → Optimization**: `best_primers_optimal.tsv` (single-peak) or `optimized_primers_multipeak.tsv` (multi-peak)
9. **All stages → Visualization**: `plot_<gene_id>.pdf`

**Data Flow Diagram**:

```
BAM + GTF + Genes
    ↓
[BAM_INDEX] → BAM.bai
    ↓
[MACS2_CALLPEAK] → narrowPeak
    ↓
[PROCESS_MACS2_PEAKS] → selected_peaks.tsv, primer_targets.fa, peaks_qc_summary.tsv
    ↓
[MAKE_PRIMER3_INPUT] → primer3_input.txt
    ↓
[RUN_PRIMER3] → primer3_output.txt
    ↓
[EXTRACT_CDNA_PRIMERS] → cdna_primers.tsv
    ↓
[PRIMERS_TO_FASTA] → primers_for_alignment.fa
    ↓
[ALIGN_PRIMERS_TRANSCRIPTOME] → primers_alignment.bam (optional)
    ↓
[ANALYZE_PRIMER_ALIGNMENTS] → primer_alignment_report.tsv, primer_alignment_summary.tsv (optional)
    ↓
[SELECT_BEST_PRIMERS] → best_primers.tsv (optional)
    ↓
[OPTIMIZE_PRIMER_ISOFORMS / OPTIMIZE_PRIMERS_MULTIPEAK] → best_primers_optimal.tsv (optional)
    ↓
[SUMMARIZE_RESULTS] → pipeline_summary.txt

Parallel branch (if --makeplots):
BAM → [MEGADEPTH_BW] → coverage.bw → [MAKEPLOTS_NEW] → plot_*.pdf
```

## 7. External Dependencies

**Table D: External Dependencies**

| Tool | Version | Purpose | License |
|------|---------|---------|---------|
| Nextflow | ≥22.04.0 | Workflow orchestration | Apache 2.0 |
| MACS2 | 2.x (version not specified) | Peak calling | BSD-3-Clause |
| Primer3 | (version not specified) | Primer design | GPLv2 |
| Bowtie2 | 2.x (version not specified) | Transcriptome alignment | GPLv3 |
| samtools | 1.16 | BAM indexing | MIT |
| megadepth | (version not specified) | BigWig generation | MIT |
| R | ≥4.0 | Data processing, QC, visualization | GPL-3 |
| Bioconductor | ≥3.14 | GenomicRanges, rtracklayer, BSgenome packages | Artistic-2.0 |
| Python | ≥3.7 | Scripting (pandas, pysam) | PSF |
| Conda/Mamba | (version not specified) | Environment management | BSD-3-Clause |

**Reference Databases**:
- BSgenome packages (e.g., `BSgenome.Hsapiens.UCSC.hg38`): Used for sequence extraction
- Transcriptome indices: User-provided Bowtie2 indices matching GTF annotation

**Hardware Requirements**:
- Memory: 8 GB (MACS2, alignment analysis), 4 GB (most other processes), 3 GB per plot (when parallelized)
- Storage: ~500 MB per sample (BigWig, BAM index, intermediate files)
- Cores: 1 per process (plotting parallelized to 4 concurrent processes)

## 8. Outputs

**Table C: Outputs**

| Name | Format | Contents | Consumer | Typical Size |
|------|--------|----------|----------|--------------|
| `*_peaks.narrowPeak` | BED6+4 | MACS2-called peaks: chr, start, end, name, score, strand, fold_enrichment, -log10(p), -log10(q), summit_offset | Peak processing, visualization | 10–500 KB |
| `selected_peaks.tsv` | TSV | Selected peak per gene: gene_id, chr, start, end, strand, peak_rank, score, qvalue, pvalue, exonic_fraction | Primer design, plotting | 5–100 KB |
| `peaks_qc_summary.tsv` | TSV | Comprehensive QC metrics: gene_id, peak_count, selected_peak_score, qvalue, pvalue, exonic_fraction, trimming_applied, failure_reason | Quality assessment, plotting | 10–200 KB |
| `primer_targets.bed` | BED6 | Target regions: chr, start, end, gene_id, score, strand | Visualization | 5–50 KB |
| `primer_targets.fa` | FASTA | Target sequences for Primer3 input | Primer design | 10–100 KB |
| `primer3_output.txt` | Primer3 format | Raw Primer3 results: primer sequences, Tm, GC%, coordinates | cDNA extraction | 50–500 KB |
| `cdna_primers.tsv` | TSV | Strand-appropriate primers: gene_id, peak_id, primer_index, sequence, tm, gc_percent, start, end, strand, length | Alignment, selection | 20–200 KB |
| `primers_for_alignment.fa` | FASTA | Top N primers per gene/peak for alignment: header format `>gene_id\|peak_id\|primer_index` | Transcriptome alignment | 10–100 KB |
| `primer_alignment_report.tsv` | TSV | Specificity classification: gene_id, primer_index, specificity_class, perfect_matches, one_mismatch_hits, two_mismatch_hits | Primer selection | 10–50 KB |
| `primer_alignment_summary.tsv` | TSV | Detailed alignments: gene_id, primer_index, transcript_id, mismatches, distance_to_3prime, alignment_start, alignment_end | Isoform optimization | 50–500 KB |
| `best_primers.tsv` | TSV | Top primer per gene: gene_id, primer_index, sequence, tm, gc_percent, specificity_class, isoforms_aligned | Optimization (optional), final selection | 5–50 KB |
| `best_primers_optimal.tsv` | TSV | Isoform-optimized primers (single-peak mode): gene_id, primer_index, isoforms_targeted, target_transcripts | Experimental design | 5–50 KB |
| `optimized_primers_multipeak.tsv` | TSV | Multi-peak optimized primers: gene_id, primer_id, peak_id, distance_to_3prime, isoforms_targeted, peak_rank, composite_score | Experimental design | 10–100 KB |
| `*.bam.bw` | BigWig | Normalized coverage tracks: per-base coverage across genome | Visualization | 50–200 MB |
| `plot_<gene_id>.pdf` | PDF | Coverage plot with gene structure, peaks, primers: vector scalable | Publication, QC review | 50–200 KB each |
| `pipeline_summary.txt` | Text | Execution summary: gene counts, success/failure tallies, parameter settings | Run documentation | 1–5 KB |

**Typical Output Ranges**:
- Genes with selected peaks: 70–95% (depends on coverage depth and q-value threshold)
- Genes with valid primers: 65–90% (depends on sequence complexity and exon structure)
- Primers with PERFECT specificity: 40–70% (depends on transcriptome complexity)

## 9. Monitoring and Logging

### Provenance and Audit Trails

Nextflow automatically generates:
- **Execution report** (via `-with-report report.html`): Process-level resource usage, runtime, exit codes
- **Timeline** (via `-with-timeline timeline.html`): Gantt chart of task execution
- **DAG visualization** (via `-with-dag dag.pdf`): Workflow dependency graph
- **Trace file** (via `-with-trace trace.txt`): Per-task execution metadata

### Process-Level Logging

Each Nextflow process writes:
- Standard output (`.command.out`)
- Standard error (`.command.err`)
- Exit status (`.exitcode`)
- Execution script (`.command.sh`)

Logs are stored in `work/` subdirectories with cryptographic hashes (e.g., `work/a1/b2c3d4.../`).

### Warning Interpretation

**Common Warnings**:
- `Gene [ENSG...] not found in GTF`: Gene ID in input list absent from annotation; verify gene ID spelling and GTF version
- `No peaks pass q-value threshold for gene [ENSG...]`: Insufficient coverage or overly stringent threshold; inspect coverage plots or relax `--macs2_qvalue_threshold`
- `Peak trimmed to zero length for gene [ENSG...]`: Peak entirely outside exons after trimming; disable `--force_exonic_trimming` or investigate peak location
- `Primer3 returned no results for gene [ENSG...]`: Template sequence too short (< 75 bp default product size) or sequence complexity issues; check `primer_targets.fa`

### Pipeline Summary

The `SUMMARIZE_RESULTS` module generates `pipeline_summary.txt` with:
- Total genes processed
- Genes with selected peaks (count, percentage)
- Genes with primers designed (count, percentage)
- Genes with uniquely aligned primers (count, percentage, if transcriptome QC enabled)
- Per-gene failure reasons (no peaks, QC fail, primer design fail, alignment fail)

## 10. Performance and Scalability

### Runtime Estimates

For a typical human RNA-seq dataset (30M reads, 1000 target genes) on a local machine (4 cores, 16 GB RAM):

| Stage | Runtime |
|-------|---------|
| BAM indexing | 30 s |
| MACS2 peak calling | 10 min |
| BigWig generation | 3 min |
| Peak processing | 1 min |
| Primer3 design | 3 min |
| cDNA extraction | 10 s |
| Transcriptome alignment | 10 min |
| Alignment analysis | 30 s |
| Optimization | 20 s |
| Visualization (1000 genes, 4 parallel) | 8 min |
| **Total** | **35 min** |

### Scalability

**Batch vs. Streaming**: The pipeline operates in batch mode (processes all genes per stage before proceeding). Assumption: Streaming would require refactoring to process genes independently through all stages, which is not currently supported.

**Parallelization**:
- MACS2, Bowtie2, and plotting modules are single-process parallel (multi-threaded within the tool)
- Plotting is parallelized across genes (`maxForks = 4` in `nextflow.config`)
- All other stages are single-threaded but parallelizable by gene via Nextflow scatter-gather patterns (not implemented by default)

**Memory Scaling**:
- MACS2: O(N) where N = reads; ~8 GB for 50M reads
- Peak processing: O(P + E) where P = peaks, E = exons; ~4 GB for human genome
- Plotting: O(L) where L = gene length; ~3 GB per plot (BigWig caching)

**Large-Scale Runs** (10,000 genes):
- Increase `maxForks` for plotting: `process.withName:MAKEPLOTS_NEW.maxForks = 16`
- Use HPC scheduler (PBS, SLURM): `-profile pbs`
- Split gene list into batches and run independently, then merge outputs

## 11. Determinism and Reproducibility

### Seeds and Random Number Generators

Assumption: None of the pipeline tools use random number generators; all processes are deterministic given identical inputs.

### Checkpoints and Caching

Nextflow's `-resume` flag enables caching: completed tasks are skipped if inputs, scripts, and parameters are unchanged. Cached results are stored in `work/` directories.

### Environment Pinning

- **Conda environments**: Each process specifies a Conda environment via `conda` directives in module definitions. Environments are defined in `env/*.yml` files with pinned package versions (e.g., `macs2=2.2.7.1`, `samtools=1.16`).
- **Container images**: Pipeline supports Singularity containers (specified via `container` directives, not enabled by default). Assumption: Container images are built from Conda environments and versioned.

### Exact Reproduction Steps

To reproduce a run:

1. **Install Nextflow**: `curl -s https://get.nextflow.io | bash`
2. **Clone repository**: `git clone <repo_url> && cd PeakPrime`
3. **Prepare inputs**: Ensure BAM, GTF, gene list, and transcriptome index are identical (verify file checksums)
4. **Run with resume**: `nextflow run main.nf --bam <bam> --gtf <gtf> --genes <genes> --transcriptome_index <index> --outdir results/ -profile local -resume -with-report -with-timeline -with-trace`
5. **Verify outputs**: Compare checksums of key outputs (`selected_peaks.tsv`, `best_primers.tsv`) against reference run

## 12. Error Handling and Edge Cases

### Malformed Inputs

| Input | Malformation | Behavior |
|-------|--------------|----------|
| BAM | Unsorted | MACS2 fails with exit code 1; Nextflow halts |
| BAM | Unindexed | Pipeline creates index automatically; proceeds |
| GTF | Invalid format (e.g., GFF3) | R import fails; error message: "GTF parsing failed" |
| GTF | Missing gene_id attributes | GenomicRanges import fails; Nextflow halts |
| Gene list | Invalid IDs (not `ENSG[0-9]+`) | Genes skipped with warnings; proceeds with valid IDs |
| Gene list | Empty file | No genes processed; pipeline exits cleanly with empty outputs |
| FASTA | Unindexed | Bioconductor auto-indexes on first use; proceeds |
| Transcriptome index | Missing files | Bowtie2 fails with "index not found"; Nextflow halts |

### Missing Metadata

- **Strand information absent from GTF**: Peaks cannot be processed; error: "Strand column required"
- **Transcript IDs missing from GTF**: Isoform optimization skipped with warning; best primers output produced

### Out-of-Range Values

- **Negative coordinates in BED files**: GenomicRanges throws error; Nextflow halts
- **`--distance_threshold < 0`**: Pipeline validates at start; error: "Parameter must be positive"
- **`--peak_rank = 0`**: Treated as invalid; error: "peak_rank must be ≥ 1"

### Retries and Fallbacks

- **Primer3 failures**: No retry logic; genes with failed Primer3 runs are logged in QC summary with `failure_reason = "primer_design_fail"`
- **Alignment failures**: No retry logic; genes without alignments proceed but skip optimization
- **Sequence extraction failures** (BSgenome): Fallback to FASTA if `--fasta` provided; otherwise fails

### Edge Cases

- **Genes with no exons in GTF**: Exonic fraction = 0; discarded if `--min_exonic_fraction` > 0
- **Peaks spanning multiple genes**: Peak assigned to gene with maximum overlap
- **Primers aligning to multiple genes**: Classified as MODERATE or POOR specificity depending on mismatch count
- **Multi-peak mode with 1 peak per gene**: Behaves identically to single-peak mode
- **Distance mode with unavailable transcripts**: Falls back to fetching from Ensembl if `--transcript_fasta` not provided

## 13. Security and Privacy Considerations

### Data Minimization

The pipeline processes genomic coordinates and sequences only; no personally identifiable information (PII) is expected in standard inputs. Assumption: BAM files are derived from cell lines, model organisms, or anonymized clinical samples.

### PII Handling

If BAM headers contain sample identifiers or metadata:
- Recommendation: Strip headers using `samtools reheader` before input
- Pipeline does not propagate BAM headers to outputs (only gene-level summaries produced)

### Sandboxing

- Conda environments isolate tool dependencies, preventing system-wide package conflicts
- Nextflow runs processes in isolated `work/` directories with restricted write access
- Assumption: Execution on shared HPC systems uses user-level permissions; no elevated privileges required

### Compliance Constraints

Assumption: Users are responsible for ensuring compliance with data governance policies (e.g., GDPR, HIPAA) when processing human genomic data. Pipeline does not implement encryption, access controls, or audit logging beyond Nextflow's standard provenance tracking.

---

## Worked Example

### Command-Line Invocation (Peak-Based Single-Peak Mode with Optimization)

```bash
nextflow run main.nf \
  --bam sample_RNA_seq.bam \
  --gtf gencode.v38.annotation.gtf \
  --genes target_genes.txt \
  --genome_package BSgenome.Hsapiens.UCSC.hg38 \
  --macs2_qvalue_threshold 0.05 \
  --peak_rank 1 \
  --force_exonic_trimming true \
  --transcriptome_index transcriptome_GRCh38.109 \
  --transcript_mapping resources/transcript_gene_mapping.GRCh38.109.tsv \
  --max_primers_per_gene 20 \
  --distance_threshold 400 \
  --makeplots \
  --outdir results_optimized/ \
  -profile local \
  -resume \
  -with-report report.html \
  -with-timeline timeline.html
```

**Expected Outputs** (1000 genes, 30M reads):

```
results_optimized/
├── sample_RNA_seq_peaks.narrowPeak          (250 KB, 1823 peaks)
├── selected_peaks.tsv                        (45 KB, 872 genes with peaks)
├── peaks_qc_summary.tsv                      (78 KB, 1000 genes with QC metrics)
├── primer_targets.bed                        (23 KB, 872 regions)
├── primer_targets.fa                         (67 KB, 872 sequences)
├── primer3_output.txt                        (412 KB, 4360 primer pairs)
├── cdna_primers.tsv                          (156 KB, 4360 primers)
├── primers_for_alignment.fa                  (89 KB, 17440 primers, 20/gene)
├── primers_alignment.bam                     (2.3 MB, 87231 alignments)
├── primer_alignment_report.tsv               (34 KB, specificity classes)
├── primer_alignment_summary.tsv              (278 KB, detailed alignments)
├── best_primers.tsv                          (28 KB, 872 primers)
├── best_primers_optimal.tsv                  (31 KB, 872 primers with isoforms)
├── sample_RNA_seq.bam.bw                     (147 MB, coverage track)
├── plot_ENSG00000000001.pdf                  (120 KB)
├── plot_ENSG00000000002.pdf                  (95 KB)
├── ...                                       (872 plots total)
└── pipeline_summary.txt                      (2 KB, execution summary)
```

**Summary Statistics** (from `pipeline_summary.txt`):

```
Total genes processed: 1000
Genes with selected peaks: 872 (87.2%)
Genes with primers designed: 872 (87.2%)
Genes with uniquely aligned primers: 623 (62.3%)

Failure breakdown:
- No peaks called: 98 (9.8%)
- QC fail (exonic fraction): 14 (1.4%)
- Primer design fail: 0 (0.0%)
- Alignment fail (no PERFECT/GOOD): 249 (28.6% of primers designed)
```

### Multi-Peak Mode with Optimization

```bash
nextflow run main.nf \
  --bam sample_RNA_seq.bam \
  --gtf gencode.v38.annotation.gtf \
  --genes target_genes.txt \
  --select_all_peaks \
  --max_peaks_per_gene 5 \
  --optimize_multipeak \
  --primers_per_gene 3 \
  --distance_weight 0.5 \
  --isoform_weight 0.3 \
  --peak_rank_weight 0.2 \
  --transcriptome_index transcriptome_GRCh38.109 \
  --outdir results_multipeak/ \
  -profile local -resume
```

**Expected Additional Outputs**:

```
results_multipeak/
├── selected_peaks.tsv                        (134 KB, 2314 peaks, avg 2.7 peaks/gene)
├── optimized_primers_multipeak.tsv           (87 KB, 2616 primers, top 3 per gene)
└── ...
```

---

## Summary

PeakPrime integrates statistical peak detection with automated primer design to enable reproducible, high-throughput targeted enrichment workflows for 3′-end RNA-seq. The pipeline supports flexible operating modes (single-peak, multi-peak, distance-based) with optional isoform-aware optimization, comprehensive QC reporting, and publication-ready visualization. All processing stages are deterministic, portable via Conda environments, and fully documented for reproducibility. The modular architecture enables incremental execution via Nextflow's caching, and execution logs provide complete provenance tracking for audit and debugging.
