# PeakPrime pipeline — Detailed steps

This document describes the main steps performed by the PeakPrime Nextflow pipeline and the postprocessing used to prepare data for the standalone Shiny app.

## Contract (what this pipeline does)
- **Inputs**: One or more BAM files (already aligned, spliced alignments required for RNA-seq), gene/target list, a GTF annotation (optional if auto-detected).
- **Outputs**: narrowPeak peak calls, primer candidates (TSV/FASTA), primer alignment results, RDS bundles for the Shiny app.
- **Success criteria**: Generated `qc_data.rds`, `coverage_index.rds`, `gtf_data.rds`, and optional `primer_alignment_summary.tsv`.
- **Note**: The pipeline does NOT perform alignment. You must provide pre-aligned BAM files.

## Step-by-step pipeline

### 1) Input discovery
- **What**: Collect user inputs and raw files.
- **Required files**: 
  - BAM file (`*.bam` + `.bam.bai`) — must be sorted and indexed
  - Gene list (text file with one Ensembl gene ID per line)
- **Optional files**: 
  - GTF annotation (auto-detected or specified via `--gtf`)
  - Parameter config (`params.config`)
- **Notes**: The pipeline requires pre-aligned BAM files. For RNA-seq, spliced alignments are expected. The pipeline auto-detects a GTF in the results directory or uses a specified path.

### 2) Coverage generation
- **What**: Create per-base coverage or binned coverage suitable for plotting and peak calling.
- **Tools**: deepTools `bamCoverage`, `bedtools genomecov`, or custom scripts
- **Key outputs**: `*.bw` bigWig files (preferred for app performance)
- **Example**: `bamCoverage -b sample.bam -o sample.bw --normalizeUsing CPM --binSize 10`
- **Notes**: Use consistent normalization (CPM/RPKM) if comparing samples.

### 3) Peak calling
- **What**: Call peaks per-sample (or combined) using MACS2 or another peak caller.
- Tools: MACS2
- Key outputs: `*_peaks.narrowPeak`, `selected_peaks.tsv`
- Parameters: `--shift`, `--extsize`, and q-value thresholds are usually configurable.
- Example: `macs2 callpeak -t sample.bam -f BAM -g hs -q 0.05 --outdir peaks/`

### 4) Peak quality control & selection
- **What**: Score peaks and generate `peaks_qc_summary.tsv` with metrics used by the app to pick candidate regions.
- **Tools**: R scripts (`pick_peaks.R`, `pick_peaks_span.R`) that compute QC metrics (peak width, fold change, summit, overlap with transcripts)
- **Key outputs**: `peaks_qc_summary.tsv`
- **Metrics**: Peak score, fold enrichment, summit position, distance to transcript start/end

### 5) Primer design
- **What**: Run primer3 (or `primer3_core`) to design candidate primers for each selected region.
- **Tools**: primer3, primer3-py wrappers, or custom Python/R scripts (`fasta_to_primer3.py`)
- **Key outputs**: `primers.tsv`, `primers.fasta` (for alignment and review)
- **Example**: `primer3_core < input.txt > output.txt` (using config from `config/primer3_settings.txt`)
- **Tips**: Keep consistent settings via `config/primer3_settings.txt` and log primer3 runtime failures. This generates multiple candidate primers per region.

### 6) Primer alignment (specificity check)
- **What**: Align designed primers back to the transcriptome/genome to check for off-targets and mismatches.
- **Tools**: Bowtie, Bowtie2, or `blastn` (for sensitive checks)
- **Key outputs**: Raw alignment files (SAM/BAM format), per-sample or combined alignment TSVs
- **Example**: `bowtie -f -v 2 -a -p 4 transcriptome_index primers.fasta > primers.alignments.sam`
- **Notes**: Index and align against the same reference used for genome coordinates. Use `-v 2` to allow up to 2 mismatches. This step is critical for identifying off-target binding.

### 7) Alignment summary & metrics
- **What**: Summarize alignments to provide stats like number of perfect matches, off-target counts, and mapping positions.
- **Tools**: R/Python scripts (`analyze_primer_alignments.py`, `analyze_primer_alignments.R`, `make_alignment_summary.py`)
- **Key outputs**: `primer_alignment_summary.tsv`, used by the Shiny app's Primer Alignment tab
- **Metrics**: Perfect matches, 1-mismatch hits, 2-mismatch hits, total alignments, off-target positions

### 8) Primer filtering & selection (best primer)
- **What**: Select the best primer per gene based on alignment specificity. Ideal primers map without mismatches to only one gene (the target gene).
- **Tools**: R/Python scripts (e.g., `select_best_primer.py`, `summarize_primers.py`)
- **Key outputs**: `selected_primers.tsv` with final primer choices per gene
- **Selection criteria**: 
  - **Primary**: Perfect match (0 mismatches) to target gene only
  - **Secondary**: Tm 58-62°C, GC% 40-60%, amplicon 80-150bp
  - **Filters**: Reject primers with perfect matches to multiple genes (off-targets)
  - **Ranking**: Among valid candidates, prefer primers with no 1-mismatch hits
- **Notes**: This is the final selection step that produces the recommended primers for experimental validation.

### 9) Postprocessing for app (preprocess_for_standalone.R)
- **What**: Convert pipeline text outputs into fast R objects (RDS) and create a `data_manifest.rds` summarizing what's available.
- **Key files produced**:
  - `qc_data.rds` — data.frame with gene-level QC and coordinates
  - `gtf_data.rds` — processed GTF for gene structures and isoforms
  - `coverage_index.rds` — optimized data structure for fast coverage extraction
  - `primer_alignment_summary.rds` — processed alignment data (if available)
  - `data_manifest.rds` — meta-info and paths
- **Command**: `Rscript preprocess_for_standalone.R results/your_run` or in R:
  ```r
  source('preprocess_for_standalone.R')
  preprocess_peakprime_standalone('results/your_run')
  ```
- **Requirements**: R packages: `data.table`, `rtracklayer`, `GenomicRanges`, `IRanges`

### 10) Run Shiny app (app_standalone.R)
- **What**: The app loads the RDS bundle and provides interactive exploration, plotting, and downloads.
- **Command**: `Rscript app_standalone.R results/your_run` or in R:
  ```r
  source('app_standalone.R')  # Enter directory when prompted
  runApp(peakprime_app)
  ```
- **Features**: Gene plot viewer, primer alignment tab, CSV/PNG export, multi-peak visualization

## Important parameters and config files
- `config/primer3_settings.txt` — primer3 tuning parameters
- Nextflow `params.config` — sample lists, genome references, and alignment options

## Quality checks & edge cases

| Issue | Impact | Solution |
|-------|--------|----------|
| Missing GTF | Gene structure plots limited | App still works for coverage and peaks; specify GTF path in preprocessing |
| Missing coverage bigWig | Plotting will fail | Ensure `bamCoverage` or equivalent was run; check `.bw` files exist |
| Primer design failures | No primers for some genes | Log primer3 errors and inspect `primers.tsv` for missing amplicons |
| Primer alignment missing | Alignment tab empty | Ensure alignment step ran and `primer_alignment_summary.tsv` exists |
| No peaks called | QC summary empty | Check MACS2 parameters (q-value may be too stringent), verify BAM quality |

## Example Nextflow command
```bash
# Basic run with BAM input
nextflow run main.nf --bam sample.bam --genes targets.txt --gtf annotation.gtf \
  --outdir results/my_run

# With transcriptome alignment QC
nextflow run main.nf --bam sample.bam --genes targets.txt \
  --transcriptome_index transcriptome_idx --transcriptome_fasta transcripts.fa \
  --outdir results/my_run --with-report report.html

# With plotting enabled
nextflow run main.nf --bam sample.bam --genes targets.txt \
  --makeplots --outdir results/my_run
```

**Note**: If you have FASTQ files, you must align them first using a splice-aware aligner like STAR or HISAT2 before running PeakPrime.

## Troubleshooting

### Issue: `preprocess_for_standalone.R` reports missing files
**Check expected outputs:**
```r
list.files('results/your_run', recursive = TRUE)
# Expected: peaks_qc_summary.tsv, *.bw, *_peaks.narrowPeak, primers.tsv
```
**Verify RDS files were created:**
```r
list.files('results/your_run', pattern = '\\.rds$')
# Expected: qc_data.rds, gtf_data.rds, coverage_index.rds, data_manifest.rds
```

### Issue: Primer alignment stats are unexpectedly empty
**Check alignment step ran:**
```bash
ls -lh results/your_run/*alignment*.tsv
```
**Verify aligner index matches reference:**
- Ensure the transcriptome index used for Bowtie matches the GTF annotation version
- Check aligner logs for errors (e.g., `bowtie.log`, `alignment.stderr`)

### Issue: App shows "GTF data: ❌ Missing"
**Re-run preprocessing with explicit GTF path:**
```r
source('preprocess_for_standalone.R')
preprocess_peakprime_standalone('results/your_run', '/path/to/annotation.gtf')
```

## Next steps / improvements
- Add optional BLAST-based primer specificity checks for higher sensitivity
- Produce a human-readable export (Excel) combining primer metrics and alignment summaries
- Add automated sanity checks in the Nextflow run to fail early on missing inputs

---

For edits or additions, tell me what specific step you'd like expanded (e.g., primer3 parameters, MACS2 options, or example Nextflow config snippets).