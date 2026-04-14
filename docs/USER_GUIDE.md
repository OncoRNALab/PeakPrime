# PeakPrime — Lab User Guide

> A step-by-step guide for designing primers with the PeakPrime pipeline.  
> Written for wet-lab users; no bioinformatics experience assumed.

---

## Table of contents

1. [What does this pipeline do?](#1-what-does-this-pipeline-do)
2. [Installation and requirements](#2-installation-and-requirements)
3. [Two ways to design primers](#3-two-ways-to-design-primers)
4. [Preparing your input files](#4-preparing-your-input-files)
5. [Running the pipeline](#5-running-the-pipeline)
6. [Key parameters you can tune](#6-key-parameters-you-can-tune)
7. [Single-peak vs multi-peak mode](#7-single-peak-vs-multi-peak-mode)
8. [Primer optimisation — what it means and why it matters](#8-primer-optimisation--what-it-means-and-why-it-matters)
9. [Understanding the output files](#9-understanding-the-output-files)
10. [Troubleshooting — what went wrong?](#10-troubleshooting--what-went-wrong)

---

## 1. What does this pipeline do?

PeakPrime automatically designs primers for your genes of interest.  
In plain terms, it does the following:

1. **Finds where in the gene the RNA is most abundantly expressed** (the "peak") using your sequencing data, or uses a fixed window near the end of the transcript.
2. **Designs primer sequences** in that region using Primer3, a standard tool for primer design.
3. **Checks whether the designed primers are specific** by aligning them against all known transcripts in the human transcriptome. Primers that hit many unintended transcripts are flagged.
4. **Selects and ranks the best primers** — those that sit close to the end of the transcript, match as many transcript versions (isoforms) of your gene as possible, and avoid hitting other genes.

---

## 2. Installation and requirements

### Software you need on your system

| Software | Purpose |
|---|---|
| **Nextflow** (≥ 22.04.0) | Runs the pipeline |
| **Conda** or **Mamba** | Manages all other software automatically |

Nextflow and Conda/Mamba are the only two things you need to install manually. All other tools (MACS2, Primer3, Bowtie2, R, Python, etc.) are installed automatically by the pipeline the first time it runs.

### Cloning the pipeline

```bash
git clone <repository_url>
cd PeakPrime
```

### Running profiles

The pipeline supports two execution profiles:

- **`local`** — runs on your current computer or login node.
- **`pbs`** — submits jobs to a PBS/Torque HPC cluster (recommended for large datasets).

Add `-profile local` or `-profile pbs` to your run command.

---

## 3. Two ways to design primers

### Mode A — Coverage-informed mode (default)

**When to use it:** You have RNA-seq data for your samples, and you want primers targeting the region where your gene is most highly expressed.

**How it works:**  
The pipeline calls MACS2 to detect peaks of RNA signal in your sequencing data. It then focuses primer design in that peak region. This ensures your primer targets an area that is genuinely expressed in your biological context.

**What you need:**
- A BAM file (your aligned RNA-seq reads)
- A GTF annotation file (gene coordinates for your genome)
- A text file listing the genes you want to design primers for (one Ensembl gene ID per line)

### Mode B — Distance mode

**When to use it:** You do not have RNA-seq data, or you want to design primers based purely on gene structure — specifically, at a fixed distance from the end of the transcript.

**How it works:**  
Instead of relying on sequencing peaks, the pipeline fetches the standard reference transcript for each gene (called the MANE transcript) from Ensembl, and then designs primers in a window of a specified length at the 3′ end of the transcript.

**What you need:**
- Your gene list (one Ensembl gene ID per line), OR a FASTA file with your transcript sequences
- A `--template_length` value (number of bases from the end of the transcript to use as the design window)

---

## 4. Preparing your input files

### Gene list file

A plain text file, one Ensembl gene ID per line:

```
ENSG00000141510
ENSG00000012048
ENSG00000139618
```

### BAM file (coverage mode only)

Your RNA-seq reads aligned to the genome. It must be sorted and indexed (a `.bai` index file must exist alongside it).

### GTF file (coverage mode only)

The gene annotation for your genome build (e.g. Ensembl GRCh38 release 109). Must match the genome build used to align your BAM.

### Bowtie2 transcriptome index (optional but recommended)

A pre-built Bowtie2 index of the human transcriptome. This is needed for the primer specificity checking step. Without it, the pipeline will design primers but will not evaluate their specificity.

---

## 5. Running the pipeline

### Coverage-informed mode (basic)

```bash
nextflow run main.nf \
  -profile local \
  --bam   your_reads.bam \
  --gtf   annotation.gtf \
  --genes gene_list.txt \
  --outdir results/
```

### Coverage-informed mode (with specificity checking)

```bash
nextflow run main.nf \
  -profile local \
  --bam   your_reads.bam \
  --gtf   annotation.gtf \
  --genes gene_list.txt \
  --transcriptome_index /path/to/transcriptome_idx \
  --outdir results/
```

### Distance mode

```bash
nextflow run main.nf \
  -profile local \
  --distance_mode \
  --genes gene_list.txt \
  --template_length 500 \
  --transcriptome_index /path/to/transcriptome_idx \
  --outdir results/
```

### HPC cluster (PBS)

Just swap `-profile local` for `-profile pbs`. Nextflow will automatically submit jobs to the queue.

---

## 6. Key parameters you can tune

### MACS2 peak-calling settings (coverage mode only)

These control how strict the pipeline is when calling peaks in your RNA-seq data.

| Parameter | Default | What it does |
|---|---|---|
| `--macs2_qvalue_threshold` | `0.05` | Statistical significance threshold for calling a peak. Lower = stricter (fewer but more confident peaks). Raise to `0.1` or higher if too few peaks are found. |
| `--macs2_min_peak_score` | `0` | Minimum peak score. Usually left at 0. |
| `--macs2_extsize` | auto | Fragment size for extending reads. Set manually if MACS2 cannot estimate it automatically from your data. |
| `--peak_selection_metric` | `score` | How to rank peaks when a gene has more than one. Options: `score` (peak height) or `qvalue` (statistical significance). |
| `--peak_rank` | `1` | Which peak to use per gene. `1` = the best-ranked peak, `2` = second-best, etc. Useful if you want to target a less dominant peak. |

### Exon trimming settings (coverage mode only)

Peaks from your sequencing data are genomic regions that may overlap both exons and introns. Since primers must work on RNA (which has introns removed), the pipeline can trim the design window to stay within exons.

| Parameter | Default | What it does |
|---|---|---|
| `--force_exonic_trimming` | `true` | **Recommended.** Trims the design window to the largest overlapping exon/UTR region. Ensures primers designed on purely exonic sequence. |
| `--min_trimmed_length` | `30` | Minimum length (in base pairs) of the trimmed region. If trimming leaves less than this, the peak is discarded. Increase if you want to guarantee a larger design window. |
| `--trim_to_exon` | `false` | A softer alternative to `force_exonic_trimming`. Only trims to the exon containing the peak summit. Ignored if `force_exonic_trimming` is on. |
| `--min_exonic_fraction` | off | Discard peaks where less than this fraction (0–1) of the window overlaps exons. Only active when `force_exonic_trimming` is off. |

### Primer3 settings

Primer3 is the tool that actually designs the primer sequences. Its settings are stored in a separate file:

```
config/primer3_settings.txt
```

You can edit this file directly, or supply your own with `--primer3_settings /path/to/your_settings.txt`.

**Default Primer3 settings:**

| Property | Min | Optimal | Max |
|---|---|---|---|
| Primer length | 18 bp | 20 bp | 22 bp |
| Melting temperature (Tm) | 59 °C | 60 °C | 61 °C |
| GC content | 35% | 50% | 65% |
| Max homopolymer run | — | — | 4 bp |
| GC clamp at 3′ end | — | 1 | — |
| Number of primers returned per gene | — | 20 | — |

> **Note:** The pipeline only uses **LEFT primers** (those matching the mRNA sequence in the 5′ → 3′ direction). Right primers are not designed.

### Primer specificity (alignment) settings

These control how the pipeline judges whether a primer is specific to your gene.

| Parameter | Default | What it does |
|---|---|---|
| `--transcriptome_index` | none | Path to a Bowtie2 transcriptome index. **Required** for specificity checking. |
| `--max_mismatches` | `3` | A primer is considered to match a transcript if it aligns with this many or fewer mismatches. Lower = more strict. |
| `--distance_threshold` | `1000` | Only count alignments within this many bases from the 3′ end of a transcript. This focuses specificity analysis on the region relevant to 3′-targeted RNA-seq. |
| `--max_primers_per_gene` | `20` | Maximum number of primers sent for alignment QC per gene. |

### Distance mode settings

| Parameter | Default | What it does |
|---|---|---|
| `--template_length` | required | How many bases from the 3′ end of the transcript to use as the design window. Larger = more design options. |
| `--transcript_fasta` | auto | If provided, uses these transcript sequences instead of fetching from Ensembl. |

---

## 7. Single-peak vs multi-peak mode

By default, the pipeline selects **one peak per gene** (the highest-scoring one) and designs primers in that window. This is called **single-peak mode**.

### Multi-peak mode

If you add `--select_all_peaks`, the pipeline instead selects **all significant peaks per gene** and designs primers in each one. This gives you:

- More primer candidates per gene
- Primers targeting different regions of the gene (useful for genes with complex expression patterns)
- A `peak_id` column in the output files identifying which peak each primer belongs to

**Limiting the number of peaks per gene:**

```bash
--select_all_peaks --max_peaks_per_gene 3
```

This takes the top 3 peaks per gene (by score or q-value, depending on `--peak_selection_metric`).

**Multi-peak optimisation:**

Adding `--optimize_multipeak` activates a scoring step that ranks primers across all peaks of a gene, taking into account:
- How close the primer is to the 3′ end of the transcript
- How many transcript variants (isoforms) it covers
- The quality ranking of the peak it was designed from

You can control the importance of each factor:

```bash
--distance_weight 0.5   # weight for 3' proximity (default 0.5)
--isoform_weight 0.3    # weight for isoform coverage (default 0.3)
--peak_rank_weight 0.2  # weight for peak quality (default 0.2)
--primers_per_gene 3    # how many final primers to select per gene (default 3)
```

---

## 8. Primer optimisation — what it means and why it matters

Once the pipeline has aligned all primers to the transcriptome, it runs an **optimisation step** to select the best primer for each gene.

### Why is this needed?

A gene often exists in several slightly different forms called **isoforms** (different transcripts produced from the same gene). A good primer should ideally match as many isoforms as possible, so that the assay works regardless of which isoform is expressed.

### What the optimiser does

The optimiser goes through all primers that passed the specificity filters and for each gene selects the primer that:

1. **Covers the most isoforms** — aligns to the highest number of transcript variants of that gene.
2. **Sits close to the 3′ end** — important for 3′-biased RNA-seq protocols where reads are enriched toward the end of the transcript.
3. **Is specific** — does not cross-react with transcripts from other genes.

The result is stored in `best_primers_optimal.tsv` (single-peak or distance mode) or `optimized_primers_multipeak.tsv` (multi-peak mode).

> If you did not provide a `--transcriptome_index`, the optimisation step is skipped entirely, and the pipeline only outputs the raw Primer3 designs.

---

## 9. Understanding the output files

All output files are saved in the directory you set with `--outdir` (default: `results/`).

### Key files at a glance

| File | What it tells you |
|---|---|
| `processed_peaks/peaks_qc_summary.tsv` | Per-gene peak calling QC — why each gene passed or failed |
| `processed_peaks/selected_peaks.tsv` | The final peak window chosen for each gene |
| `processed_peaks/selected_peaks.fa` | The DNA sequences fed to Primer3 |
| `primer3_output.txt` | Raw Primer3 output — useful for debugging |
| `cdna_primers.tsv` | All primers returned by Primer3, before specificity checking |
| `primer_alignment_report.tsv` | One row per primer — how many times it aligned to the transcriptome |
| `primer_alignment_summary.tsv` | One row per alignment hit — where exactly each primer mapped |
| `best_primers.tsv` | Primers that passed the specificity filters |
| `best_primers_optimal.tsv` | The final, isoform-optimised primer for each gene |
| `optimized_primers/optimized_primers_multipeak.tsv` | Multi-peak mode: final selected primers |
| `alignment_stats.txt` | Overall Bowtie2 alignment statistics |
| `macs2_peaks/` | Folder containing all raw MACS2 output files |

### The two alignment output files in detail

**`primer_alignment_report.tsv`** — the complete inventory

Every primer Primer3 designed appears here. The key column is `num_alignments`:
- `0` = the primer did not map to any transcript in the transcriptome (see [troubleshooting](#10-troubleshooting--what-went-wrong))
- `1` = the primer maps to exactly one location (ideal for specific primers)
- `>1` = the primer maps to multiple locations (may indicate cross-reactivity)

**`primer_alignment_summary.tsv`** — the detailed alignment view

Only primers that successfully mapped to at least one transcript appear here. Each row is one alignment hit. Key columns:
- `aligned_transcript` — which transcript it mapped to
- `aligned_gene_name` — which gene that transcript belongs to
- `mismatches` — number of sequence mismatches (0 = perfect match)
- `distance_to_end` — how many bases from the 3′ end of the transcript
- `alignment_strand` — which strand it mapped to

---

## 10. Troubleshooting — what went wrong?

### How to check if MACS2 failed to find peaks

**File to open:** `macs2_peaks/<sample>_peaks.narrowPeak`

If this file is empty (no rows of data), MACS2 found no peaks in your data. This usually means:
- The BAM file has very low coverage for your genes of interest
- The `--macs2_qvalue_threshold` is too strict → try increasing it (e.g. `--macs2_qvalue_threshold 0.1`)
- The BAM is not from RNA-seq or is aligned to a different genome build

---

### How to check if peaks didn't pass the filtering criteria

**File to open:** `processed_peaks/peaks_qc_summary.tsv`

Open this file in Excel or a similar program. Each row is one of your target genes. Look at these columns:

| Column | What it means |
|---|---|
| `final_selection` | `TRUE` = this gene has a valid peak and will get primers; `FALSE` = failed |
| `failure_reason` | Plain-text explanation of why the gene was excluded |
| `has_macs2_peaks` | Were any MACS2 peaks found for this gene at all? |
| `has_significant_peaks` | Did any peaks pass the q-value threshold? |
| `has_overlapping_peaks` | Did any significant peaks overlap the gene body? |
| `passes_exonic_filter` | Did the peak region contain enough exonic sequence? |
| `exonic_fraction` | What fraction (0–1) of the peak overlaps exons |

**Common failure reasons and what to do:**

| Failure reason | What it means | What to try |
|---|---|---|
| `No MACS2 peaks detected` | No RNA signal peaks found for this gene | Increase `--macs2_qvalue_threshold`; check the BAM has coverage for this gene |
| `Failed significance filter (q>…)` | Peaks found but below significance threshold | Increase `--macs2_qvalue_threshold` |
| `Significant peaks found but no overlap with gene boundaries` | Peaks exist nearby but not within the gene coordinates in the GTF | Check GTF matches your BAM genome build |
| `Failed exonic fraction filter` | Peak is mostly in intronic sequence | Lower `--min_exonic_fraction` or switch to `--force_exonic_trimming` |
| `Failed forced exonic trimming … too short` | The exonic part of the peak is shorter than `--min_trimmed_length` | Lower `--min_trimmed_length` (e.g. `--min_trimmed_length 100`) |

---

### How to check if Primer3 failed to design primers

Even when a peak is found, Primer3 may be unable to design primers if the region has unfavourable sequence (e.g. too many repeated bases, very low GC content, etc.).

**File to open:** `cdna_primers.tsv`

If a gene is missing from this file entirely, Primer3 found no valid primers for it.

**For more detail:** `primer3_output.txt`  
Search for your gene ID. If you see `PRIMER_PAIR_NUM_RETURNED=0`, Primer3 tried but failed to design any primer meeting the constraints.

**What to try:**
- Relax the Primer3 constraints in `config/primer3_settings.txt` (e.g. widen the Tm range: `PRIMER_MIN_TM=57`, `PRIMER_MAX_TM=63`)
- Increase `--min_trimmed_length` is *not* a solution here — if anything, the design window may be too short. Lower `--min_trimmed_length` to give Primer3 a longer template.
- In distance mode, increase `--template_length` to give Primer3 a larger region to work with.

---

### How to check if primers passed or failed the alignment filtering

**File to open:** `primer_alignment_report.tsv`

This file contains every primer Primer3 designed, one per row. Look at `num_alignments`:

- **`0` alignments** — the primer did not map to any transcript. It will not appear in `primer_alignment_summary.tsv` and will be excluded from all downstream steps. This can happen if:
  - The peak/design window overlapped intronic sequence (even after trimming, the primer may sit near an exon boundary)
  - There is a version mismatch between the genome FASTA used for design and the transcriptome Bowtie2 index
  - The primer sequence is genuinely absent from all transcripts in the index

- **`1` alignment, to the correct gene** — ideal. This primer is specific.

- **Multiple alignments to the same gene** — fine; the primer covers multiple isoforms of your target.

- **Alignments to other genes** — problematic. The primer may not be specific. These are flagged during the optimisation step.

**File to open:** `best_primers.tsv`

Genes present here have at least one primer that:
1. Aligned with ≤ `max_mismatches` mismatches
2. Aligned within `distance_threshold` bases of the 3′ end
3. Shows specificity for the correct gene

If a gene is in `primer_alignment_report.tsv` but **absent** from `best_primers.tsv`, all its primers were filtered out because they either:
- Had too many mismatches
- Mapped too far from the 3′ end
- Cross-reacted with other genes

**What to try:**
- Increase `--max_mismatches` (e.g. from `0` to `1` or `2`) to allow slightly imperfect matches
- Increase `--distance_threshold` (e.g. to `2000`) if your target region is farther from the 3′ end
- Check `primer_alignment_summary.tsv` for that gene to understand exactly where its primers are mapping

---

### Quick summary of the diagnostic flow

```
Gene of interest
    │
    ▼
peaks_qc_summary.tsv ──► final_selection = FALSE? → check failure_reason
    │
    ▼ (peak found)
cdna_primers.tsv ──────► gene missing? → Primer3 failed → relax primer3_settings.txt
    │
    ▼ (primers designed)
primer_alignment_report.tsv ──► num_alignments = 0? → primer didn't map to transcriptome
    │
    ▼ (primers mapped)
best_primers.tsv ──────► gene missing? → failed specificity/distance filter
    │
    ▼ (passed filters)
best_primers_optimal.tsv ──► final primer recommendation
```

---

## Appendix — Running the summary report manually

The pipeline includes a script that produces a concise text summary of results. Run it after the pipeline completes:

```bash
python bin/summarize_pipeline_results.py results/
```

This prints how many genes got primers, how many failed at each step, and lists the failure reasons.
