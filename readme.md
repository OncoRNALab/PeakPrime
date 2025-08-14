# RNA‑seq Primer Peak Pipeline (Isoform‑agnostic)

This Nextflow DSL2 pipeline picks a robust, isoform‑agnostic high‑coverage coordinate per gene from RNA‑seq data, extracts a ±60 bp window, and designs primers with Primer3.

## TL;DR

```bash
# Basic usage with FASTA file (recommended)
nextflow run main.nf -prof**Q: Can I use a genome source that differs from the BAM/GTF build?**

No. Always match builds (e.g., hg38 with hg38). Whether using FASTA or BSgenome, the genome build must be consistent with your BAM and GTF files. Mismatches will yield incorrect sequences.

**Q: Should I use FASTA or BSgenome for sequence extraction?**

FASTA files are generally recommended as they're more flexible and faster for large-scale analyses. BSgenome packages are convenient for standard genomes but require installation of specific Bioconductor packages.

**Q: Can I run multiple BAMs?**

This version expects one BAM at a time. A multi‑sample mode (e.g., pick peaks on a merged coverage or per sample) can be added.

**Q: How do I choose QC thresholds?**

- Start without thresholds to see your data's coverage distribution in `qc_coverage_summary.tsv`
- For `min_window_mean_pct`: 25% is a good starting point, 50% is conservative, 10% is permissive
- For `max_gap`: 10bp is reasonable for most primers, 5bp is strict

**Q: When should I use the sliding window feature?**

Enable `--sliding_window true` when you want the pipeline to automatically find better primer locations if the initial peak-centered window has low coverage or large gaps.-bam sample.sorted.bam \
  --gtf genes.gtf \
  --genes genes.txt \
  --fasta genome.fa \
  --outdir results

# With advanced QC (recommended for production)
nextflow run main.nf -profile conda \
  --bam sample.sorted.bam \
  --gtf genes.gtf \
  --genes genes.txt \
  --fasta genome.fa \
  --sliding_window true \
  --min_window_mean_pct 25 \
  --trim_to_exon true \
  --outdir results
```

---

## Contents

- [Overview](#overview)
- [Key features](#key-features)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Parameters](#parameters)
- [Primer3 configuration](#primer3-configuration)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Cite](#cite)

---

## Overview

The pipeline computes per‑base coverage over the **union of exons per gene** (all isoforms), smooths the signal, selects the highest‑coverage coordinate, extracts a 120‑nt template (±60 bp), and feeds it to Primer3 to propose PCR primers.

**Why isoform‑agnostic?** Any base that is exonic in *any* isoform is eligible, so you land primers in the most supported portion of the gene in your data — without committing to a specific transcript.

> If you need isoform‑aware logic (dominant transcript or constitutive exons only), open an issue — the pipeline can be extended with optional modes.

---

## Key features

- ⚡ **Fast**: converts BAM → BigWig via `megadepth` and queries coverage efficiently.
- 🧬 **RNA‑aware**: restricts to **exons** only (union across isoforms).
- 🧽 **Noise‑robust**: windowed smoothing before peak picking.
- 🧪 **Primer‑ready**: emits FASTA targets and runs Primer3.
- 📊 **QC**: per‑gene coverage summaries.
- 🔧 **Robust**: handles complex genomes with spike-ins and non-standard chromosomes automatically using direct BigWig reading.

---

## Inputs

- `--bam` (**required**) – Spliced RNA‑seq BAM, sorted (`.bam`). Index will be created automatically if missing.
- `--gtf` (**required**) – Gene annotation in GTF for the **same genome build** as the BAM.
- `--genes` (**required**) – Text file with one Ensembl gene ID per line (see example below).
- `--fasta` (**optional**) – Genome FASTA file for sequence extraction (preferred method).
- `--genome_package` (**optional**) – Bioconductor BSgenome package name (fallback if no FASTA provided).

> **Sequence source priority:** If both `--fasta` and `--genome_package` are provided, FASTA takes priority. If neither is provided, defaults to `BSgenome.Hsapiens.NCBI.GRCh38`.

> **Build consistency is critical:** BAM ↔ GTF ↔ genome source must match (e.g., all hg38). Mismatches will yield wrong sequences.

**Example gene list file (`genes.txt`):**
```
ENSG00000141510
ENSG00000155657
ENSG00000157764
```
*Note: Use Ensembl gene IDs that match your GTF annotation. One ID per line, no headers or extra formatting.*

---

## Outputs (in `--outdir`)

- `*.bw` – BigWig coverage track.
- `primer_targets.fa` – FASTA with one 120‑nt sequence per gene (`>GENE|chr:start-end(strand)`).
- `primer_targets.bed` – BED coordinates of each target window.
- `peaks.tsv` – Chosen window per gene (gene, chr, start, end, strand).
- `qc_coverage_summary.tsv` – Per‑gene coverage summary (total exonic bases, max/mean/median).
- `primer3_input.txt` / `primer3_output.txt` – Primer3 request/response.

---

## Installation

### Requirements

- **Nextflow** ≥ 22
- **Java** ≥ 11
- **Conda/Mamba** (recommended) or your own environment manager

### Get the code

```
# from the repo root
ls
# main.nf, nextflow.config, envs/, bin/, config/
```

### Environments

The pipeline ships with Conda envs for all steps. Use `-profile conda` to auto‑create them.

---

## Quick start

1. **Prepare inputs**

   - Sort your BAM if needed:
     ```bash
     samtools sort -o sample.sorted.bam sample.bam
     ```
     (The pipeline will create the index automatically)
   - Ensure your GTF matches the BAM genome build.
   - Create `genes.txt` with Ensembl IDs (one per line).

2. **Run**

   ```bash
   # Using FASTA file (recommended)
   nextflow run main.nf -profile conda \
     --bam sample.sorted.bam \
     --gtf genes.gtf \
     --genes genes.txt \
     --fasta genome.fa \
     --pad 60 \
     --smooth_k 31 \
     --primer3_settings config/primer3_settings.txt \
     --outdir results
   ```

   ```bash
   # Using BSgenome package (alternative)
   nextflow run main.nf -profile conda \
     --bam sample.sorted.bam \
     --gtf genes.gtf \
     --genes genes.txt \
     --genome_package BSgenome.Hsapiens.UCSC.hg38 \
     --pad 60 \
     --smooth_k 31 \
     --primer3_settings config/primer3_settings.txt \
     --outdir results
   ```

3. **Inspect results**

   - `results/primer3_output.txt` for designed primers
   - `results/peaks.tsv` & `results/primer_targets.bed` to see chosen windows

---

## Parameters

### Core Parameters

| Name                 | Type   | Default                          | Description                                                                    |
| -------------------- | ------ | -------------------------------- | ------------------------------------------------------------------------------ |
| `--bam`              | path   | *none*                           | Input RNA‑seq BAM (sorted). Index created automatically if missing.           |
| `--gtf`              | path   | *none*                           | Gene annotation GTF (same build).                                              |
| `--genes`            | path   | *none*                           | File with Ensembl gene IDs (one per line).                                     |
| `--fasta`            | path   | *none*                           | Genome FASTA file for sequence extraction (preferred method).                  |
| `--genome_package`   | string | `BSgenome.Hsapiens.NCBI.GRCh38`  | Bioconductor genome used for sequence extraction (fallback if no FASTA). Must match build of BAM/GTF. |
| `--pad`              | int    | 60                               | Half‑window size around peak (±bp). Total sequence length = `2*pad + 1`.       |
| `--smooth_k`         | int    | 31                               | Window size for running mean smoothing (odd integer).                          |
| `--primer3_settings` | path   | `config/primer3_settings.txt`    | Primer3 constraints.                                                           |
| `--outdir`           | path   | `results`                        | Output directory.                                                              |

### Window Coverage QC Parameters (Advanced)

| Name                      | Type   | Default | Description                                                                    |
| ------------------------- | ------ | ------- | ------------------------------------------------------------------------------ |
| `--sliding_window`        | bool   | false   | Enable sliding-window selection if initial peak window fails QC thresholds.   |
| `--min_window_mean`       | float  | null    | Minimum absolute mean coverage across final window. Leave null to disable.    |
| `--min_window_mean_pct`   | float  | null    | **Dynamic threshold** - minimum window mean as % (0-100) of gene's peak coverage. Overrides `min_window_mean`. |
| `--max_gap`               | int    | null    | Maximum allowed longest zero-coverage run within final window.                 |
| `--search_slop`           | int    | 1000    | Extra bases around exon span for BigWig import cache.                         |

### Exon Overlap QC Parameters (Advanced)

| Name                      | Type   | Default | Description                                                                    |
| ------------------------- | ------ | ------- | ------------------------------------------------------------------------------ |
| `--trim_to_exon`          | bool   | false   | Trim final window to boundaries of exon containing the peak.                   |
| `--min_exonic_fraction`   | float  | null    | Minimum required fraction (0-1) of window overlapping exons. Flags failures but still reports. |

### Usage Examples

**Basic usage:**
```bash
nextflow run main.nf -profile conda \
  --bam sample.bam --gtf genes.gtf --genes gene_list.txt --fasta genome.fa
```

**With dynamic QC (recommended for production):**
```bash
nextflow run main.nf -profile conda \
  --bam sample.bam --gtf genes.gtf --genes gene_list.txt --fasta genome.fa \
  --sliding_window true --min_window_mean_pct 25 --trim_to_exon true
```

**Conservative quality filtering:**
```bash
nextflow run main.nf -profile conda \
  --bam sample.bam --gtf genes.gtf --genes gene_list.txt --fasta genome.fa \
  --sliding_window true --min_window_mean_pct 50 --max_gap 5 --min_exonic_fraction 0.9
```

---

## Primer3 configuration

Edit `config/primer3_settings.txt` to fit your assay. Example values are provided:

```
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=35.0
PRIMER_MAX_GC=65.0
PRIMER_PRODUCT_SIZE_RANGE=120-180
PRIMER_NUM_RETURN=5
PRIMER_EXPLAIN_FLAG=1
```

> You can restrict to amplicon lengths that fit your downstream assay; Primer3 will design inward‑facing primers within the 120‑nt template.

---

## Reproducibility

- Use the `conda` profile to pin tool versions (`megadepth`, `primer3`, R/Bioconductor packages).
- Record the commit hash of this repo and your `nextflow` version.
- Archive `nextflow.config` and any changes to `envs/*.yml`.

---

## Troubleshooting

**No peaks found**

- Check that gene IDs in `genes.txt` exist in the GTF.
- Confirm BAM, GTF, and genome source (FASTA or `--genome_package`) use the **same genome build**.
- Coverage may be too low — try a larger smoothing window (`--smooth_k 51`) or verify expression.

**Sequence extraction errors**

- Ensure your FASTA file is indexed (`.fai` file should exist alongside the FASTA).
- If using BSgenome, install the correct package in the conda environment or change `--genome_package` to a package you have.
- Verify chromosome naming consistency between BAM/GTF and genome source (e.g., "chr1" vs "1").

**BigWig coverage processing with spike-ins**

- The pipeline uses direct BigWig reading via `rtracklayer` for robust coverage extraction.
- If your genome/GTF contains spike-in sequences (ERCC controls, custom chromosomes), they are automatically filtered out.
- Only genes on standard chromosomes (1-22, X, Y, MT) will be processed.
- This approach handles complex genomes more reliably than previous methods.

**No coverage detected**

- If the pipeline reports "max coverage: 0" for all genes, check that:
  - Your BAM file actually contains reads aligned to the target chromosomes
  - The BigWig conversion was successful (check for non-zero file size)
  - Chromosome names match between BAM, GTF, and reference genome (especially "chr" prefixes)
  - Your genes are actually expressed in the sample (try different gene IDs as a test)

**Primer3 returns no primers**

- Loosen constraints in `config/primer3_settings.txt` (e.g., GC%, Tm, product size range).

**Conda env creation is slow**

- Prefer `mamba` if available: `NXF_CONDA_CLI=mamba nextflow run ...`

**Cluster/scheduler usage**

- Add a profile in `nextflow.config` similar to:
  ```groovy
  profiles {
    slurm {
      process.executor = 'slurm'
      process.queue = 'your-queue'
      process.cpus = 4
      process.memory = '8 GB'
      process.time = '2 h'
      process.conda = true
    }
  }
  ```

---

## FAQ

**Q: Can I use a Bioconductor genome that differs from the BAM/GTf build?**

> No. Always match builds (e.g., hg38 with hg38). Otherwise you’ll extract incorrect sequences.

**Q: Can I run multiple BAMs?**

> This version expects one BAM at a time. A multi‑sample mode (e.g., pick peaks on a merged coverage or per sample) can be added.

**Q: How do I avoid isoform‑specific exons?**

> Use a constitutive‑exon mode (not included here yet) or design per‑transcript.

**Q: Can I mask SNPs?**

> Yes — extend `pick_peaks.R` to soft‑mask bases overlapping a SNP BED (replace with `N`) before writing FASTA.

---

## Cite

If this pipeline helps your work, please cite the tools it wraps:

- **rtracklayer** (Bioconductor) - for BigWig coverage extraction
- **megadepth** - for BAM to BigWig conversion  
- **Primer3** - for primer design
- **GenomicRanges/GenomicFeatures** (Bioconductor) - for genomic interval operations

> Add exact references/versions from your environment lockfiles for manuscripts.

