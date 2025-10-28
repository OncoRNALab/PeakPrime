# Transcript-to-Gene Mapping

## Overview

To improve pipeline performance, we use a **pre-built transcript-to-gene mapping file** instead of parsing the large transcriptome FASTA file during every run.

## Default Mapping File

**Location**: `resources/transcript_gene_mapping.GRCh38.109.tsv`

This file maps Ensembl transcript IDs to gene names for the GRCh38.109 transcriptome and is automatically used by the pipeline.

## How It's Used

When primers are aligned to the transcriptome, the results show transcript IDs (e.g., `ENST00000456328`). The mapping file converts these to readable gene names (e.g., `DDX11L2`) in the QC reports.

### Benefits:
- **Much faster**: No need to parse large FASTA files during each run
- **Cached**: Mapping is done once and reused for all runs
- **Automatic**: Works out-of-the-box with default parameters

## Generating the Mapping File

### For the Default Transcriptome

The mapping file for GRCh38.109 can be generated once:

```bash
# Submit as a job (recommended)
qsub scripts/generate_transcript_mapping.pbs

# Or run interactively (if R module is loaded)
ml R/4.3.2-gfbf-2023a
bash scripts/generate_transcript_mapping.sh
```

This creates: `resources/transcript_gene_mapping.GRCh38.109.tsv`

### For a Different Transcriptome

If you're using a different transcriptome, generate a custom mapping:

```bash
Rscript bin/create_transcript_gene_mapping.R \
  --transcriptome_fasta /path/to/your/transcriptome.fa \
  --output /path/to/your/mapping.tsv
```

Then specify it in your workflow:

```bash
nextflow run main.nf \
  --transcriptome_index /path/to/your/index \
  --transcript_mapping /path/to/your/mapping.tsv \
  ...
```

## File Format

The mapping file is a simple TSV with two columns:

```
transcript_id	gene_name
ENST00000456328	DDX11L2
ENST00000450305	DDX11L2
ENST00000488147	WASH7P
...
```

## Fallback Behavior

If the mapping file doesn't exist, the pipeline will:
1. Try to parse `--transcriptome_fasta` if provided (slow)
2. Otherwise, show gene names as "Unknown" in reports

## Parameters

- `--transcript_mapping` (default: `resources/transcript_gene_mapping.GRCh38.109.tsv`)
  - Path to pre-built mapping file
  - Much faster than parsing FASTA
  
- `--transcriptome_fasta` (optional)
  - Only needed if mapping file doesn't exist
  - Slower fallback for gene name extraction
