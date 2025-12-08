# Distance-Based Primer Design Workflow

## Overview

The distance-based primer design workflow is an alternative to the peak-based primer design. Instead of using RNA-seq coverage peaks to select primer targets, this workflow designs primers at a fixed distance from the 3' end of transcripts.

This is particularly useful for:
- Designing primers for 3' RNA-seq protocols (e.g., QuantSeq, 3'-Tag-seq)
- Standardizing primer positions across different samples
- Targeting specific regions without requiring BAM files

## Workflow Modes

The workflow supports two input modes:

### Mode 1: Gene ID List (Automatic MANE Transcript Fetching)
Provide a list of Ensembl gene IDs, and the pipeline will:
1. Fetch MANE Select transcript IDs and sequences from Ensembl REST API
2. Extract N bases from the 3' end
3. Design primers using Primer3

### Mode 2: Transcript FASTA File
Provide a FASTA file with transcript sequences, and the pipeline will:
1. Extract N bases from the 3' end of each sequence
2. Design primers using Primer3

## Usage

### Basic Usage with Gene IDs

```bash
nextflow run main.nf \
  --distance_mode \
  --genes gene_ids.txt \
  --template_length 300 \
  --outdir results_distance
```

### Usage with Transcript FASTA

```bash
nextflow run main.nf \
  --distance_mode \
  --transcript_fasta transcripts.fasta \
  --template_length 300 \
  --outdir results_distance
```

### With Transcriptome Alignment (Specificity Check)

```bash
nextflow run main.nf \
  --distance_mode \
  --genes gene_ids.txt \
  --template_length 300 \
  --transcriptome_index /path/to/transcriptome_idx \
  --outdir results_distance
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--distance_mode` | Enable distance-based primer design mode |
| `--template_length` | Number of bases to extract from 3' end (e.g., 300) |
| `--genes` OR `--transcript_fasta` | Input: gene ID list or transcript FASTA |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--transcriptome_index` | Bowtie2 index for primer specificity checking | null |
| `--primer3_settings` | Primer3 configuration file | `config/primer3_settings.txt` |
| `--outdir` | Output directory | `results` |

## Input File Formats

### Gene ID List (`--genes`)

Plain text file with one Ensembl gene ID per line:

```
ENSG00000139618
ENSG00000157764
ENSG00000141510
```

### Transcript FASTA (`--transcript_fasta`)

Standard FASTA format:

```
>ENST00000380152 BRCA2
ATGGCGATTCCCGAAGATGGTGATGCGGATCGGG...
>ENST00000269305 TP53
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAG...
```

## Output Files

The workflow generates the following outputs in the `--outdir`:

### Mode 1 (Gene IDs) Output:
- `mane_transcripts.fasta` - MANE Select transcript sequences from Ensembl
- `mane_mapping.tsv` - Mapping of gene IDs to transcript IDs
- `3prime_sequences.fasta` - Extracted 3' end sequences (templates)
- `primer3_input.txt` - Primer3 input file
- `primer3_output.txt` - Designed primers from Primer3
- `primers.fasta` - Primer sequences in FASTA format

### Mode 2 (Transcript FASTA) Output:
- `3prime_sequences.fasta` - Extracted 3' end sequences (templates)
- `primer3_input.txt` - Primer3 input file
- `primer3_output.txt` - Designed primers from Primer3
- `primers.fasta` - Primer sequences in FASTA format

### Additional Output (if `--transcriptome_index` provided):
- `primers_alignment.bam` - Primer alignments to transcriptome
- `primers_alignment.bam.bai` - BAM index
- `alignment_stats.txt` - Alignment statistics

## Workflow Steps

### Step 1: Input Processing

**Option 1 (Gene IDs):**
```
Gene IDs → Ensembl REST API → MANE transcripts
```

**Option 2 (Transcript FASTA):**
```
Transcript FASTA → Direct use
```

### Step 2: 3' End Extraction

Extract N bases from the 3' end of each transcript:

```python
# If sequence length > template_length:
template = sequence[-template_length:]
# Otherwise:
template = sequence  # Use entire sequence
```

### Step 3: Primer Design

Run Primer3 with the extracted templates using the same configuration as the peak-based workflow.

### Step 4: Specificity Check (Optional)

If `--transcriptome_index` is provided, align primers to the transcriptome using Bowtie2 to check for off-target binding.

## Example Workflows

### Example 1: Basic Gene-Based Design

```bash
# Create gene list
echo -e "ENSG00000139618\nENSG00000157764\nENSG00000141510" > genes.txt

# Run pipeline
nextflow run main.nf \
  --distance_mode \
  --genes genes.txt \
  --template_length 300 \
  --outdir results_distance_basic
```

### Example 2: Custom Transcripts with QC

```bash
# Run with custom transcripts and alignment QC
nextflow run main.nf \
  --distance_mode \
  --transcript_fasta my_transcripts.fasta \
  --template_length 400 \
  --transcriptome_index /data/indices/human_transcriptome \
  --outdir results_distance_qc
```


## Comparison: Distance-Based vs Peak-Based

| Feature | Distance-Based | Peak-Based |
|---------|---------------|------------|
| Input Required | Gene IDs or Transcript FASTA | BAM + GTF + Gene IDs |
| Primer Location | Fixed distance from 3' end | Coverage peak regions |
| Coverage Analysis | Not required | Required (MACS2) |
| Use Case | Standardized 3' targeting | Data-driven peak targeting |
| Speed | Fast (no BAM processing) | Slower (peak calling) |

## Technical Details

### MANE Select Transcripts

MANE (Matched Annotation from NCBI and EMBL-EBI) Select represents a collaboration between NCBI and Ensembl to identify a default transcript per protein-coding gene that is:
- Well-supported clinically
- Identical between RefSeq and Ensembl/GENCODE
- Suitable as a reference for variant interpretation

### Ensembl REST API

The workflow uses the Ensembl REST API with rate limiting:
- Base URL: `https://rest.ensembl.org`
- Rate limit: ~0.15 seconds between requests
- Endpoints used:
  - `/lookup/id/{gene_id}?expand=1` - Get gene/transcript info
  - `/sequence/id/{transcript_id}?type=cdna` - Get cDNA sequence

### Primer3 Configuration

The workflow reuses the same Primer3 settings from the peak-based workflow (`config/primer3_settings.txt`), ensuring consistency across both modes.

## Troubleshooting

### Issue: "No MANE Select transcript found"

Some genes may not have MANE Select transcripts. Consider:
- Using canonical transcripts instead
- Providing a custom transcript FASTA with Mode 2

### Issue: "Template length exceeds transcript length"

If your `--template_length` is larger than some transcripts:
- The workflow will use the entire transcript sequence
- Check output headers for `original_length` vs `3prime_extracted` values

### Issue: API rate limiting

If you have many genes (>100):
- The workflow includes built-in rate limiting
- Consider running overnight for large gene lists
- Alternative: Use Mode 2 with pre-downloaded transcripts

## Integration with Main Pipeline

The distance-based workflow is fully integrated into the main pipeline:

```groovy
workflow {
  if (params.distance_mode) {
    distance_primer_design()  // This workflow
  } else if (params.makeplots && !params.bam) {
    makeplots()
  } else {
    primer_design()  // Original peak-based workflow
  }
}
```

## Future Enhancements

Potential improvements for future versions:
- Support for canonical transcript fallback when MANE not available
- Batch API requests for improved performance
- Integration with primer alignment analysis
- Support for selecting best primers based on alignment results
- Custom transcript selection rules (e.g., longest, most expressed)

## References

- Ensembl REST API: https://rest.ensembl.org/
- MANE Project: https://www.ncbi.nlm.nih.gov/refseq/MANE/
- Primer3: https://primer3.org/
