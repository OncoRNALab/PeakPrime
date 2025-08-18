#!/bin/bash

# Wrapper script to create transcript mapping if it doesn't exist
# and then run the analysis with the optimized approach

set -euo pipefail

# Input parameters
TRANSCRIPTOME_FASTA="$1"
OUTPUT_MAPPING="$2"
ALIGNMENT_BAM="$3"
PRIMERS_TSV="$4"
OUT_REPORT="$5"
OUT_SUMMARY="$6"

# Check if mapping file already exists
if [ ! -f "$OUTPUT_MAPPING" ]; then
    echo "Creating transcript-to-gene mapping file..."
    Rscript create_transcript_mapping.R \
        --transcriptome_fasta "$TRANSCRIPTOME_FASTA" \
        --output_mapping "$OUTPUT_MAPPING"
else
    echo "Using existing transcript mapping: $OUTPUT_MAPPING"
fi

# Run the analysis with the mapping file
echo "Running primer alignment analysis with pre-built mapping..."
Rscript analyze_primer_alignments.R \
    --alignment_bam "$ALIGNMENT_BAM" \
    --primers_tsv "$PRIMERS_TSV" \
    --transcript_mapping "$OUTPUT_MAPPING" \
    --out_report "$OUT_REPORT" \
    --out_summary "$OUT_SUMMARY"

echo "Analysis complete!"
