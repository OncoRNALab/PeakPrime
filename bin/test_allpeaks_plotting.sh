#!/bin/bash
# Test script for new all-peaks plotting feature

# This script demonstrates the new plotting feature that shows all MACS2 peaks
# within the gene span as horizontal bars in the coverage plots

echo "Testing new all-peaks plotting feature..."

# Example usage (adjust paths as needed):
# Rscript bin/MakePlots_new.R \
#   --gene ENSG00000013306 \
#   --bw path/to/sample.bam.bw \
#   --gtf path/to/annotations.gtf \
#   --peaks path/to/peaks.tsv \
#   --qc path/to/qc_summary.tsv \
#   --primer path/to/primer_targets.bed \
#   --narrowpeak path/to/sample_peaks.narrowPeak \
#   --out test_plot_with_all_peaks.png

echo "New feature capabilities:"
echo "1. Reads MACS2 narrowPeak file to show ALL peaks in gene"
echo "2. Displays peaks as horizontal bars colored by intensity"
echo "3. Adds new 'All Peaks' track at bottom of plot"
echo "4. Maintains backward compatibility (optional --narrowpeak parameter)"
echo "5. Legend shows peak score intensity when peaks are present"

echo ""
echo "To use with Nextflow makeplots workflow:"
echo "nextflow run main.nf --makeplots \\"
echo "  --bw ./results/sample.bam.bw \\"
echo "  --gtf annotations.gtf \\"
echo "  --genes target_genes.txt \\"
echo "  --peaks_tsv ./results/peaks.tsv \\"
echo "  --qc_tsv ./results/qc_summary.tsv \\"
echo "  --primer ./results/primer_targets.bed \\"
echo "  --narrowpeak ./results/macs2_peaks/sample_peaks.narrowPeak"
