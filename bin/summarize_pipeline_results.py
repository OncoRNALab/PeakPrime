#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd

def find_file(pattern, outdir):
    files = glob.glob(os.path.join(outdir, '**', pattern), recursive=True)
    return files[0] if files else None

def main(outdir):
    # First, get all target genes from QC summary
    qc_file = find_file('peaks_qc_summary.tsv', outdir)
    all_target_genes = set()
    genes_without_peak = set()
    genes_with_peak = set()
    
    if qc_file and os.path.getsize(qc_file) > 0:
        qc_df = pd.read_csv(qc_file, sep='\t')
        all_target_genes = set(str(gid) for gid in qc_df['gene_id'])
        # Genes with selected peaks (final_selection == True)
        genes_with_peak = set(str(gid) for gid in qc_df[qc_df['final_selection'] == True]['gene_id'])
        # Genes without selected peaks (final_selection == False)
        genes_without_peak = set(str(gid) for gid in qc_df[qc_df['final_selection'] == False]['gene_id'])
    else:
        # Fallback: try to get genes with peaks from selected_peaks.tsv
        peaks_file = find_file('selected_peaks.tsv', outdir)
        if peaks_file and os.path.getsize(peaks_file) > 0:
            peaks_df = pd.read_csv(peaks_file, sep='\t')
            genes_with_peak = set(str(gid) for gid in peaks_df['gene'])
    
    print(f"Total target genes: {len(all_target_genes)}")
    print(f"Genes with selected peak: {len(genes_with_peak)}")
    print(f"Genes WITHOUT selected peak: {len(genes_without_peak)}")

    # 2. Genes with at least one Primer3-generated primer
    primers_file = find_file('cdna_primers.tsv', outdir) or find_file('primer3_output.txt', outdir)
    genes_with_primer = set()
    if primers_file and os.path.getsize(primers_file) > 0:
        if primers_file.endswith('.tsv'):
            primers_df = pd.read_csv(primers_file, sep='\t')
            if 'gene_id' in primers_df.columns:
                genes_with_primer = set(str(gid) for gid in primers_df['gene_id'])
        else:
            # Fallback: parse Primer3 output for SEQUENCE_ID
            with open(primers_file) as f:
                for line in f:
                    if line.startswith('SEQUENCE_ID='):
                        genes_with_primer.add(str(line.strip().split('=',1)[1]))
    print(f"Genes with at least one Primer3 primer: {len(genes_with_primer)}")

    # 3. Genes with at least one aligned primer
    align_file = find_file('primer_alignment_summary.tsv', outdir)
    genes_with_aligned = set()
    if align_file and os.path.getsize(align_file) > 0:
        align_df = pd.read_csv(align_file, sep='\t')
        if 'gene_id' in align_df.columns:
            genes_with_aligned = set(str(gid) for gid in align_df['gene_id'])
    print(f"Genes with at least one aligned primer: {len(genes_with_aligned)}")

    # 4. Genes with at least one uniquely aligned or best primer
    # Use best_primers.tsv which contains primers that passed selection criteria
    best_file = find_file('best_primers.tsv', outdir)
    genes_with_unique = set()
    if best_file and os.path.getsize(best_file) > 0:
        best_df = pd.read_csv(best_file, sep='\t')
        if 'gene_id' in best_df.columns:
            genes_with_unique = set(str(gid) for gid in best_df['gene_id'])
    print(f"Genes with at least one uniquely aligned or best primer: {len(genes_with_unique)}")

    print("\nSummary Table:")
    print(f"{'Category':45} | {'Gene count'}")
    print(f"{'-'*45}-+{'-'*10}")
    print(f"Total target genes{'':24} | {len(all_target_genes)}")
    print(f"Genes with selected peak{'':18} | {len(genes_with_peak)}")
    print(f"Genes WITHOUT selected peak{'':16} | {len(genes_without_peak)}")
    print(f"Genes with at least one Primer3 primer{'':4} | {len(genes_with_primer)}")
    print(f"Genes with at least one aligned primer{'':6} | {len(genes_with_aligned)}")
    print(f"Genes with at least one uniquely aligned/best primer | {len(genes_with_unique)}")

    # Output gene IDs with selected peak but without a uniquely aligned/best primer
    missing_best_primer = sorted(genes_with_peak - genes_with_unique)
    print(f"\nGenes with selected peak but WITHOUT a uniquely aligned/best primer: {len(missing_best_primer)}")
    if missing_best_primer:
        print("Gene IDs:")
        for gid in missing_best_primer:
            print(gid)

    # Output gene IDs without selected peaks
    genes_without_peak_sorted = sorted(genes_without_peak)
    print(f"\nGenes WITHOUT selected peak: {len(genes_without_peak_sorted)}")
    if genes_without_peak_sorted:
        print("Gene IDs:")
        for gid in genes_without_peak_sorted:
            print(gid)
        
        # If QC data is available, show failure reasons
        if qc_file and os.path.getsize(qc_file) > 0:
            failed_df = qc_df[qc_df['final_selection'] == False]
            if not failed_df.empty:
                print("\nFailure reasons summary:")
                failure_counts = failed_df['failure_reason'].value_counts()
                for reason, count in failure_counts.items():
                    print(f"  {reason}: {count} genes")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <pipeline_output_dir>")
        sys.exit(1)
    main(sys.argv[1])
