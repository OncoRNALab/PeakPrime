#!/usr/bin/env python3
"""
Select best primers from alignment reports.

3-stage filtering strategy for 3' RNA-seq protocols:
1. Require mismatches == 0 (perfect alignment quality)
2. Require distance_to_end <= threshold (protocol-specific relevance)
3. Require unique gene mapping (specificity, only for relevant hits)

This allows primers with off-target hits to other genes as long as 
those hits are far from the 3' end (where amplification occurs).

Writes out 'best_primers.tsv' in the current working directory.
"""
import argparse
import pandas as pd
import sys


def main():
    p = argparse.ArgumentParser(description='Select best primers with 3-stage filtering')
    p.add_argument('--report', required=True)
    p.add_argument('--summary', required=True)
    p.add_argument('--out', default='best_primers.tsv')
    p.add_argument('--distance_threshold', type=float, default=1000.0,
                   help='Maximum distance from 3\' end for alignments to be considered relevant (bp). '
                        'Default: 1000. Use larger values (e.g., 2000-5000) for more permissive filtering.')
    args = p.parse_args()

    report = pd.read_csv(args.report, sep='\t', dtype=str)
    try:
        summary = pd.read_csv(args.summary, sep='\t', dtype=str)
    except Exception:
        # fallback: use quick heuristic from report (unique alignments only)
        out = report[report.get('num_alignments', '0').astype(int) == 1].copy()
        # Remove alignment_quality from fallback output if present
        if 'alignment_quality' in out.columns:
            out = out.drop('alignment_quality', axis=1)
        out.to_csv(args.out, sep='\t', index=False)
        print(f'Wrote {args.out} (fallback quick filter)')
        return

    # Normalize columns
    if 'primer_id' not in summary.columns:
        # Construct primer_id from components if not present
        if all(col in summary.columns for col in ['gene_id', 'primer_index', 'primer_type', 'gene_strand']):
            summary['primer_id'] = (summary['gene_id'].astype(str) + "|idx" + 
                                  summary['primer_index'].astype(str) + "|" + 
                                  summary['primer_type'].astype(str) + "|strand" + 
                                  summary['gene_strand'].astype(str))
        else:
            print('detailed summary missing required columns for primer_id construction; falling back to report heuristic', file=sys.stderr)
            out = report[report.get('num_alignments', '0').astype(int) == 1].copy()
            # Remove alignment_quality from fallback output
            if 'alignment_quality' in out.columns:
                out = out.drop('alignment_quality', axis=1)
            out.to_csv(args.out, sep='\t', index=False)
            print(f'Wrote {args.out} (fallback quick filter)')
            return

    # ===== STAGE 1: Filter for perfect matches (zero mismatches) =====
    print("\n=== 3-STAGE FILTERING FOR 3' RNA-SEQ ===")
    print(f"Total alignments in summary: {len(summary)}")
    
    summary['mismatches'] = pd.to_numeric(summary.get('mismatches', 9999), errors='coerce').fillna(9999).astype(int)
    stage1_filtered = summary[summary['mismatches'] == 0].copy()
    
    print(f"STAGE 1 (mismatches=0): {len(stage1_filtered)} alignments retained ({len(stage1_filtered)/max(len(summary),1)*100:.1f}%)")
    
    if stage1_filtered.empty:
        print("No perfect-match alignments found. Creating empty output.")
        cols = list(report.columns) + ['aligned_gene_name', 'zero_mismatch_alignments', 'distance_to_end_min']
        if 'alignment_quality' in cols:
            cols.remove('alignment_quality')
        if 'best_mapq' in cols:
            cols.remove('best_mapq')
        pd.DataFrame(columns=cols).to_csv(args.out, sep='\t', index=False)
        print(f'Wrote {args.out} (no perfect matches)')
        return
    
    # ===== STAGE 2: Filter by distance to 3' end =====
    # Convert distance_to_end to numeric, treating "NA" as a very large value (far from end)
    stage1_filtered['distance_to_end_numeric'] = pd.to_numeric(
        stage1_filtered.get('distance_to_end', 'NA'), 
        errors='coerce'
    ).fillna(float('inf'))  # Treat NA as infinitely far (will be filtered out)
    
    stage2_filtered = stage1_filtered[stage1_filtered['distance_to_end_numeric'] <= args.distance_threshold].copy()
    
    print(f"STAGE 2 (distance≤{args.distance_threshold}bp from 3' end): {len(stage2_filtered)} alignments retained ({len(stage2_filtered)/len(stage1_filtered)*100:.1f}%)")
    print(f"  Filtered out: {len(stage1_filtered) - len(stage2_filtered)} alignments too far from 3' end")
    
    if stage2_filtered.empty:
        print(f"No alignments within {args.distance_threshold}bp of 3' end. Creating empty output.")
        cols = list(report.columns) + ['aligned_gene_name', 'zero_mismatch_alignments', 'distance_to_end_min']
        if 'alignment_quality' in cols:
            cols.remove('alignment_quality')
        if 'best_mapq' in cols:
            cols.remove('best_mapq')
        pd.DataFrame(columns=cols).to_csv(args.out, sep='\t', index=False)
        print(f'Wrote {args.out} (no alignments near 3\' end)')
        return
    
    # ===== STAGE 3: Require unique gene mapping (only for relevant alignments) =====

    # ===== STAGE 3: Require unique gene mapping (only for relevant alignments) =====
    selected = []
    primers_with_multi_gene = []
    
    for pid, group in stage2_filtered.groupby('primer_id'):
        genes = group['aligned_gene_name'].dropna().unique()
        # require unique gene (only considering alignments near 3' end)
        if len(genes) == 1:
            rpt = report[report['primer_id'] == pid]
            if not rpt.empty:
                row = rpt.iloc[0].to_dict()
                if 'alignment_quality' in row:
                    del row['alignment_quality']
                if 'best_mapq' in row:
                    del row['best_mapq']  # Remove meaningless MAPQ
                row['aligned_gene_name'] = genes[0]
                row['zero_mismatch_alignments'] = len(group)
                # Add minimum distance to 3' end for this primer
                row['distance_to_end_min'] = int(group['distance_to_end_numeric'].min())
                selected.append(row)
        else:
            primers_with_multi_gene.append((pid, genes))
    
    print(f"STAGE 3 (unique gene mapping): {len(selected)} primers retained")
    print(f"  Filtered out: {len(stage2_filtered.groupby('primer_id')) - len(selected)} primers with multi-gene hits near 3' end")
    
    # Show examples of rejected primers
    if primers_with_multi_gene and len(primers_with_multi_gene) <= 5:
        print("\n  Examples of multi-gene primers (rejected):")
        for pid, genes in primers_with_multi_gene:
            print(f"    {pid}: maps to {len(genes)} genes near 3' end: {', '.join(genes[:3])}")
    elif len(primers_with_multi_gene) > 5:
        print(f"\n  {len(primers_with_multi_gene)} primers rejected for multi-gene hits near 3' end")
    
    print(f"\n=== FINAL RESULT: {len(selected)} primers passed all filters ===\n")

    if not selected:
        # write empty with headers to be predictable
        cols = list(report.columns) + ['aligned_gene_name', 'zero_mismatch_alignments', 'distance_to_end_min']
        if 'alignment_quality' in cols:
            cols.remove('alignment_quality')
        if 'best_mapq' in cols:
            cols.remove('best_mapq')
        pd.DataFrame(columns=cols).to_csv(args.out, sep='\t', index=False)
        print(f'Wrote {args.out} (no primers passed 3-stage filtering)')
        return

    outdf = pd.DataFrame(selected)
    outdf.to_csv(args.out, sep='\t', index=False)
    print(f'Wrote {args.out} with {len(outdf)} primers')
    print(f'  Column "distance_to_end_min" shows closest distance to 3\' end (bp)')
    print(f'  Column "zero_mismatch_alignments" shows number of perfect alignments near 3\' end')


if __name__ == '__main__':
    main()
