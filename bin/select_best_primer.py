#!/usr/bin/env python3
"""
Select best primers from alignment reports.

Rules used:
- Prefer high-quality (MAPQ==255) alignments
- Require that among MAPQ==255 alignments a primer maps to a single unique gene name
- Require mismatches == 0 for those MAPQ==255 alignments

Writes out 'best_primers.tsv' in the current working directory.
"""
import argparse
import pandas as pd
import sys


def main():
    p = argparse.ArgumentParser(description='Select best primers')
    p.add_argument('--report', required=True)
    p.add_argument('--summary', required=True)
    p.add_argument('--out', default='best_primers.tsv')
    args = p.parse_args()

    report = pd.read_csv(args.report, sep='\t', dtype=str)
    try:
        summary = pd.read_csv(args.summary, sep='\t', dtype=str)
    except Exception:
                # fallback: use quick heuristic from report
        out = report[(report.get('num_alignments', '0').astype(int) == 1) &
                     (report.get('best_mapq', '0').astype(int) == 255)].copy()
        # Remove alignment_quality from fallback output too
        if 'alignment_quality' in out.columns:
            out = out.drop('alignment_quality', axis=1)
        out.to_csv(args.out, sep='	', index=False)
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
            out = report[(report.get('num_alignments', '0').astype(int) == 1) &
                         (report.get('best_mapq', '0').astype(int) == 255)].copy()
            # Remove alignment_quality from fallback output
            if 'alignment_quality' in out.columns:
                out = out.drop('alignment_quality', axis=1)
            out.to_csv(args.out, sep='\t', index=False)
            print(f'Wrote {args.out} (fallback quick filter)')
            return

    # Keep only MAPQ==255 alignments
    summary['alignment_mapq'] = pd.to_numeric(summary.get('alignment_mapq', 0), errors='coerce').fillna(0).astype(int)
    summary['mismatches'] = pd.to_numeric(summary.get('mismatches', 9999), errors='coerce').fillna(9999).astype(int)
    hq = summary[summary['alignment_mapq'] == 255].copy()

    selected = []
    for pid, group in hq.groupby('primer_id'):
        genes = group['aligned_gene_name'].dropna().unique()
        # require unique gene and all mismatches == 0
        if len(genes) == 1 and (group['mismatches'] == 0).all():
            # find report row for this primer
            rpt = report[report['primer_id'] == pid]
            if not rpt.empty:
                row = rpt.iloc[0].to_dict()
                # Remove alignment_quality from final output - keep only numerical metrics
                if 'alignment_quality' in row:
                    del row['alignment_quality']
                row['aligned_gene_name'] = genes[0]
                row['hq_alignments'] = len(group)
                row['zero_mismatch_alignments'] = len(group[group['mismatches'] == 0])
                selected.append(row)

    if not selected:
        # write empty with headers to be predictable
        cols = list(report.columns) + ['aligned_gene_name', 'hq_alignments']
        pd.DataFrame(columns=cols).to_csv(args.out, sep='\t', index=False)
        print(f'Wrote {args.out} (no strict matches found)')
        return

    outdf = pd.DataFrame(selected)
    outdf.to_csv(args.out, sep='\t', index=False)
    print(f'Wrote {args.out} with {len(outdf)} entries')


if __name__ == '__main__':
    main()
