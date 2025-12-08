#!/usr/bin/env python3
"""
Summarize PeakPrime pipeline outputs and plot number of primers per gene after each step.

Usage:
  python summarize_primers.py <output_dir>

Requires: pandas, matplotlib
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: python summarize_primers.py <output_dir>")
    sys.exit(1)

outdir = sys.argv[1]
report_file = os.path.join(outdir, 'primer_alignment_report.tsv')
summary_file = os.path.join(outdir, 'primer_alignment_summary.tsv')
best_file = os.path.join(outdir, 'best_primers.tsv')

# Read files
report = pd.read_csv(report_file, sep='\t', dtype=str)
summary = pd.read_csv(summary_file, sep='\t', dtype=str)
best = pd.read_csv(best_file, sep='\t', dtype=str)

# Helper for unique primer count by gene_id and primer_index
primer_index_col = 'primer_index' if 'primer_index' in best.columns else 'primer_id'
gene_col = 'gene_id' if 'gene_id' in best.columns else 'gene'

# Best primers: unique primer_index per gene_id
best_counts = best.groupby(gene_col)[primer_index_col].nunique().reset_index()
best_counts.columns = ['gene_id', 'best_primers']

# Alignment summary: unique primer_index per gene_id
summary_counts = summary.groupby(gene_col)[primer_index_col].nunique().reset_index()
summary_counts.columns = ['gene_id', 'alignment_summary_primers']

# Alignment report: number of primers mapped per gene_id
report_counts = report.groupby(gene_col)[primer_index_col].nunique().reset_index()
report_counts.columns = ['gene_id', 'alignment_report_primers']

# Merge all counts
df = best_counts.merge(summary_counts, on='gene_id', how='outer')
df = df.merge(report_counts, on='gene_id', how='outer')
df = df.fillna(0)
for col in ['best_primers', 'alignment_summary_primers', 'alignment_report_primers']:
    df[col] = df[col].astype(int)



# Plot (points, integer y-ticks, improved legend, and lines connecting dots for each gene)
plt.figure(figsize=(10,6))
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5, alpha=0.7)
labels = {
    'alignment_report_primers': 'Primer3 Primers',
    'alignment_summary_primers': 'Aligned Primers',
    'best_primers': 'Best Primers'
}
steps = ['alignment_report_primers', 'alignment_summary_primers', 'best_primers']
dot_sizes = [160, 70, 30]  # Different sizes for each step
alphas = [0.7, 0.6, 0.5]   # Different transparency for each step
for i, step in enumerate(steps):
    plt.scatter(df['gene_id'], df[step], label=labels[step], s=dot_sizes[i], alpha=alphas[i])
# Add lines connecting the three dots for each gene
for idx, row in df.iterrows():
    plt.plot([row['gene_id']]*3, [row['alignment_report_primers'], row['alignment_summary_primers'], row['best_primers']], color='gray', alpha=0.5, zorder=0)
plt.xticks(rotation=45, ha='right', fontsize=6)
plt.ylabel('Number of Primers')
plt.xlabel('Gene ID')
plt.title('Number of Primers per Gene After Each Pipeline Step')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.savefig(os.path.join(outdir, 'primer_counts_per_gene.pdf'))
plt.close()
print(f"Wrote {os.path.join(outdir, 'primer_counts_per_gene.pdf')}")

# Write summary TSV
tsv_out = os.path.join(outdir, 'primer_counts_per_gene.tsv')
df.to_csv(tsv_out, sep='\t', index=False)
print(f"Wrote {tsv_out}")
