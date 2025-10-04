#!/usr/bin/env python3
"""
Select optimal primer combination that maximizes isoform coverage.

This script implements a greedy algorithm to select primers that:
1. Map perfectly (0 mismatches) to only one gene (target specificity)
2. Are within a distance threshold from the transcript 3' end
3. Maximize the number of distinct isoforms covered

Algorithm:
- Filter primers by specificity and distance criteria
- Use greedy set cover: iteratively pick the primer that covers the most uncovered isoforms
- Continue until no more isoforms can be covered or max primers reached
"""

import pandas as pd
import argparse
from collections import defaultdict
import sys


def load_alignment_data(alignment_file):
    """Load primer alignment summary TSV file."""
    print(f"Loading alignment data from: {alignment_file}")
    df = pd.read_csv(alignment_file, sep='\t')
    print(f"  Loaded {len(df)} alignment records")
    print(f"  Unique genes: {df['gene_id'].nunique()}")
    print(f"  Unique primers: {df.groupby(['gene_id', 'primer_index']).ngroups}")
    return df


def filter_specific_primers(df, target_gene):
    """
    Filter primers that map perfectly to only the target gene.
    
    Returns primers where:
    - All perfect-match alignments (mismatches=0) are to the target gene
    - Ensures no off-target perfect matches
    """
    print(f"\nFiltering primers for gene: {target_gene}")
    
    # Get all alignments for this gene's primers
    gene_primers = df[df['gene_id'] == target_gene].copy()
    print(f"  Total alignments for {target_gene}: {len(gene_primers)}")
    
    # Get unique primer indices
    primer_indices = gene_primers['primer_index'].unique()
    print(f"  Unique primers: {len(primer_indices)}")
    
    specific_primers = []
    
    for primer_idx in primer_indices:
        # Get all alignments for this primer across all genes
        primer_alignments = df[
            (df['gene_id'] == target_gene) & 
            (df['primer_index'] == primer_idx)
        ]
        
        # Check perfect matches (0 mismatches)
        perfect_matches = primer_alignments[primer_alignments['mismatches'] == 0]
        
        # Get unique genes with perfect matches
        perfect_match_genes = perfect_matches['aligned_gene_name'].unique()
        
        # Keep only if ALL perfect matches are to target gene
        # (primer might align to multiple isoforms of the same gene - that's OK)
        if len(perfect_match_genes) == 1 and perfect_match_genes[0] == perfect_matches.iloc[0]['aligned_gene_name']:
            specific_primers.append(primer_idx)
    
    print(f"  Specific primers (no off-target perfect matches): {len(specific_primers)}")
    
    # Return filtered dataframe
    result = gene_primers[gene_primers['primer_index'].isin(specific_primers)]
    return result


def filter_by_distance(df, max_distance):
    """
    Filter primers within distance threshold from transcript end.
    
    Args:
        df: DataFrame with primer alignments
        max_distance: Maximum allowed distance_to_end value
    
    Returns:
        Filtered DataFrame
    """
    print(f"\nFiltering by distance to end (max={max_distance}nt)")
    before = len(df['primer_index'].unique())
    
    # Filter by distance_to_end
    result = df[df['distance_to_end'] <= max_distance].copy()
    
    after = len(result['primer_index'].unique())
    print(f"  Primers before: {before}, after: {after}, removed: {before - after}")
    
    return result


def select_optimal_primers(df, max_primers=5):
    """
    Greedy algorithm to select primer set that maximizes isoform coverage.
    
    Algorithm:
    1. Start with empty primer set and empty covered isoforms set
    2. For each remaining primer, calculate how many NEW isoforms it would cover
    3. Pick the primer that covers the most new isoforms
    4. Add to selected set and update covered isoforms
    5. Repeat until no new isoforms can be covered or max_primers reached
    
    Args:
        df: Filtered DataFrame with specific primers within distance threshold
        max_primers: Maximum number of primers to select
    
    Returns:
        List of selected primer indices and coverage statistics
    """
    print(f"\nSelecting optimal primer combination (max={max_primers} primers)")
    
    # Build primer -> isoforms mapping (only perfect matches)
    primer_to_isoforms = defaultdict(set)
    primer_info = {}
    
    for _, row in df[df['mismatches'] == 0].iterrows():
        primer_idx = row['primer_index']
        transcript = row['aligned_transcript']
        primer_to_isoforms[primer_idx].add(transcript)
        
        # Store primer metadata
        if primer_idx not in primer_info:
            primer_info[primer_idx] = {
                'sequence': row['primer_sequence'],
                'gene_id': row['gene_id'],
                'gene_name': row['aligned_gene_name'],
                'gene_strand': row['gene_strand'],
                'min_distance': row['distance_to_end']
            }
        else:
            # Keep minimum distance across isoforms
            primer_info[primer_idx]['min_distance'] = min(
                primer_info[primer_idx]['min_distance'],
                row['distance_to_end']
            )
    
    print(f"  Total candidate primers: {len(primer_to_isoforms)}")
    print(f"  Total isoforms to cover: {len(set.union(*primer_to_isoforms.values())) if primer_to_isoforms else 0}")
    
    # Greedy selection
    selected_primers = []
    covered_isoforms = set()
    all_isoforms = set.union(*primer_to_isoforms.values()) if primer_to_isoforms else set()
    
    for iteration in range(max_primers):
        if not primer_to_isoforms:
            print(f"  Iteration {iteration + 1}: No more candidate primers")
            break
        
        # Calculate marginal coverage for each remaining primer
        best_primer = None
        best_new_coverage = 0
        best_new_isoforms = set()
        
        for primer_idx, isoforms in primer_to_isoforms.items():
            new_isoforms = isoforms - covered_isoforms
            new_coverage = len(new_isoforms)
            
            if new_coverage > best_new_coverage:
                best_primer = primer_idx
                best_new_coverage = new_coverage
                best_new_isoforms = new_isoforms
        
        # Check if we can add more coverage
        if best_new_coverage == 0:
            print(f"  Iteration {iteration + 1}: No additional isoforms can be covered")
            break
        
        # Add best primer to selection
        selected_primers.append(best_primer)
        covered_isoforms.update(best_new_isoforms)
        del primer_to_isoforms[best_primer]
        
        total_isoforms = len(all_isoforms)
        covered_count = len(covered_isoforms)
        coverage_pct = (covered_count / total_isoforms * 100) if total_isoforms > 0 else 0
        
        print(f"  Iteration {iteration + 1}: Selected primer {best_primer}")
        print(f"    - Covers {best_new_coverage} new isoforms: {', '.join(sorted(best_new_isoforms))}")
        print(f"    - Total coverage: {covered_count}/{total_isoforms} ({coverage_pct:.1f}%)")
        print(f"    - Distance to end: {primer_info[best_primer]['min_distance']}nt")
    
    # Prepare results
    results = []
    for rank, primer_idx in enumerate(selected_primers, 1):
        info = primer_info[primer_idx]
        isoforms_covered = primer_to_isoforms.get(primer_idx, set())
        
        # Get isoforms covered by THIS primer (including those already covered by previous primers)
        primer_alignments = df[
            (df['primer_index'] == primer_idx) & 
            (df['mismatches'] == 0)
        ]
        all_primer_isoforms = set(primer_alignments['aligned_transcript'].unique())
        
        results.append({
            'rank': rank,
            'primer_index': primer_idx,
            'primer_sequence': info['sequence'],
            'gene_id': info['gene_id'],
            'gene_name': info['gene_name'],
            'gene_strand': info['gene_strand'],
            'min_distance_to_end': info['min_distance'],
            'isoforms_covered': len(all_primer_isoforms),
            'isoform_ids': ','.join(sorted(all_primer_isoforms))
        })
    
    return results, covered_isoforms, all_isoforms


def main():
    parser = argparse.ArgumentParser(
        description='Select optimal primer set for maximum isoform coverage',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Process single gene with default parameters
  python select_optimal_primer_set.py -i primer_alignment_summary.tsv -g ENSG00000104687 -o optimal_primers.tsv
  
  # Custom distance threshold and max primers
  python select_optimal_primer_set.py -i primer_alignment_summary.tsv -g ENSG00000104687 -d 300 -m 3 -o optimal_primers.tsv
  
  # Process all genes in the file
  python select_optimal_primer_set.py -i primer_alignment_summary.tsv --all-genes -o optimal_primers_all.tsv
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input primer alignment summary TSV file')
    parser.add_argument('-g', '--gene', 
                        help='Target gene ID (e.g., ENSG00000104687)')
    parser.add_argument('--all-genes', action='store_true',
                        help='Process all genes in the input file')
    parser.add_argument('-d', '--max-distance', type=int, default=400,
                        help='Maximum distance from transcript end (default: 400)')
    parser.add_argument('-m', '--max-primers', type=int, default=5,
                        help='Maximum number of primers to select (default: 5)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file with selected primers')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.gene and not args.all_genes:
        parser.error("Either --gene or --all-genes must be specified")
    
    # Load data
    df = load_alignment_data(args.input)
    
    # Determine which genes to process
    if args.all_genes:
        genes_to_process = df['gene_id'].unique()
        print(f"\nProcessing all genes: {len(genes_to_process)} genes found")
    else:
        genes_to_process = [args.gene]
    
    # Process each gene
    all_results = []
    
    for gene_id in genes_to_process:
        print(f"\n{'='*80}")
        print(f"Processing gene: {gene_id}")
        print('='*80)
        
        # Filter by specificity
        specific = filter_specific_primers(df, gene_id)
        
        if len(specific) == 0:
            print(f"  WARNING: No specific primers found for {gene_id}")
            continue
        
        # Filter by distance
        filtered = filter_by_distance(specific, args.max_distance)
        
        if len(filtered) == 0:
            print(f"  WARNING: No primers within distance threshold for {gene_id}")
            continue
        
        # Select optimal set
        results, covered, total = select_optimal_primers(filtered, args.max_primers)
        
        if results:
            all_results.extend(results)
            
            coverage_pct = (len(covered) / len(total) * 100) if total else 0
            print(f"\n  FINAL RESULTS for {gene_id}:")
            print(f"    Selected {len(results)} primers")
            print(f"    Coverage: {len(covered)}/{len(total)} isoforms ({coverage_pct:.1f}%)")
            print(f"    Uncovered isoforms: {', '.join(sorted(total - covered)) if total - covered else 'None'}")
        else:
            print(f"  WARNING: No primers selected for {gene_id}")
    
    # Write results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(args.output, sep='\t', index=False)
        print(f"\n{'='*80}")
        print(f"Results written to: {args.output}")
        print(f"Total primers selected: {len(all_results)}")
        print(f"Genes with primers: {results_df['gene_id'].nunique()}")
    else:
        print(f"\nWARNING: No primers selected for any gene")
        sys.exit(1)


if __name__ == '__main__':
    main()
