#!/usr/bin/env python3
"""
Multi-Peak Primer Optimization

Selects the best primer combination across multiple peaks per gene.

Scoring criteria:
1. Distance to 3' end (shorter = better)
2. Isoform coverage (more distinct isoforms = better)  
3. Peak rank (lower rank number = better)

Strategy:
- For each gene with multiple peaks, score all primers across all peaks
- Select top N primers per gene (default: 3)
- Resolve isoform conflicts across genes (similar to optimize_primer_isoforms.py)

Usage:
    python bin/optimize_primers_multipeak.py \\
        --best_primers best_primers.tsv \\
        --alignment_summary primer_alignment_summary.tsv \\
        --primers_per_gene 3 \\
        --distance_weight 0.5 \\
        --isoform_weight 0.3 \\
        --peak_rank_weight 0.2 \\
        --out optimized_primers_multi.tsv
"""

import argparse
import pandas as pd
import sys
from collections import defaultdict


def parse_gene_peak_id(gene_peak_id):
    """
    Extract base gene_id and peak_rank from combined identifier.
    
    Args:
        gene_peak_id: String in format 'ENSG00000123456_peak_1' or 'ENSG00000123456'
    
    Returns:
        tuple: (base_gene_id, peak_rank)
        
    Examples:
        'ENSG00000197756_peak_1' -> ('ENSG00000197756', 1)
        'ENSG00000197756_peak_2' -> ('ENSG00000197756', 2)
        'ENSG00000197756' -> ('ENSG00000197756', 1)  # backward compatible
    """
    if '_peak_' in gene_peak_id:
        parts = gene_peak_id.rsplit('_peak_', 1)
        return parts[0], int(parts[1])
    else:
        # Backward compatible: single peak mode
        return gene_peak_id, 1


def extract_primer_info(primer_id):
    """
    Extract components from primer_id.
    
    Args:
        primer_id: String like 'ENSG00000197756_peak_1|1|F'
    
    Returns:
        dict: {gene_peak_id, base_gene_id, peak_rank, primer_index, primer_type}
    """
    parts = primer_id.split('|')
    if len(parts) < 3:
        return None
    
    gene_peak_id = parts[0]
    primer_index = parts[1]
    primer_type = parts[2]
    
    base_gene_id, peak_rank = parse_gene_peak_id(gene_peak_id)
    
    return {
        'gene_peak_id': gene_peak_id,
        'base_gene_id': base_gene_id,
        'peak_rank': peak_rank,
        'primer_index': primer_index,
        'primer_type': primer_type
    }


class MultiPeakPrimerOptimizer:
    """
    Optimizes primer selection across multiple peaks per gene.
    """
    
    def __init__(self, distance_weight=0.5, isoform_weight=0.3, peak_rank_weight=0.2):
        """
        Initialize optimizer with scoring weights.
        
        Args:
            distance_weight: Weight for distance to 3' end score (0-1)
            isoform_weight: Weight for isoform coverage score (0-1)
            peak_rank_weight: Weight for peak rank score (0-1)
        """
        self.w_dist = distance_weight
        self.w_iso = isoform_weight
        self.w_rank = peak_rank_weight
        
        # Normalize weights to sum to 1
        total = self.w_dist + self.w_iso + self.w_rank
        if total > 0:
            self.w_dist /= total
            self.w_iso /= total
            self.w_rank /= total
    
    def score_primer(self, distance, isoforms, peak_rank, max_distance, max_isoforms):
        """
        Calculate composite score for a primer.
        
        Args:
            distance: Distance to 3' end (bp)
            isoforms: Number of isoforms targeted
            peak_rank: Peak rank (1=best, 2=second, etc.)
            max_distance: Maximum distance in dataset (for normalization)
            max_isoforms: Maximum isoforms in dataset (for normalization)
        
        Returns:
            float: Composite score (0-1, higher is better)
        """
        # Distance score (closer to 3' end = better, so invert)
        if max_distance > 0:
            dist_score = 1.0 - (distance / max_distance)
            dist_score = max(0.0, min(1.0, dist_score))
        else:
            dist_score = 1.0
        
        # Isoform score (more isoforms = better)
        if max_isoforms > 0:
            iso_score = isoforms / max_isoforms
        else:
            iso_score = 0.0
        
        # Peak rank score (lower rank number = better)
        # Rank 1 = 1.0, Rank 2 = 0.5, Rank 3 = 0.33, etc.
        rank_score = 1.0 / peak_rank if peak_rank > 0 else 0.0
        
        # Weighted combination
        composite = (
            self.w_dist * dist_score +
            self.w_iso * iso_score +
            self.w_rank * rank_score
        )
        
        return composite
    
    def load_isoform_mapping(self, alignment_summary_path):
        """
        Load primer-to-isoform mapping from alignment summary.
        
        Args:
            alignment_summary_path: Path to primer_alignment_summary.tsv
        
        Returns:
            dict: {primer_id: set of transcript IDs}
        """
        try:
            summary = pd.read_csv(alignment_summary_path, sep='\t')
        except Exception as e:
            print(f"Warning: Could not load alignment summary: {e}", file=sys.stderr)
            return {}
        
        primer_isoforms = defaultdict(set)
        
        # Build mapping of primer_id -> set of transcript IDs with perfect alignment
        # Check for both possible column names for compatibility
        transcript_col = None
        if 'aligned_transcript' in summary.columns:
            transcript_col = 'aligned_transcript'
        elif 'aligned_transcript_id' in summary.columns:
            transcript_col = 'aligned_transcript_id'
        
        if 'primer_id' in summary.columns and transcript_col:
            # Filter for perfect matches (mismatches == 0)
            if 'mismatches' in summary.columns:
                perfect = summary[summary['mismatches'] == 0]
                print(f"Found {len(perfect)} perfect alignment records (0 mismatches)", file=sys.stderr)
            else:
                perfect = summary
                print("No 'mismatches' column found, using all alignments", file=sys.stderr)
            
            for _, row in perfect.iterrows():
                primer_id = row['primer_id']
                transcript_id = row.get(transcript_col, None)
                if transcript_id and pd.notna(transcript_id):
                    primer_isoforms[primer_id].add(transcript_id)
            
            print(f"Built isoform mapping for {len(primer_isoforms)} primers", file=sys.stderr)
        else:
            print(f"Warning: Required columns not found. primer_id: {'primer_id' in summary.columns}, transcript_col: {transcript_col}", file=sys.stderr)
        
        return dict(primer_isoforms)
    
    def filter_multi_gene_primers(self, alignment_summary_path, primer_ids):
        """
        Filter out primers that align to multiple distinct genes (non-specific primers).
        
        Args:
            alignment_summary_path: Path to primer_alignment_summary.tsv
            primer_ids: Set or list of primer IDs to check
        
        Returns:
            set: primer IDs that are gene-specific (single gene only)
        """
        if not alignment_summary_path:
            print("Warning: No alignment summary provided, cannot filter multi-gene primers", file=sys.stderr)
            return set(primer_ids)
        
        try:
            summary = pd.read_csv(alignment_summary_path, sep='\t')
        except Exception as e:
            print(f"Warning: Could not load alignment summary for filtering: {e}", file=sys.stderr)
            return set(primer_ids)
        
        # Check for required columns
        if 'primer_id' not in summary.columns or 'aligned_gene_name' not in summary.columns:
            print("Warning: Required columns not found for multi-gene filtering", file=sys.stderr)
            return set(primer_ids)
        
        # Filter for perfect matches only
        if 'mismatches' in summary.columns:
            summary = summary[summary['mismatches'] == 0].copy()
        
        # For each primer, count unique genes
        specific_primers = set()
        multi_gene_primers = []
        
        for pid in primer_ids:
            primer_aligns = summary[summary['primer_id'] == pid]
            if primer_aligns.empty:
                # No alignments found - keep the primer (might be a data issue)
                specific_primers.add(pid)
                continue
            
            genes = primer_aligns['aligned_gene_name'].dropna().unique()
            if len(genes) == 1:
                specific_primers.add(pid)
            else:
                multi_gene_primers.append((pid, len(genes), genes[:5]))  # Store first 5 gene names
        
        # Report filtering results
        total_primers = len(primer_ids)
        filtered_count = len(multi_gene_primers)
        
        if filtered_count > 0:
            print(f"\n=== MULTI-GENE PRIMER FILTERING ===")
            print(f"Filtered out {filtered_count}/{total_primers} primers that align to multiple genes")
            
            # Show examples
            examples = multi_gene_primers[:5]
            if examples:
                print(f"\nExamples of filtered multi-gene primers:")
                for pid, gene_count, genes in examples:
                    gene_list = ', '.join(str(g) for g in genes)
                    more = f" (+{gene_count - len(genes)} more)" if gene_count > len(genes) else ""
                    print(f"  {pid}: {gene_count} genes ({gene_list}{more})")
        else:
            print(f"\nAll {total_primers} primers are gene-specific (single gene only)")
        
        return specific_primers
    
    def optimize(self, best_primers_df, primer_isoforms, primers_per_gene=3, alignment_summary_path=None):
        """
        Select best primers per gene using greedy set cover algorithm.
        
        Strategy:
        1. Find primer with most NEW isoforms not yet covered
        2. Add to selected set
        3. Remove those isoforms from uncovered set
        4. Repeat until desired number reached or no new isoforms
        5. Ties broken by mean distance to 3' end (shorter = better)
        
        Args:
            best_primers_df: DataFrame from select_best_primer.py
            primer_isoforms: Dict of {primer_id: set of transcript IDs}
            primers_per_gene: Number of primers to select per gene
        
        Returns:
            DataFrame: Optimized primer set with composite scores
        """
        # Add primer info columns
        primer_info = best_primers_df['primer_id'].apply(extract_primer_info)
        
        # Filter out any invalid primer IDs
        valid_mask = primer_info.notna()
        if not valid_mask.all():
            print(f"Warning: {(~valid_mask).sum()} invalid primer IDs found", file=sys.stderr)
            primer_info = primer_info[valid_mask]
            best_primers_df = best_primers_df[valid_mask].copy()
        
        # Extract components
        best_primers_df['base_gene_id'] = primer_info.apply(lambda x: x['base_gene_id'] if x else None)
        best_primers_df['peak_rank'] = primer_info.apply(lambda x: x['peak_rank'] if x else 1)
        best_primers_df['gene_peak_id'] = primer_info.apply(lambda x: x['gene_peak_id'] if x else None)
        
        # Add isoform counts and sets
        best_primers_df['isoforms_targeted'] = best_primers_df['primer_id'].apply(
            lambda pid: len(primer_isoforms.get(pid, set()))
        )
        best_primers_df['isoform_set'] = best_primers_df['primer_id'].apply(
            lambda pid: primer_isoforms.get(pid, set())
        )
        
        # Ensure distance_to_end_min is numeric
        if 'distance_to_end_min' in best_primers_df.columns:
            best_primers_df['distance_to_end_min'] = pd.to_numeric(
                best_primers_df['distance_to_end_min'], 
                errors='coerce'
            ).fillna(999999)  # Large value for missing distances
        else:
            print("Warning: distance_to_end_min column not found. Using default distance=0 for all primers.", file=sys.stderr)
            best_primers_df['distance_to_end_min'] = 0
        
        # Filter out multi-gene primers
        if alignment_summary_path:
            all_primer_ids = set(best_primers_df['primer_id'].unique())
            specific_primer_ids = self.filter_multi_gene_primers(alignment_summary_path, all_primer_ids)
            best_primers_df = best_primers_df[best_primers_df['primer_id'].isin(specific_primer_ids)].copy()
        
        print(f"\n=== MULTI-PEAK PRIMER OPTIMIZATION ===")
        print(f"Total primers after filtering: {len(best_primers_df)}")
        print(f"Using greedy set cover algorithm for isoform maximization")
        print(f"Tie-breaking: mean distance to 3' end (shorter=better)")
        
        # Select top primers per gene using greedy algorithm
        results = []
        gene_stats = []
        
        for base_gene, gene_group in best_primers_df.groupby('base_gene_id'):
            # Initialize for this gene
            uncovered_isoforms = set()
            for isoform_set in gene_group['isoform_set']:
                uncovered_isoforms.update(isoform_set)
            
            total_isoforms = len(uncovered_isoforms)
            selected_primers = []
            remaining_primers = gene_group.copy()
            
            print(f"\nGene {base_gene}: {total_isoforms} isoforms, {len(gene_group)} primers")
            
            # Greedy selection
            for rank in range(1, primers_per_gene + 1):
                if len(remaining_primers) == 0:
                    print(f"  Round {rank}: No more primers available")
                    break
                
                if len(uncovered_isoforms) == 0:
                    print(f"  Round {rank}: All isoforms covered, stopping")
                    break
                
                # Calculate new isoforms each primer would add
                best_primer_idx = None
                best_new_isoforms = 0
                best_mean_dist = float('inf')
                best_primer_data = None
                
                for idx, row in remaining_primers.iterrows():
                    primer_isoforms_set = row['isoform_set']
                    new_isoforms = primer_isoforms_set & uncovered_isoforms
                    num_new = len(new_isoforms)
                    
                    # Calculate mean distance for isoforms this primer targets
                    # (used for tie-breaking)
                    mean_dist = row['distance_to_end_min']
                    
                    # Select if: more new isoforms, OR same new isoforms but closer distance
                    if (num_new > best_new_isoforms or 
                        (num_new == best_new_isoforms and mean_dist < best_mean_dist)):
                        best_primer_idx = idx
                        best_new_isoforms = num_new
                        best_mean_dist = mean_dist
                        best_primer_data = row
                
                if best_primer_idx is None or best_new_isoforms == 0:
                    print(f"  Round {rank}: No primer adds new isoforms, stopping")
                    break
                
                # Add selected primer
                selected_row = best_primer_data.copy()
                selected_row['rank_within_gene'] = rank
                selected_row['new_isoforms_added'] = best_new_isoforms
                selected_row['isoforms_covered_cumulative'] = total_isoforms - len(uncovered_isoforms) + best_new_isoforms
                selected_row['selection_reason'] = (
                    f"Rank {rank}/{len(gene_group)}: "
                    f"added {best_new_isoforms} new isoforms "
                    f"(total {selected_row['isoforms_covered_cumulative']}/{total_isoforms}), "
                    f"dist={int(best_mean_dist)}bp, peak={int(selected_row['peak_rank'])}"
                )
                
                selected_primers.append(selected_row)
                
                # Remove covered isoforms and this primer from consideration
                uncovered_isoforms -= best_primer_data['isoform_set']
                remaining_primers = remaining_primers.drop(best_primer_idx)
                
                print(f"  Round {rank}: Selected {selected_row['gene_peak_id']} "
                      f"(+{best_new_isoforms} isoforms, {len(uncovered_isoforms)} remaining)")
            
            # Calculate composite scores for selected primers
            # (for backward compatibility with reporting)
            max_dist = gene_group['distance_to_end_min'].max()
            max_iso = total_isoforms
            
            for primer_row in selected_primers:
                primer_row['composite_score'] = self.score_primer(
                    primer_row['distance_to_end_min'],
                    primer_row['isoforms_targeted'],
                    primer_row['peak_rank'],
                    max_dist,
                    max_iso
                )
                results.append(primer_row.to_dict())
            
            # Track stats
            coverage = ((total_isoforms - len(uncovered_isoforms)) / max(total_isoforms, 1)) * 100
            gene_stats.append({
                'base_gene_id': base_gene,
                'total_peaks': gene_group['peak_rank'].nunique(),
                'total_primers': len(gene_group),
                'selected_primers': len(selected_primers),
                'total_isoforms': total_isoforms,
                'isoforms_covered': total_isoforms - len(uncovered_isoforms),
                'coverage_pct': coverage
            })
        
        # Create output DataFrame
        optimized = pd.DataFrame(results)
        stats_df = pd.DataFrame(gene_stats)
        
        # Print summary
        print(f"\n=== OPTIMIZATION SUMMARY ===")
        print(f"Selected {len(optimized)} primers across {len(stats_df)} genes")
        
        # Only print detailed stats if we have results
        if len(stats_df) > 0:
            print(f"  Genes with multiple peaks: {(stats_df['total_peaks'] > 1).sum()}")
            print(f"  Average primers per gene: {len(optimized)/max(len(stats_df), 1):.1f}")
            print(f"  Average isoform coverage: {stats_df['coverage_pct'].mean():.1f}%")
            
            # Show examples
            multi_peak_genes = stats_df[stats_df['total_peaks'] > 1].head(3)
            if not multi_peak_genes.empty:
                print(f"\nExample multi-peak genes:")
                for _, gene_stat in multi_peak_genes.iterrows():
                    gene_id = gene_stat['base_gene_id']
                    gene_primers = optimized[optimized['base_gene_id'] == gene_id]
                    cov_pct = gene_stat['coverage_pct']
                    print(f"  {gene_id}: {gene_stat['total_peaks']} peaks, "
                          f"{len(gene_primers)} primers selected, "
                          f"{int(gene_stat['isoforms_covered'])}/{int(gene_stat['total_isoforms'])} isoforms ({cov_pct:.0f}%)")
        else:
            print(f"  WARNING: No primers passed filtering. All primers were either non-specific or failed quality checks.")
        
        print(f"\n=== OPTIMIZATION COMPLETE ===\n")
        
        # Add target_transcripts column (comma-separated list of transcript IDs)
        if 'isoform_set' in optimized.columns:
            optimized['target_transcripts'] = optimized['isoform_set'].apply(
                lambda isoform_set: ','.join(sorted(isoform_set)) if isoform_set else ''
            )
            # Drop temporary column after creating target_transcripts
            optimized = optimized.drop(columns=['isoform_set'])
        else:
            optimized['target_transcripts'] = ''
        
        # Reorder columns to match single peak mode format
        # Core columns that should be first (matching single peak mode)
        core_columns = [
            'gene_id', 'primer_index', 'primer_type', 'primer_sequence', 
            'gene_strand', 'rationale', 'primer_id', 'num_alignments', 
            'aligned_gene_name', 'zero_mismatch_alignments', 'distance_to_end_min',
            'isoforms_targeted', 'target_transcripts'
        ]
        
        # Multi-peak specific columns (optional, at the end)
        multipeak_columns = [
            'base_gene_id', 'peak_rank', 'gene_peak_id', 'rank_within_gene',
            'new_isoforms_added', 'isoforms_covered_cumulative', 
            'selection_reason', 'composite_score'
        ]
        
        # Build final column order: core columns + any multi-peak columns that exist
        final_columns = [col for col in core_columns if col in optimized.columns]
        final_columns += [col for col in multipeak_columns if col in optimized.columns]
        
        # Add any remaining columns not in the predefined lists
        final_columns += [col for col in optimized.columns if col not in final_columns]
        
        optimized = optimized[final_columns]
        
        return optimized


def main():
    parser = argparse.ArgumentParser(
        description='Multi-peak primer optimization with distance-based scoring'
    )
    parser.add_argument('--best_primers', required=True,
                        help='TSV file from select_best_primer.py')
    parser.add_argument('--alignment_summary', required=False,
                        help='TSV file with primer-isoform alignments (optional)')
    parser.add_argument('--out', default='optimized_primers_multi.tsv',
                        help='Output TSV file')
    parser.add_argument('--primers_per_gene', type=int, default=3,
                        help='Number of primers to select per gene (default: 3)')
    parser.add_argument('--distance_weight', type=float, default=0.5,
                        help='Weight for distance scoring (default: 0.5)')
    parser.add_argument('--isoform_weight', type=float, default=0.3,
                        help='Weight for isoform coverage (default: 0.3)')
    parser.add_argument('--peak_rank_weight', type=float, default=0.2,
                        help='Weight for peak rank (default: 0.2)')
    
    args = parser.parse_args()
    
    # Load input data
    try:
        best_primers = pd.read_csv(args.best_primers, sep='\t')
    except Exception as e:
        print(f"Error loading best_primers file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Load isoform mapping if available
    primer_isoforms = {}
    if args.alignment_summary:
        optimizer = MultiPeakPrimerOptimizer(
            distance_weight=args.distance_weight,
            isoform_weight=args.isoform_weight,
            peak_rank_weight=args.peak_rank_weight
        )
        primer_isoforms = optimizer.load_isoform_mapping(args.alignment_summary)
        print(f"Loaded isoform mappings for {len(primer_isoforms)} primers")
    else:
        print("No alignment_summary provided, isoform scoring will be 0")
        optimizer = MultiPeakPrimerOptimizer(
            distance_weight=args.distance_weight,
            isoform_weight=0.0,  # No isoform data
            peak_rank_weight=args.peak_rank_weight
        )
    
    # Run optimization
    optimized = optimizer.optimize(
        best_primers, 
        primer_isoforms, 
        args.primers_per_gene,
        alignment_summary_path=args.alignment_summary
    )
    
    # Write output (even if empty, write headers for pipeline compatibility)
    if len(optimized) == 0:
        print("\n⚠️  WARNING: No primers passed filtering!")
        print("Possible reasons:")
        print("  - All primers align to multiple genes (non-specific)")
        print("  - All primers failed quality thresholds")
        print("  - No primers found for target genes")
        print("\nCreating empty output file with headers for pipeline compatibility...")
        
    optimized.to_csv(args.out, sep='\t', index=False)
    print(f"Wrote {args.out} with {len(optimized)} optimized primers")
    
    if len(optimized) == 0:
        print("\n💡 Suggestions:")
        print("  - Check best_primers.tsv to see if primers were generated")
        print("  - Review primer_alignment_summary.tsv for multi-gene alignments")
        print("  - Consider relaxing filtering criteria or adjusting MACS2 parameters")
        sys.exit(0)  # Exit successfully to allow pipeline to continue


if __name__ == '__main__':
    main()
