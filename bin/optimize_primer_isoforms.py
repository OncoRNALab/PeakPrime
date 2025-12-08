#!/usr/bin/env python3
"""Optimize primer selection to maximize distinct isoform coverage.

Reads best_primers.tsv and primer_alignment_summary.tsv, then selects
one optimal primer per gene to maximize the total number of unique
transcripts (isoforms) targeted, ensuring no isoform is targeted by
multiple primers.

Strategy:
1. Within each gene: Select primer targeting most isoforms
2. Across genes: Resolve conflicts where multiple primers target same isoform
3. Output: One primer per gene with maximal unique isoform coverage

Note: Uses same max_mismatches threshold as select_best_primer.py to ensure
      consistent primer specificity filtering.
"""

import argparse
import pandas as pd
import sys
from typing import Dict, Set, List, Tuple
from collections import defaultdict


def load_data(best_primers_file: str, alignment_summary_file: str, 
              distance_threshold: float = 1000.0,
              max_mismatches: int = 0) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate input files.
    
    Args:
        max_mismatches: Maximum mismatches for primer specificity (0-3).
                       Must match value used in select_best_primer.py.
    """
    print(f"\n{'='*60}")
    print("LOADING INPUT DATA")
    print(f"{'='*60}")
    
    # Load best primers
    best_primers = pd.read_csv(best_primers_file, sep='\t', dtype=str)
    print(f"✓ Loaded {len(best_primers)} primers from {best_primers_file}")
    print(f"  Genes: {best_primers['gene_id'].nunique()}")
    
    # Load alignment summary
    try:
        alignment_summary = pd.read_csv(alignment_summary_file, sep='\t', dtype=str)
        print(f"✓ Loaded {len(alignment_summary)} alignments from {alignment_summary_file}")
    except FileNotFoundError:
        print(f"✗ Alignment summary file not found: {alignment_summary_file}")
        print("  Cannot perform isoform optimization without transcript information.")
        sys.exit(1)
    
    # Validate required columns
    required_best = ['gene_id', 'primer_index']
    required_align = ['primer_index', 'aligned_transcript', 'aligned_gene_name', 'mismatches']
    
    missing_best = [col for col in required_best if col not in best_primers.columns]
    missing_align = [col for col in required_align if col not in alignment_summary.columns]
    
    if missing_best or missing_align:
        if missing_best:
            print(f"✗ Missing columns in best_primers: {missing_best}")
        if missing_align:
            print(f"✗ Missing columns in alignment_summary: {missing_align}")
        sys.exit(1)
    
    # Filter alignment summary to only best primers
    # CRITICAL: Must match on BOTH gene_id AND primer_index, as primer_index alone is not unique
    # Create composite key for filtering
    best_primers_keys = set(
        zip(best_primers['gene_id'], best_primers['primer_index'].astype(str))
    )
    alignment_summary['composite_key'] = list(
        zip(alignment_summary['gene_id'], alignment_summary['primer_index'].astype(str))
    )
    alignment_summary = alignment_summary[
        alignment_summary['composite_key'].isin(best_primers_keys)
    ].copy()
    alignment_summary = alignment_summary.drop('composite_key', axis=1)
    print(f"  Filtered to {len(alignment_summary)} alignments for validated primers")
    
    # Apply quality filters to alignments
    # Convert numeric columns
    alignment_summary['mismatches'] = pd.to_numeric(
        alignment_summary['mismatches'], errors='coerce'
    ).fillna(9999).astype(int)
    
    if 'distance_to_end' in alignment_summary.columns:
        alignment_summary['distance_to_end'] = pd.to_numeric(
            alignment_summary['distance_to_end'], errors='coerce'
        ).fillna(float('inf'))
    
    # Filter: within mismatch threshold
    initial_count = len(alignment_summary)
    alignment_summary = alignment_summary[alignment_summary['mismatches'] <= max_mismatches]
    print(f"  Applied quality filters:")
    print(f"    • Alignments (mismatches≤{max_mismatches}): {len(alignment_summary)}/{initial_count}")
    
    # Filter: within distance threshold (if distance data available)
    if 'distance_to_end' in alignment_summary.columns:
        before = len(alignment_summary)
        alignment_summary = alignment_summary[
            alignment_summary['distance_to_end'] <= distance_threshold
        ]
        print(f"    • Distance ≤ {distance_threshold}bp: {len(alignment_summary)}/{before}")
    
    return best_primers, alignment_summary


def build_primer_isoform_mapping(best_primers: pd.DataFrame, 
                                  alignment_summary: pd.DataFrame) -> Dict[Tuple[str, int], Set[str]]:
    """
    Build mapping of (gene_id, primer_index) -> set of target transcript IDs.
    Only includes on-target transcripts (gene_id from both sources must match).
    """
    print(f"\n{'='*60}")
    print("BUILDING PRIMER-ISOFORM MAPPING")
    print(f"{'='*60}")
    
    # Create composite key mapping: primer_index -> gene_id
    # Since primer_index alone is not unique, we need gene_id context
    primer_to_gene = dict(zip(best_primers['primer_index'].astype(str), best_primers['gene_id']))
    
    # Also need to match gene_id from alignment_summary
    # Ensure alignment_summary has gene_id column (should be there)
    if 'gene_id' not in alignment_summary.columns:
        print(f"✗ Missing 'gene_id' column in alignment_summary")
        print(f"  Available columns: {list(alignment_summary.columns)}")
        return {}
    
    primer_isoforms = defaultdict(set)
    
    for _, row in alignment_summary.iterrows():
        gene_id = row['gene_id']
        primer_index = str(row['primer_index'])
        transcript_id = row['aligned_transcript']
        
        # Skip if missing data
        if pd.isna(gene_id) or pd.isna(primer_index) or pd.isna(transcript_id):
            continue
        
        # Create composite key (gene_id, primer_index)
        composite_key = (gene_id, int(primer_index))
        
        # Add transcript to this primer's set
        # Note: We already filtered to gene_id matches in alignment_summary
        # and applied quality filters (mismatches≤threshold, distance threshold)
        primer_isoforms[composite_key].add(transcript_id)
    
    # Convert to regular dict
    primer_isoforms = dict(primer_isoforms)
    
    # Statistics
    primers_with_isoforms = len(primer_isoforms)
    total_isoforms = sum(len(isoforms) for isoforms in primer_isoforms.values())
    avg_isoforms = total_isoforms / primers_with_isoforms if primers_with_isoforms > 0 else 0
    
    print(f"✓ Mapped {primers_with_isoforms} primers to target isoforms")
    print(f"  Total distinct primer-isoform pairs: {total_isoforms}")
    print(f"  Average isoforms per primer: {avg_isoforms:.1f}")
    
    # Show primers with no isoform data
    # Create composite keys from best_primers for comparison
    best_primers_keys = set(zip(best_primers['gene_id'], best_primers['primer_index'].astype(int)))
    primers_without_isoforms = best_primers_keys - set(primer_isoforms.keys())
    if primers_without_isoforms:
        print(f"  ⚠ {len(primers_without_isoforms)} primers have no isoform data")
    
    return primer_isoforms


def select_best_primer_per_gene(best_primers: pd.DataFrame,
                                 primer_isoforms: Dict[Tuple[str, int], Set[str]]) -> Dict[str, Tuple[str, int]]:
    """
    For each gene, select the primer that targets the most distinct isoforms.
    Returns: {gene_id: (gene_id, primer_index)}
    """
    print(f"\n{'='*60}")
    print("PHASE 1: WITHIN-GENE OPTIMIZATION")
    print(f"{'='*60}")
    
    # Group primers by gene
    gene_primers = best_primers.groupby('gene_id')['primer_index'].apply(lambda x: x.astype(int).tolist()).to_dict()
    
    selected_primers = {}
    multi_primer_genes = 0
    
    for gene_id, primer_indices in gene_primers.items():
        if len(primer_indices) == 1:
            # Only one primer for this gene - automatically select it
            selected_primers[gene_id] = (gene_id, primer_indices[0])
        else:
            multi_primer_genes += 1
            # Multiple primers - select one with most isoforms
            primer_scores = []
            for primer_index in primer_indices:
                composite_key = (gene_id, primer_index)
                isoforms = primer_isoforms.get(composite_key, set())
                isoform_count = len(isoforms)
                
                # Get additional metrics for tie-breaking
                primer_data = best_primers[
                    (best_primers['gene_id'] == gene_id) & 
                    (best_primers['primer_index'].astype(int) == primer_index)
                ].iloc[0]
                distance = float(primer_data.get('distance_to_end_min', 999999))
                num_alignments = int(primer_data.get('zero_mismatch_alignments', 0))
                
                primer_scores.append((
                    isoform_count,      # Primary: maximize isoforms
                    -distance,          # Secondary: minimize distance (closer to 3' end)
                    -num_alignments,    # Tertiary: minimize alignments (more specific)
                    primer_index        # Quaternary: deterministic tie-breaker
                ))
            
            # Select primer with best score (reverse sort, take first)
            primer_scores.sort(reverse=True)
            best_primer_index = primer_scores[0][3]  # primer_index is 4th element
            selected_primers[gene_id] = (gene_id, best_primer_index)
            
            if len(set(score[0] for score in primer_scores)) == 1:
                # All primers have same isoform count - used tie-breakers
                pass
    
    print(f"✓ Selected {len(selected_primers)} primers (one per gene)")
    print(f"  Genes with multiple primers: {multi_primer_genes}")
    print(f"  Genes with single primer: {len(gene_primers) - multi_primer_genes}")
    
    # Count isoforms in initial selection
    total_isoforms = set()
    for composite_key in selected_primers.values():
        total_isoforms.update(primer_isoforms.get(composite_key, set()))
    
    print(f"  Initial total distinct isoforms: {len(total_isoforms)}")
    
    return selected_primers


def resolve_isoform_conflicts(selected_primers: Dict[str, Tuple[str, int]],
                               primer_isoforms: Dict[Tuple[str, int], Set[str]],
                               best_primers: pd.DataFrame,
                               gene_primers: Dict[str, List[int]]) -> Dict[str, Tuple[str, int]]:
    """
    Iteratively resolve conflicts where multiple primers target the same isoform.
    Returns updated {gene_id: (gene_id, primer_index)} mapping with no shared isoforms.
    """
    print(f"\n{'='*60}")
    print("PHASE 2: CROSS-GENE DEDUPLICATION")
    print(f"{'='*60}")
    
    # Build gene_primers mapping if not provided
    if not gene_primers:
        gene_primers = best_primers.groupby('gene_id')['primer_index'].apply(lambda x: x.astype(int).tolist()).to_dict()
    
    iteration = 0
    max_iterations = 100  # Safety limit
    
    while iteration < max_iterations:
        iteration += 1
        
        # Build isoform -> primers mapping for current selection
        isoform_to_primers = defaultdict(set)
        for gene_id, composite_key in selected_primers.items():
            isoforms = primer_isoforms.get(composite_key, set())
            for isoform in isoforms:
                isoform_to_primers[isoform].add((gene_id, composite_key))
        
        # Find conflicting isoforms (targeted by multiple primers)
        conflicts = {
            isoform: primers 
            for isoform, primers in isoform_to_primers.items() 
            if len(primers) > 1
        }
        
        if not conflicts:
            print(f"\n✓ No conflicts remaining after {iteration-1} iteration(s)")
            break
        
        if iteration == 1:
            print(f"  Found {len(conflicts)} conflicting isoforms")
        
        # Resolve conflicts: for each conflicting isoform, keep best primer
        primers_to_remove = set()
        
        for isoform, primers_set in conflicts.items():
            # Score each primer in the conflict
            primer_scores = []
            for gene_id, composite_key in primers_set:
                isoform_count = len(primer_isoforms.get(composite_key, set()))
                
                # Get additional metrics
                gene_id_key, primer_index_key = composite_key
                primer_data = best_primers[
                    (best_primers['gene_id'] == gene_id_key) & 
                    (best_primers['primer_index'].astype(int) == primer_index_key)
                ].iloc[0]
                distance = float(primer_data.get('distance_to_end_min', 999999))
                
                # Check how many alternative primers this gene has
                alternatives = len(gene_primers.get(gene_id, []))
                
                primer_scores.append((
                    isoform_count,        # Primary: maximize total isoforms
                    -alternatives,        # Secondary: genes with fewer alternatives have priority
                    -distance,            # Tertiary: closer to 3' end
                    gene_id,              # Quaternary: deterministic
                    composite_key
                ))
            
            # Sort and keep best, remove others
            primer_scores.sort(reverse=True)
            winners = [primer_scores[0]]  # Keep the best
            
            # Remove all except winner
            for score in primer_scores[1:]:
                gene_id, composite_key = score[3], score[4]
                primers_to_remove.add(gene_id)
        
        if not primers_to_remove:
            # No changes made, break to avoid infinite loop
            print(f"  ⚠ Could not resolve conflicts (iteration {iteration})")
            break
        
        # Remove losing primers and try to replace with alternatives
        genes_without_primers = set()
        
        for gene_id in primers_to_remove:
            del selected_primers[gene_id]
            
            # Try to select alternative primer for this gene
            available_primers = gene_primers.get(gene_id, [])
            already_tried = set()
            
            for primer_index in available_primers:
                if primer_index in already_tried:
                    continue
                
                already_tried.add(primer_index)
                composite_key = (gene_id, primer_index)
                
                # Check if this primer's isoforms conflict with currently selected
                primer_isoforms_set = primer_isoforms.get(composite_key, set())
                
                # Check for conflicts with other selected primers
                has_conflict = False
                for other_gene, other_composite_key in selected_primers.items():
                    if other_gene == gene_id:
                        continue
                    other_isoforms = primer_isoforms.get(other_composite_key, set())
                    if primer_isoforms_set & other_isoforms:  # Intersection
                        has_conflict = True
                        break
                
                if not has_conflict:
                    selected_primers[gene_id] = composite_key
                    break
            
            if gene_id not in selected_primers:
                genes_without_primers.add(gene_id)
        
        print(f"  Iteration {iteration}: Removed {len(primers_to_remove)} primers, "
              f"{len(genes_without_primers)} genes without alternatives")
    
    if iteration >= max_iterations:
        print(f"  ⚠ Reached maximum iterations ({max_iterations})")
    
    # Final statistics
    final_isoforms = set()
    for composite_key in selected_primers.values():
        final_isoforms.update(primer_isoforms.get(composite_key, set()))
    
    print(f"\n✓ Final selection:")
    print(f"  Genes with primers: {len(selected_primers)}")
    print(f"  Total distinct isoforms: {len(final_isoforms)}")
    
    return selected_primers


def create_output(best_primers: pd.DataFrame,
                  selected_primers: Dict[str, Tuple[str, int]],
                  primer_isoforms: Dict[Tuple[str, int], Set[str]],
                  output_file: str):
    """Create optimized output file with isoform information."""
    print(f"\n{'='*60}")
    print("GENERATING OUTPUT")
    print(f"{'='*60}")
    
    # Extract (gene_id, primer_index) tuples from selected_primers
    selected_keys = set(selected_primers.values())
    
    # Filter best_primers to only selected primers
    output_rows = []
    for gene_id, primer_index in selected_keys:
        matching_rows = best_primers[
            (best_primers['gene_id'] == gene_id) & 
            (best_primers['primer_index'].astype(int) == primer_index)
        ]
        if not matching_rows.empty:
            output_rows.append(matching_rows.iloc[0])
    
    output_df = pd.DataFrame(output_rows)
    
    # Add isoform information
    output_df['isoforms_targeted'] = output_df.apply(
        lambda row: len(primer_isoforms.get((row['gene_id'], int(row['primer_index'])), set())),
        axis=1
    )
    
    output_df['target_transcripts'] = output_df.apply(
        lambda row: ','.join(sorted(primer_isoforms.get((row['gene_id'], int(row['primer_index'])), set()))),
        axis=1
    )
    
    # Sort by gene_id for consistency
    output_df = output_df.sort_values('gene_id')
    
    # Write output
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Written {len(output_df)} optimized primers to: {output_file}")
    
    # Summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Input primers: {len(best_primers)}")
    print(f"Input genes: {best_primers['gene_id'].nunique()}")
    print(f"Output primers: {len(output_df)}")
    print(f"Output genes: {output_df['gene_id'].nunique()}")
    print(f"Genes removed: {best_primers['gene_id'].nunique() - output_df['gene_id'].nunique()}")
    
    total_isoforms = output_df['isoforms_targeted'].sum()
    unique_isoforms = len(set().union(*[primer_isoforms.get((row['gene_id'], int(row['primer_index'])), set()) 
                                        for _, row in output_df.iterrows()]))
    
    print(f"\nIsoform coverage:")
    print(f"  Total primer-isoform pairs: {total_isoforms}")
    print(f"  Distinct isoforms targeted: {unique_isoforms}")
    if len(output_df) > 0:
        print(f"  Average isoforms per primer: {total_isoforms / len(output_df):.1f}")
    
    # Check for any remaining conflicts (should be zero)
    all_isoforms_list = []
    for _, row in output_df.iterrows():
        composite_key = (row['gene_id'], int(row['primer_index']))
        all_isoforms_list.extend(primer_isoforms.get(composite_key, set()))
    
    if len(all_isoforms_list) != len(set(all_isoforms_list)):
        duplicates = len(all_isoforms_list) - len(set(all_isoforms_list))
        print(f"\n⚠ WARNING: {duplicates} isoform conflicts remain!")
    else:
        print(f"\n✓ VERIFIED: All isoforms are unique (no conflicts)")


def main():
    parser = argparse.ArgumentParser(
        description='Optimize primer selection to maximize distinct isoform coverage',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--best_primers', required=True,
                       help='Input file: best_primers.tsv (validated primers)')
    parser.add_argument('--alignment_summary', required=True,
                       help='Input file: primer_alignment_summary.tsv (with transcript IDs)')
    parser.add_argument('--output', default='best_primers_optimal.tsv',
                       help='Output file: optimized primer selection (default: best_primers_optimal.tsv)')
    parser.add_argument('--distance_threshold', type=float, default=1000.0,
                       help='Maximum distance from 3\' end for isoform consideration (bp). Default: 1000')
    parser.add_argument('--max_mismatches', type=int, default=0,
                       help='Maximum mismatches for primer specificity (0-3). Must match select_best_primer.py. Default: 0')
    
    args = parser.parse_args()
    
    # Load data
    best_primers, alignment_summary = load_data(
        args.best_primers, 
        args.alignment_summary,
        args.distance_threshold,
        args.max_mismatches
    )
    
    # Build primer-isoform mapping
    primer_isoforms = build_primer_isoform_mapping(best_primers, alignment_summary)
    
    # Check if we have isoform data
    if not primer_isoforms:
        print("\n✗ No isoform data available. Cannot perform optimization.")
        print("  Copying best_primers.tsv to output (no optimization applied)")
        best_primers.to_csv(args.output, sep='\t', index=False)
        return
    
    # Phase 1: Select best primer per gene (within-gene optimization)
    selected_primers = select_best_primer_per_gene(best_primers, primer_isoforms)
    
    # Build gene_primers mapping for conflict resolution
    gene_primers = best_primers.groupby('gene_id')['primer_index'].apply(lambda x: x.astype(int).tolist()).to_dict()
    
    # Phase 2: Resolve cross-gene isoform conflicts
    selected_primers = resolve_isoform_conflicts(
        selected_primers, 
        primer_isoforms,
        best_primers,
        gene_primers
    )
    
    # Create output
    create_output(best_primers, selected_primers, primer_isoforms, args.output)
    
    print(f"\n{'='*60}")
    print("OPTIMIZATION COMPLETE")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
