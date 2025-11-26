#!/usr/bin/env python3
"""
Find primers by sequence in best_primers.tsv file.

Given a text file with primer sequences (one per line) and a best_primers.tsv file,
reports which sequences are found and their associated gene information.

Usage:
    python find_primers_by_sequence.py \\
        --primers best_primers.tsv \\
        --sequences sequences.txt \\
        --out matches.tsv

Output:
    TSV file with columns: query_sequence, gene_id, peak_id, primer_index, primer_sequence
    Also prints summary to stdout.
"""

import argparse
import pandas as pd
import sys


def load_query_sequences(sequences_file):
    """
    Load sequences from text file (one per line).
    
    Args:
        sequences_file: Path to text file with one sequence per line
    
    Returns:
        list: List of sequences (stripped and uppercased)
    """
    sequences = []
    with open(sequences_file, 'r') as f:
        for line in f:
            seq = line.strip().upper()
            if seq:  # Skip empty lines
                sequences.append(seq)
    return sequences


def find_matching_primers(primers_df, query_sequences):
    """
    Find primers in dataframe that match query sequences.
    
    Args:
        primers_df: DataFrame from best_primers.tsv
        query_sequences: List of sequences to search for
    
    Returns:
        DataFrame: Matched primers with query_sequence column added
    """
    # Normalize primer_sequence column to uppercase for comparison
    primers_df = primers_df.copy()
    primers_df['primer_sequence_upper'] = primers_df['primer_sequence'].str.upper()
    
    # Create a set of query sequences for fast lookup
    query_set = set(query_sequences)
    
    # Filter primers that match any query sequence
    matched = primers_df[primers_df['primer_sequence_upper'].isin(query_set)].copy()
    
    # Add the original query sequence (in case of mixed case in input)
    matched['query_sequence'] = matched['primer_sequence_upper']
    
    # Drop temporary uppercase column
    matched = matched.drop(columns=['primer_sequence_upper'])
    
    return matched


def main():
    parser = argparse.ArgumentParser(
        description='Find primers by sequence in best_primers.tsv file'
    )
    parser.add_argument('--primers', required=True,
                        help='Path to best_primers.tsv file')
    parser.add_argument('--sequences', required=True,
                        help='Path to text file with sequences (one per line)')
    parser.add_argument('--out', default='primer_matches.tsv',
                        help='Output TSV file (default: primer_matches.tsv)')
    
    args = parser.parse_args()
    
    # Load input files
    print(f"Loading primers from: {args.primers}")
    try:
        primers_df = pd.read_csv(args.primers, sep='\t')
    except Exception as e:
        print(f"Error loading primers file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading query sequences from: {args.sequences}")
    try:
        query_sequences = load_query_sequences(args.sequences)
    except Exception as e:
        print(f"Error loading sequences file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nTotal primers in file: {len(primers_df)}")
    print(f"Total query sequences: {len(query_sequences)}")
    print(f"Unique query sequences: {len(set(query_sequences))}")
    
    # Find matches
    matched = find_matching_primers(primers_df, query_sequences)
    
    print(f"\n{'='*60}")
    print(f"MATCHES FOUND: {len(matched)}")
    print(f"{'='*60}")
    
    if len(matched) == 0:
        print("\n❌ No matches found!")
        print("\nPossible reasons:")
        print("  - Sequences are not in the primers file")
        print("  - Case mismatch (script normalizes to uppercase)")
        print("  - Whitespace differences (script strips whitespace)")
        print("\nCreating empty output file...")
        
        # Create empty output with headers
        output_columns = ['query_sequence', 'gene_id', 'peak_id', 'primer_index', 'primer_sequence']
        pd.DataFrame(columns=output_columns).to_csv(args.out, sep='\t', index=False)
        print(f"Wrote empty output to: {args.out}")
        return
    
    # Select and order output columns
    output_columns = ['query_sequence', 'gene_id', 'peak_id', 'primer_index', 'primer_sequence']
    
    # Check which columns exist in the matched dataframe
    available_columns = [col for col in output_columns if col in matched.columns]
    
    # Add any additional columns not in the predefined list
    additional_columns = [col for col in matched.columns if col not in output_columns]
    final_columns = available_columns + additional_columns
    
    output_df = matched[final_columns]
    
    # Write output
    output_df.to_csv(args.out, sep='\t', index=False)
    print(f"\n✅ Wrote {len(output_df)} matches to: {args.out}")
    
    # Print detailed summary
    print(f"\n{'='*60}")
    print("MATCH DETAILS")
    print(f"{'='*60}")
    
    for idx, row in output_df.iterrows():
        gene_id = row.get('gene_id', 'N/A')
        peak_id = row.get('peak_id', 'N/A')
        primer_index = row.get('primer_index', 'N/A')
        sequence = row.get('primer_sequence', 'N/A')
        
        print(f"\nGene: {gene_id}")
        print(f"  Peak ID: {peak_id}")
        print(f"  Primer index: {primer_index}")
        print(f"  Sequence: {sequence}")
    
    # Check for query sequences that were NOT found
    matched_sequences = set(matched['query_sequence'].str.upper())
    query_sequences_upper = set(seq.upper() for seq in query_sequences)
    not_found = query_sequences_upper - matched_sequences
    
    if not_found:
        print(f"\n{'='*60}")
        print(f"SEQUENCES NOT FOUND: {len(not_found)}")
        print(f"{'='*60}")
        for seq in sorted(not_found):
            print(f"  {seq}")
    else:
        print(f"\n✅ All query sequences were found!")


if __name__ == '__main__':
    main()
