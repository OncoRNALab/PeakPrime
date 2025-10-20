#!/usr/bin/env python3
"""
Check if specific primer sequences appear in Primer3 output.

This script:
1. Extracts all designed primers from Primer3 output file
2. Reads a list of query sequences from a text file (one per line)
3. Checks which query sequences match the designed primers
4. Reports matches and mismatches

Usage:
    python check_primer_sequences.py --primer3_output primer3_output.txt --query_sequences sequences.txt --output matches.tsv
"""

import argparse
import sys
from collections import defaultdict


def parse_primer3_output(primer3_file):
    """
    Parse Primer3 output file and extract all primer sequences with position info.
    
    Note: In Primer3 output, PRIMER_LEFT_0_SEQUENCE comes BEFORE PRIMER_LEFT_0=start,length
    So we need to do a two-pass approach or store sequences temporarily.
    
    Returns:
        dict: {sequence_id: {
            'LEFT': [(seq, start, length), ...], 
            'RIGHT': [(seq, start, length), ...], 
            'INTERNAL': [(seq, start, length), ...],
            'template_length': int
        }}
    """
    primers = defaultdict(lambda: {'LEFT': {}, 'RIGHT': {}, 'INTERNAL': {}, 'template_length': None})
    current_seq_id = None
    
    with open(primer3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line == '=':
                current_seq_id = None
                continue
            
            if '=' in line:
                key, value = line.split('=', 1)
                
                if key == 'SEQUENCE_ID':
                    current_seq_id = value
                
                elif key == 'SEQUENCE_TEMPLATE' and current_seq_id:
                    # Store template length
                    primers[current_seq_id]['template_length'] = len(value)
                
                elif current_seq_id and key.startswith('PRIMER_'):
                    # Extract primer sequences
                    if '_SEQUENCE' in key:
                        if key.startswith('PRIMER_LEFT_') and '_SEQUENCE' in key:
                            idx = key.replace('PRIMER_LEFT_', '').replace('_SEQUENCE', '')
                            if idx not in primers[current_seq_id]['LEFT']:
                                primers[current_seq_id]['LEFT'][idx] = {'seq': None, 'start': None, 'length': None}
                            primers[current_seq_id]['LEFT'][idx]['seq'] = value
                        
                        elif key.startswith('PRIMER_RIGHT_') and '_SEQUENCE' in key:
                            idx = key.replace('PRIMER_RIGHT_', '').replace('_SEQUENCE', '')
                            if idx not in primers[current_seq_id]['RIGHT']:
                                primers[current_seq_id]['RIGHT'][idx] = {'seq': None, 'start': None, 'length': None}
                            primers[current_seq_id]['RIGHT'][idx]['seq'] = value
                        
                        elif key.startswith('PRIMER_INTERNAL_') and '_SEQUENCE' in key:
                            idx = key.replace('PRIMER_INTERNAL_', '').replace('_SEQUENCE', '')
                            if idx not in primers[current_seq_id]['INTERNAL']:
                                primers[current_seq_id]['INTERNAL'][idx] = {'seq': None, 'start': None, 'length': None}
                            primers[current_seq_id]['INTERNAL'][idx]['seq'] = value
                    
                    # Extract primer positions (e.g., PRIMER_LEFT_0=42,20)
                    # Check if this is a position line (key is just PRIMER_LEFT_0, not PRIMER_LEFT_0_TM, etc.)
                    elif ',' in value and not any(x in key for x in ['_TM', '_GC', '_SELF', '_HAIRPIN', '_END', '_PENALTY']):
                        if key.startswith('PRIMER_LEFT_'):
                            idx = key.replace('PRIMER_LEFT_', '')
                            if idx.isdigit():
                                start, length = map(int, value.split(','))
                                if idx not in primers[current_seq_id]['LEFT']:
                                    primers[current_seq_id]['LEFT'][idx] = {'seq': None, 'start': None, 'length': None}
                                primers[current_seq_id]['LEFT'][idx]['start'] = start
                                primers[current_seq_id]['LEFT'][idx]['length'] = length
                        
                        elif key.startswith('PRIMER_RIGHT_'):
                            idx = key.replace('PRIMER_RIGHT_', '')
                            if idx.isdigit():
                                start, length = map(int, value.split(','))
                                if idx not in primers[current_seq_id]['RIGHT']:
                                    primers[current_seq_id]['RIGHT'][idx] = {'seq': None, 'start': None, 'length': None}
                                primers[current_seq_id]['RIGHT'][idx]['start'] = start
                                primers[current_seq_id]['RIGHT'][idx]['length'] = length
                        
                        elif key.startswith('PRIMER_INTERNAL_'):
                            idx = key.replace('PRIMER_INTERNAL_', '')
                            if idx.isdigit():
                                start, length = map(int, value.split(','))
                                if idx not in primers[current_seq_id]['INTERNAL']:
                                    primers[current_seq_id]['INTERNAL'][idx] = {'seq': None, 'start': None, 'length': None}
                                primers[current_seq_id]['INTERNAL'][idx]['start'] = start
                                primers[current_seq_id]['INTERNAL'][idx]['length'] = length
    
    # Convert dict format to list format
    result = defaultdict(lambda: {'LEFT': [], 'RIGHT': [], 'INTERNAL': [], 'template_length': None})
    for seq_id, primer_data in primers.items():
        result[seq_id]['template_length'] = primer_data['template_length']
        for primer_type in ['LEFT', 'RIGHT', 'INTERNAL']:
            for idx in sorted(primer_data[primer_type].keys(), key=lambda x: int(x) if x.isdigit() else 0):
                info = primer_data[primer_type][idx]
                if info['seq'] is not None:  # Only include if we have a sequence
                    result[seq_id][primer_type].append((info['seq'], info['start'], info['length']))
    
    return result


def read_query_sequences(query_file):
    """
    Read query sequences from a text file (one sequence per line).
    Strips whitespace and ignores empty lines.
    
    Returns:
        list: List of query sequences
    """
    sequences = []
    with open(query_file, 'r') as f:
        for line in f:
            seq = line.strip()
            if seq and not seq.startswith('#'):  # Skip empty lines and comments
                sequences.append(seq.upper())  # Convert to uppercase for comparison
    return sequences


def find_matches(primers, query_sequences):
    """
    Find which query sequences match designed primers.
    
    Args:
        primers: Dict from parse_primer3_output
        query_sequences: List of query sequences
    
    Returns:
        dict: {query_seq: [(seq_id, primer_type, primer_index, start, length, template_length), ...]}
    """
    matches = defaultdict(list)
    
    # Create a lookup dictionary: sequence -> [(seq_id, type, index, start, length, template_length), ...]
    primer_lookup = defaultdict(list)
    
    for seq_id, primer_data in primers.items():
        template_length = primer_data['template_length']
        for primer_type, primer_list in primer_data.items():
            if primer_type == 'template_length':
                continue
            for idx, primer_info in enumerate(primer_list):
                if len(primer_info) == 3:
                    primer_seq, start, length = primer_info
                else:
                    primer_seq = primer_info
                    start, length = None, None
                primer_lookup[primer_seq.upper()].append((seq_id, primer_type, idx, start, length, template_length))
    
    # Check each query sequence
    for query_seq in query_sequences:
        if query_seq in primer_lookup:
            matches[query_seq] = primer_lookup[query_seq]
    
    return matches, primer_lookup


def main():
    parser = argparse.ArgumentParser(
        description='Check if query sequences match primers in Primer3 output'
    )
    parser.add_argument(
        '--primer3_output',
        required=True,
        help='Primer3 output file (e.g., primer3_output.txt)'
    )
    parser.add_argument(
        '--query_sequences',
        required=True,
        help='Text file with query sequences (one per line)'
    )
    parser.add_argument(
        '--output',
        default='primer_matches.tsv',
        help='Output file for matches (default: primer_matches.tsv)'
    )
    parser.add_argument(
        '--summary',
        action='store_true',
        help='Print summary statistics to stdout'
    )
    
    args = parser.parse_args()
    
    # Parse Primer3 output
    print(f"Parsing Primer3 output: {args.primer3_output}", file=sys.stderr)
    primers = parse_primer3_output(args.primer3_output)
    
    # Count total primers
    total_primers = 0
    total_left = 0
    total_right = 0
    total_internal = 0
    
    for seq_id, primer_data in primers.items():
        total_left += len(primer_data['LEFT'])
        total_right += len(primer_data['RIGHT'])
        total_internal += len(primer_data['INTERNAL'])
    
    total_primers = total_left + total_right + total_internal
    
    print(f"Found {total_primers} primers across {len(primers)} sequences:", file=sys.stderr)
    print(f"  - LEFT primers: {total_left}", file=sys.stderr)
    print(f"  - RIGHT primers: {total_right}", file=sys.stderr)
    print(f"  - INTERNAL primers: {total_internal}", file=sys.stderr)
    
    # Read query sequences
    print(f"\nReading query sequences: {args.query_sequences}", file=sys.stderr)
    query_sequences = read_query_sequences(args.query_sequences)
    print(f"Found {len(query_sequences)} query sequences", file=sys.stderr)
    
    # Find matches
    print("\nSearching for matches...", file=sys.stderr)
    matches, primer_lookup = find_matches(primers, query_sequences)
    
    # Write output
    with open(args.output, 'w') as out:
        # Write header
        out.write('query_sequence\tstatus\tsequence_id\tprimer_type\tprimer_index\tprimer_start\tprimer_length\ttemplate_length\tdistance_to_end\n')
        
        # Write matches
        matched_queries = set()
        for query_seq in query_sequences:
            if query_seq in matches:
                matched_queries.add(query_seq)
                for seq_id, primer_type, primer_idx, start, length, template_length in matches[query_seq]:
                    # Calculate distance to end
                    if start is not None and template_length is not None:
                        distance_to_end = template_length - start
                    else:
                        distance_to_end = 'NA'
                    
                    start_str = str(start) if start is not None else 'NA'
                    length_str = str(length) if length is not None else 'NA'
                    template_str = str(template_length) if template_length is not None else 'NA'
                    
                    out.write(f'{query_seq}\tMATCH\t{seq_id}\t{primer_type}\t{primer_idx}\t{start_str}\t{length_str}\t{template_str}\t{distance_to_end}\n')
            else:
                out.write(f'{query_seq}\tNO_MATCH\t-\t-\t-\t-\t-\t-\t-\n')
    
    print(f"\nResults written to: {args.output}", file=sys.stderr)
    
    # Print summary
    matched_count = len(matched_queries)
    unmatched_count = len(query_sequences) - matched_count
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY:", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total query sequences: {len(query_sequences)}", file=sys.stderr)
    print(f"Matched: {matched_count} ({matched_count/len(query_sequences)*100:.1f}%)", file=sys.stderr)
    print(f"Unmatched: {unmatched_count} ({unmatched_count/len(query_sequences)*100:.1f}%)", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    if args.summary:
        print("\nDETAILED MATCHES:", file=sys.stderr)
        for query_seq in sorted(matched_queries):
            print(f"\n{query_seq}:", file=sys.stderr)
            for seq_id, primer_type, primer_idx, start, length, template_length in matches[query_seq]:
                distance = template_length - start if (start is not None and template_length is not None) else 'NA'
                print(f"  - {seq_id} | {primer_type} | index {primer_idx} | start={start} | distance_to_end={distance}", file=sys.stderr)
        
        if unmatched_count > 0:
            print("\nUNMATCHED QUERIES:", file=sys.stderr)
            for query_seq in query_sequences:
                if query_seq not in matched_queries:
                    print(f"  - {query_seq}", file=sys.stderr)


if __name__ == '__main__':
    main()
