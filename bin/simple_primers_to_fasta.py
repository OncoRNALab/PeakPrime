#!/usr/bin/env python3
"""
Simple conversion of Primer3 output to FASTA format for distance-based primer design.
Unlike primers_to_fasta.R, this doesn't require gene information.
"""

import argparse
import sys
import re


def parse_primer3_output(primer3_file):
    """
    Parse Primer3 output file and extract primer sequences.
    
    Returns:
        list of dicts with primer information
    """
    primers = []
    current_seq = {}
    
    with open(primer3_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line or line == '=':
                if current_seq:
                    primers.append(current_seq)
                    current_seq = {}
                continue
            
            if '=' in line:
                key, value = line.split('=', 1)
                current_seq[key] = value
    
    # Don't forget the last sequence
    if current_seq:
        primers.append(current_seq)
    
    return primers


def extract_primers_to_fasta(primer3_file, output_fasta):
    """
    Extract primer sequences from Primer3 output and write to FASTA.
    """
    sequences = parse_primer3_output(primer3_file)
    primer_count = 0
    
    with open(output_fasta, 'w') as out:
        for seq_dict in sequences:
            seq_id = seq_dict.get('SEQUENCE_ID', 'unknown')
            
            # Extract all LEFT primers for this sequence
            i = 0
            while f'PRIMER_LEFT_{i}_SEQUENCE' in seq_dict:
                primer_seq = seq_dict[f'PRIMER_LEFT_{i}_SEQUENCE']
                primer_penalty = seq_dict.get(f'PRIMER_LEFT_{i}_PENALTY', 'NA')
                primer_tm = seq_dict.get(f'PRIMER_LEFT_{i}_TM', 'NA')
                primer_gc = seq_dict.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 'NA')
                
                # Write FASTA entry
                header = f">{seq_id}_primer_{i} penalty={primer_penalty} tm={primer_tm} gc={primer_gc}"
                out.write(f"{header}\n")
                out.write(f"{primer_seq}\n")
                
                primer_count += 1
                i += 1
    
    return primer_count


def main():
    parser = argparse.ArgumentParser(
        description='Convert Primer3 output to FASTA format (simple version for distance mode)'
    )
    parser.add_argument(
        '--primer3-output',
        required=True,
        help='Input Primer3 output file'
    )
    parser.add_argument(
        '--output-fasta',
        required=True,
        help='Output FASTA file with primer sequences'
    )
    
    args = parser.parse_args()
    
    print(f"Reading Primer3 output from: {args.primer3_output}", file=sys.stderr)
    
    try:
        primer_count = extract_primers_to_fasta(args.primer3_output, args.output_fasta)
        print(f"Extracted {primer_count} primers to {args.output_fasta}", file=sys.stderr)
        
        if primer_count == 0:
            print("WARNING: No primers were extracted", file=sys.stderr)
            # Create empty file so workflow can continue
            with open(args.output_fasta, 'w') as f:
                f.write("")
        
        print("Success!", file=sys.stderr)
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
