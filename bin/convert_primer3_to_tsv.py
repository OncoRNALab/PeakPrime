#!/usr/bin/env python3
"""
Convert Primer3 output to TSV format compatible with analyze_primer_alignments.py

This script reads Primer3 output and creates a TSV with columns needed for alignment analysis:
- gene_id: extracted from SEQUENCE_ID
- primer_index: primer number
- primer_type: LEFT (for our left-only primers)
- gene_strand: default to + (unknown in distance mode)
- primer_sequence: the actual primer sequence
"""

import sys
import argparse


def parse_primer3_output(primer3_file):
    """
    Parse Primer3 output file and extract primer information
    
    Returns list of primer dictionaries
    """
    primers = []
    current_seq = {}
    
    with open(primer3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line == '=':
                # End of record
                if current_seq:
                    current_seq = {}
                continue
            
            if '=' in line:
                key, value = line.split('=', 1)
                
                if key == 'SEQUENCE_ID':
                    current_seq['sequence_id'] = value
                    
                elif key.startswith('PRIMER_LEFT_') and '_SEQUENCE' in key:
                    # Extract primer number: PRIMER_LEFT_0_SEQUENCE -> 0
                    primer_num = key.replace('PRIMER_LEFT_', '').replace('_SEQUENCE', '')
                    
                    # Get gene_id from sequence_id
                    # After our changes, sequence_id should be the gene_id directly
                    gene_id = current_seq.get('sequence_id', 'unknown')
                    
                    primers.append({
                        'gene_id': gene_id,
                        'primer_index': primer_num,
                        'primer_type': 'LEFT',
                        'gene_strand': '+',  # Default - unknown in distance mode
                        'primer_sequence': value
                    })
    
    return primers


def main():
    parser = argparse.ArgumentParser(
        description='Convert Primer3 output to TSV for alignment analysis'
    )
    parser.add_argument(
        '--primer3_output',
        required=True,
        help='Input Primer3 output file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file'
    )
    parser.add_argument(
        '--gene_mapping',
        help='Optional file mapping sequence IDs to gene IDs (two columns: seq_id, gene_id)'
    )
    
    args = parser.parse_args()
    
    # Load gene mapping if provided
    gene_map = {}
    if args.gene_mapping:
        with open(args.gene_mapping, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene_map[parts[0]] = parts[1]
    
    # Parse primers
    primers = parse_primer3_output(args.primer3_output)
    
    # Apply gene mapping if available
    if gene_map:
        for p in primers:
            if p['gene_id'] in gene_map:
                p['gene_id'] = gene_map[p['gene_id']]
    
    print(f"Converted {len(primers)} primers", file=sys.stderr)
    
    # Write TSV
    with open(args.output, 'w') as out:
        # Write header - include primer_id column to match what analyze_primer_alignments.py expects
        out.write('gene_id\tprimer_index\tprimer_type\tgene_strand\tprimer_sequence\tprimer_id\n')
        
        # Write primers
        for p in primers:
            # Create primer_id in the same format as simple_primers_to_fasta.py: GENEID_primer_INDEX
            primer_id = f"{p['gene_id']}_primer_{p['primer_index']}"
            out.write(f"{p['gene_id']}\t{p['primer_index']}\t{p['primer_type']}\t"
                     f"{p['gene_strand']}\t{p['primer_sequence']}\t{primer_id}\n")


if __name__ == '__main__':
    main()
