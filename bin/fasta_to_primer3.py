#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

def parse_fasta(fasta_file):
    """Parse FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences

def load_primer3_settings(settings_file):
    """Load Primer3 settings from file."""
    settings = {}
    with open(settings_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and '=' in line and not line.startswith('#'):
                key, value = line.split('=', 1)
                settings[key.strip()] = value.strip()
    return settings

def create_primer3_input(sequences, settings, output_file):
    """Create Primer3 input file from FASTA sequences and settings."""
    
    with open(output_file, 'w') as f:
        for i, (header, sequence) in enumerate(sequences):
            # Extract sequence ID from header
            # Priority: 1) Use part before '|' if present (for gene|info format)
            #          2) Use first word of header (for simple IDs like "ERCC_00002")
            #          3) Fall back to seq_N if header is empty
            if '|' in header:
                gene_id = header.split('|')[0]
            elif header:
                gene_id = header.split()[0]  # First word
            else:
                gene_id = f"seq_{i+1}"
            
            # Write sequence identifier
            f.write(f"SEQUENCE_ID={gene_id}\n")
            
            # Write sequence template
            f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
            
            # Write all settings
            for key, value in settings.items():
                f.write(f"{key}={value}\n")
            
            # Add separator
            f.write("=\n")

def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA sequences to Primer3 input format"
    )
    parser.add_argument(
        "--fasta", 
        required=True, 
        help="Input FASTA file with target sequences"
    )
    parser.add_argument(
        "--settings", 
        required=True, 
        help="Primer3 settings file"
    )
    parser.add_argument(
        "--out", 
        required=True, 
        help="Output Primer3 input file"
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file {args.fasta} not found", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.settings).exists():
        print(f"Error: Settings file {args.settings} not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Parse inputs
        sequences = parse_fasta(args.fasta)
        settings = load_primer3_settings(args.settings)
        
        if not sequences:
            print("Error: No sequences found in FASTA file", file=sys.stderr)
            sys.exit(1)
        
        # Create Primer3 input
        create_primer3_input(sequences, settings, args.out)
        
        print(f"Successfully created Primer3 input file: {args.out}")
        print(f"Processed {len(sequences)} sequences")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()