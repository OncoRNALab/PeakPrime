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
            # MODIFIED: Preserve peak identifiers in multi-peak mode
            # Format: ENSG00000123456|peak_N chr1:1000-1300(+)
            # We want to keep "ENSG00000123456|peak_N" as the SEQUENCE_ID
            
            if '|' in header:
                # Split by '|' to separate gene_id and additional info
                parts = header.split('|')
                gene_id = parts[0]
                
                # Check if this is multi-peak format (has peak_N identifier)
                if len(parts) > 1 and parts[1].strip().split()[0].startswith('peak_'):
                    # Extract peak identifier (e.g., "peak_1" from "peak_1 chr1:1000-1300(+)")
                    peak_info = parts[1].strip().split()[0]  # Get "peak_1"
                    # Preserve both gene ID and peak identifier
                    sequence_id = f"{gene_id}|{peak_info}"
                else:
                    # Not multi-peak format, just use gene_id
                    sequence_id = gene_id
            elif header:
                sequence_id = header.split()[0]  # First word
            else:
                sequence_id = f"seq_{i+1}"
            
            # Write sequence identifier
            f.write(f"SEQUENCE_ID={sequence_id}\n")
            
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