#!/usr/bin/env python3
"""
Extract N bases from the 3' end of transcript sequences in a FASTA file.
This creates primer design templates for 3' end-targeted primers.
"""

import argparse
import sys
from typing import Iterator, Tuple


def parse_fasta(fasta_file: str) -> Iterator[Tuple[str, str, str]]:
    """
    Parse FASTA file and yield (header, sequence_id, sequence) tuples.
    
    Args:
        fasta_file: Path to FASTA file
        
    Yields:
        Tuples of (full_header, sequence_id, sequence)
    """
    header = None
    sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            if line.startswith('>'):
                # If we have a previous sequence, yield it
                if header is not None:
                    seq_id = header.split()[0]
                    yield (header, seq_id, ''.join(sequence))
                
                # Start new sequence
                header = line[1:]  # Remove '>'
                sequence = []
            else:
                sequence.append(line)
        
        # Don't forget the last sequence
        if header is not None:
            seq_id = header.split()[0]
            yield (header, seq_id, ''.join(sequence))


def extract_3prime(sequence: str, length: int) -> str:
    """
    Extract the last N bases from a sequence.
    
    Args:
        sequence: DNA sequence string
        length: Number of bases to extract from 3' end
        
    Returns:
        3' end sequence of specified length (or entire sequence if shorter)
    """
    if len(sequence) <= length:
        return sequence
    else:
        return sequence[-length:]


def main():
    parser = argparse.ArgumentParser(
        description='Extract N bases from the 3\' end of sequences in a FASTA file'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output FASTA file with 3\' sequences'
    )
    parser.add_argument(
        '--length',
        type=int,
        required=True,
        help='Number of bases to extract from 3\' end'
    )
    
    args = parser.parse_args()
    
    if args.length <= 0:
        print(f"ERROR: Length must be positive (got {args.length})", file=sys.stderr)
        sys.exit(1)
    
    print(f"Extracting {args.length} bases from 3' end of sequences...", file=sys.stderr)
    
    processed_count = 0
    skipped_count = 0
    
    with open(args.output, 'w') as out:
        for header, seq_id, sequence in parse_fasta(args.input):
            if not sequence:
                print(f"Warning: Empty sequence for {seq_id}, skipping", file=sys.stderr)
                skipped_count += 1
                continue
            
            # Extract 3' end
            threeprime_seq = extract_3prime(sequence, args.length)
            
            # Write to output with updated header indicating 3' extraction
            original_len = len(sequence)
            extracted_len = len(threeprime_seq)
            
            out.write(f">{header} 3prime_extracted={extracted_len}bp original_length={original_len}bp\n")
            out.write(f"{threeprime_seq}\n")
            
            processed_count += 1
            
            if processed_count <= 5:  # Log first few for debugging
                print(f"  {seq_id}: {original_len}bp -> {extracted_len}bp (3' end)", file=sys.stderr)
    
    print(f"\nCompleted: Processed {processed_count} sequences", file=sys.stderr)
    if skipped_count > 0:
        print(f"  Skipped {skipped_count} empty sequences", file=sys.stderr)
    
    if processed_count == 0:
        print("ERROR: No sequences were processed", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
