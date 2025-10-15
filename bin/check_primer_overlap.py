#!/usr/bin/env python3
"""
Check if primers (aligned to genome) overlap with genomic regions.
"""

import re
from Bio import SeqIO

def parse_fasta_regions(fasta_file):
    """
    Parse FASTA file and extract genomic regions.
    Returns dict: {gene_id: {'chrom': str, 'start': int, 'end': int, 'strand': str, 'sequence': str}}
    
    Header format: >ENSG00000068383|10:132783177-132783301(+)
    """
    regions = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.id
        parts = header.split('|')
        gene_id = parts[0]
        
        # Extract coordinates and strand from format: 10:132783177-132783301(+)
        coord_match = re.match(r'(\w+):(\d+)-(\d+)(\([+-]\))', parts[1])
        if coord_match:
            chrom = coord_match.group(1)
            start = int(coord_match.group(2))
            end = int(coord_match.group(3))
            strand = coord_match.group(4).strip('()')
            
            regions[gene_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'sequence': str(record.seq)
            }
    return regions

def parse_sam_file(sam_file):
    """
    Parse SAM file with primers aligned to genome.
    Returns list of dicts with alignment details.
    
    SAM format fields:
    0: QNAME (primer name)
    1: FLAG
    2: RNAME (reference/chromosome)
    3: POS (position, 1-based)
    5: CIGAR
    9: SEQ
    """
    alignments = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            
            primer_name = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            sequence = fields[9]
            
            # Extract gene ID from primer name (format: primer_ENSGXXXXXXXX)
            gene_id = primer_name.replace('primer_', '')
            
            # Calculate end position from CIGAR
            match = re.match(r'(\d+)M', cigar)
            if match:
                length = int(match.group(1))
                end_pos = pos + length - 1  # Convert to inclusive end (0-based to 1-based SAM)
            else:
                end_pos = pos
            
            # Determine strand from FLAG
            # FLAG & 16 indicates reverse complement
            strand = '-' if (flag & 16) else '+'
            
            alignments.append({
                'primer_name': primer_name,
                'gene_id': gene_id,
                'chrom': chrom,
                'start': pos,
                'end': end_pos,
                'strand': strand,
                'sequence': sequence,
                'cigar': cigar,
                'flag': flag
            })
    
    return alignments

def check_overlap(query_start, query_end, region_start, region_end):
    """
    Check if query interval overlaps with region interval.
    Both intervals are 1-based inclusive.
    """
    return not (query_end < region_start or query_start > region_end)

def check_primers_in_regions(sam_file, fasta_file, output_file='primer_region_check.txt'):
    """
    Main function to check if primers fall within genomic regions.
    """
    regions = parse_fasta_regions(fasta_file)
    alignments = parse_sam_file(sam_file)
    
    results = {
        'inside': [],
        'outside': [],
        'no_region': []
    }
    
    for alignment in alignments:
        gene_id = alignment['gene_id']
        
        if gene_id not in regions:
            results['no_region'].append(alignment)
            continue
        
        region = regions[gene_id]
        
        # Check if primer is on the correct chromosome
        if alignment['chrom'] != region['chrom']:
            results['outside'].append({**alignment, 'region': region})
            continue
        
        # Check overlap
        if check_overlap(alignment['start'], alignment['end'], region['start'], region['end']):
            results['inside'].append({**alignment, 'region': region})
        else:
            results['outside'].append({**alignment, 'region': region})
    
    # Write results
    with open(output_file, 'w') as f:
        f.write("PRIMER-GENOMIC REGION OVERLAP CHECK\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"PRIMERS INSIDE REGIONS ({len(results['inside'])})\n")
        f.write("-" * 80 + "\n")
        for item in results['inside']:
            f.write(f"Primer: {item['primer_name']}\n")
            f.write(f"  Gene ID: {item['gene_id']}\n")
            f.write(f"  Genomic location: chr{item['chrom']}:{item['start']}-{item['end']} ({item['strand']})\n")
            f.write(f"  Region: chr{item['region']['chrom']}:{item['region']['start']}-{item['region']['end']} ({item['region']['strand']})\n")
            f.write(f"  Overlap: YES\n\n")
        
        f.write(f"\nPRIMERS OUTSIDE REGIONS ({len(results['outside'])})\n")
        f.write("-" * 80 + "\n")
        for item in results['outside']:
            f.write(f"Primer: {item['primer_name']}\n")
            f.write(f"  Gene ID: {item['gene_id']}\n")
            f.write(f"  Genomic location: chr{item['chrom']}:{item['start']}-{item['end']} ({item['strand']})\n")
            f.write(f"  Region: chr{item['region']['chrom']}:{item['region']['start']}-{item['region']['end']} ({item['region']['strand']})\n")
            f.write(f"  Overlap: NO\n\n")
        
        f.write(f"\nPRIMERS WITH NO MATCHING REGION ({len(results['no_region'])})\n")
        f.write("-" * 80 + "\n")
        for item in results['no_region']:
            f.write(f"Primer: {item['primer_name']}\n")
            f.write(f"  Gene ID: {item['gene_id']}\n")
            f.write(f"  Genomic location: chr{item['chrom']}:{item['start']}-{item['end']} ({item['strand']})\n")
            f.write(f"  Status: No genomic region defined for this gene\n\n")
    
    # Print summary
    total = len(alignments)
    inside = len(results['inside'])
    outside = len(results['outside'])
    no_region = len(results['no_region'])
    
    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print(f"{'='*80}")
    print(f"Total primers: {total}")
    print(f"Inside regions: {inside} ({100*inside/total:.1f}%)")
    print(f"Outside regions: {outside} ({100*outside/total:.1f}%)")
    print(f"No matching region: {no_region} ({100*no_region/total:.1f}%)")
    print(f"\nResults written to {output_file}")
    
    return results

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python primer_check.py <sam_file> <fasta_file> [output_file]")
        sys.exit(1)
    
    sam_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'primer_region_check.txt'
    
    check_primers_in_regions(sam_file, fasta_file, output_file)