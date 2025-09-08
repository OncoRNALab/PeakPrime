#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyze primer alignment results from bowtie2 - Python replacement for R script
This script processes BAM files and produces comprehensive alignment analysis
including all metrics from make_alignment_summary.py
"""

import argparse
import pandas as pd
import pysam
from Bio import SeqIO
import sys
import os

def parse_gene_names(transcriptome_fasta=None, transcript_mapping=None):
    """
    Parse gene names from transcriptome FASTA headers or load from mapping file
    Returns dictionary mapping transcript_id -> gene_name
    """
    
    if transcript_mapping and os.path.exists(transcript_mapping):
        print("Loading pre-built transcript-to-gene mapping...")
        try:
            mapping_data = pd.read_csv(transcript_mapping, sep='\t')
            
            if 'transcript_id' not in mapping_data.columns or 'gene_name' not in mapping_data.columns:
                raise ValueError("Mapping file must contain 'transcript_id' and 'gene_name' columns")
            
            transcript_to_gene = dict(zip(mapping_data['transcript_id'], mapping_data['gene_name']))
            print(f"Loaded mapping for {len(transcript_to_gene)} transcripts")
            return transcript_to_gene
            
        except Exception as e:
            print(f"Error loading mapping file: {e}")
            print("Falling back to FASTA parsing...")
    
    # Fallback to FASTA parsing if no mapping file
    if not transcriptome_fasta or not os.path.exists(transcriptome_fasta):
        print("Warning: Transcriptome FASTA not provided or not found. Gene names will not be available.")
        return {}
    
    print("Parsing gene names from transcriptome FASTA (this may take time)...")
    transcript_to_gene = {}
    
    try:
        with open(transcriptome_fasta) as fasta_file:
            for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
                # Extract transcript ID (first part before space)
                transcript_id = record.id.split()[0] if ' ' in record.id else record.id
                
                # Extract gene name from header (assuming format like "ENST00000456328 gene=DDX11L2")
                gene_name = "Unknown"
                if "gene=" in record.description:
                    try:
                        gene_part = record.description.split("gene=")[1].split()[0]
                        gene_name = gene_part
                    except IndexError:
                        pass
                
                transcript_to_gene[transcript_id] = gene_name
                
                # Progress indicator for large files
                if (i + 1) % 10000 == 0:
                    print(f"Processed {i + 1} sequences")
        
        print(f"Parsed {len(transcript_to_gene)} transcript-to-gene mappings")
        
        # Debug: show a few example mappings
        test_transcripts = ["ENST00000521262", "ENST00000524349", "ENST00000371437"]
        found_count = 0
        for test_id in test_transcripts:
            if test_id in transcript_to_gene:
                found_count += 1
                if found_count <= 2:
                    print(f"  Example mapping: {test_id} -> {transcript_to_gene[test_id]}")
        
        return transcript_to_gene
        
    except Exception as e:
        print(f"Error parsing transcriptome FASTA: {e}")
        return {}

def analyze_alignments(alignment_bam, primers_tsv, transcriptome_fasta=None, 
                      transcript_mapping=None, out_report="primer_alignment_report.tsv", 
                      out_summary="primer_alignment_summary.tsv"):
    """
    Main function to analyze primer alignments
    """
    
    # Read original primer data
    print("Loading primer data...")
    primers = pd.read_csv(primers_tsv, sep='\t')
    primers['primer_id'] = (primers['gene_id'].astype(str) + "|" +
                           "idx" + primers['primer_index'].astype(str) + "|" +
                           primers['primer_type'].astype(str) + "|" +
                           "strand" + primers['gene_strand'].astype(str))
    
    # Create a lookup dictionary for gene strand information
    primer_to_gene_strand = dict(zip(primers['primer_id'], primers['gene_strand']))
    
    print(f"Loaded {len(primers)} primers")
    
    # Parse transcriptome FASTA for gene name mapping
    transcript_to_gene = parse_gene_names(transcriptome_fasta, transcript_mapping)
    
    # Read alignment results
    print("Reading BAM file...")
    try:
        samfile = pysam.AlignmentFile(alignment_bam, "r")
        
        # Get reference sequence lengths for distance calculations
        ref_lengths = {sq['SN']: sq['LN'] for sq in samfile.header['SQ']}
        
        alignments_data = []
        
        for read in samfile:
            if read.is_unmapped:
                continue
            
            # Only process forward-aligned primers (5' to 3' orientation)
            if read.is_reverse:
                continue
                
            # Extract alignment information including all metrics from make_alignment_summary.py
            gene_id = read.query_name.split("_")[0]
            primer_id_parts = "_".join(read.query_name.split("_")[2:])
            transcript_id = read.reference_name
            transcript_length = ref_lengths.get(transcript_id, None)
            
            # Get gene name from transcript mapping
            gene_name = transcript_to_gene.get(transcript_id, "Unknown")
            
            # Alignment metrics
            alignment_length = read.query_alignment_length if read.query_alignment_length else 0
            start_pos = read.reference_start
            end_pos = read.reference_end if read.reference_end is not None else start_pos + alignment_length
            
            # Distance to end of transcript - strand-specific calculation
            # This represents the total distance from the primer start to the transcript end
            # (includes the primer length itself)
            distance_to_end = "NA"
            if transcript_length is not None:
                # Get gene strand information for this primer
                gene_strand = primer_to_gene_strand.get(read.query_name, "+")  # default to positive
                
                if gene_strand == "-":
                    # For negative strand genes: primers align in reverse orientation
                    # Distance from primer start (leftmost position = 3' end of primer) to transcript start (5' end of gene)
                    # This includes the primer length
                    distance_to_end = start_pos
                else:
                    # For positive strand genes: 
                    # Distance from primer start (leftmost position = 5' end of primer) to transcript end (3' end of gene)
                    # This includes the primer length
                    distance_to_end = transcript_length - start_pos
            
            # Alignment direction - should always be forward after filtering
            direction = "forward"  # Only forward alignments are processed
            
            # Extract mismatches from NM tag
            mismatches = read.get_tag("NM") if read.has_tag("NM") else "NA"
            
            # MAPQ score
            mapq = read.mapping_quality if read.mapping_quality is not None else 255
            
            alignments_data.append({
                'primer_id': read.query_name,
                'gene_id': gene_id,
                'primer_index_parts': primer_id_parts,
                'transcript_id': transcript_id,
                'gene_name': gene_name,
                'gene_strand': primer_to_gene_strand.get(read.query_name, "+"),
                'alignment_start': start_pos,
                'alignment_end': end_pos,
                'alignment_length': alignment_length,
                'mismatches': mismatches,
                'direction': direction,
                'distance_to_end': distance_to_end,
                'mapq': mapq,
                'transcript_length': transcript_length if transcript_length else "NA"
            })
        
        samfile.close()
        
        print(f"Read {len(alignments_data)} alignments from BAM file")
        
    except Exception as e:
        print(f"Error reading BAM file: {e}")
        alignments_data = []
    
    # Convert to DataFrame
    if alignments_data:
        alignment_df = pd.DataFrame(alignments_data)
    else:
        # Create empty DataFrame with correct columns
        alignment_df = pd.DataFrame(columns=['primer_id', 'gene_id', 'primer_index_parts', 
                                           'transcript_id', 'gene_name', 'gene_strand', 'alignment_start', 
                                           'alignment_end', 'alignment_length', 'mismatches', 
                                           'direction', 'distance_to_end', 'mapq', 'transcript_length'])
    
    # Count alignments per primer
    if not alignment_df.empty:
        alignment_counts = alignment_df['primer_id'].value_counts().to_dict()
    else:
        alignment_counts = {}
    
    # Create comprehensive report
    print("Creating alignment report...")
    report = primers.copy()
    report['num_alignments'] = report['primer_id'].map(alignment_counts).fillna(0).astype(int)
    
    # Write report
    report.to_csv(out_report, sep='\t', index=False)
    print(f"Alignment report written to: {out_report}")
    
    # Create detailed alignment summary (all alignments)
    print("Creating detailed alignment summary...")
    
    if not alignment_df.empty:
        # Use all alignments since MAPQ is meaningless with bowtie2 -a
        print(f"Processing {len(alignment_df)} total alignments")
        
        # Merge with primer information
        primer_cols = ['primer_id', 'gene_id', 'primer_index', 'primer_type', 'primer_sequence', 'gene_strand']
        detailed_summary = alignment_df.merge(
            primers[primer_cols], 
            on='primer_id', 
            how='left',
            suffixes=('_align', '_primer')
        )
        
        # Use primer gene_id and gene_strand if alignment versions are missing
        detailed_summary['gene_id'] = detailed_summary['gene_id_primer'].fillna(detailed_summary['gene_id_align'])
        detailed_summary['gene_strand'] = detailed_summary['gene_strand_primer'].fillna(detailed_summary['gene_strand_align'])
        
        # Select and rename columns (remove MAPQ references)
        detailed_summary = detailed_summary.rename(columns={
            'transcript_id': 'aligned_transcript',
            'gene_name': 'aligned_gene_name',
            'direction': 'alignment_strand'
        })
        
        # Reorder columns (remove alignment_mapq)
        output_columns = [
            'gene_id', 'primer_index', 'primer_type', 'primer_sequence', 'gene_strand',
            'aligned_transcript', 'aligned_gene_name', 'alignment_start', 'alignment_end',
            'alignment_length', 'alignment_strand', 'mismatches', 
            'distance_to_end', 'transcript_length'
        ]
        
        # Debug: check which columns are available
        available_columns = detailed_summary.columns.tolist()
        missing_columns = [col for col in output_columns if col not in available_columns]
        if missing_columns:
            print(f"Warning: Missing columns: {missing_columns}")
            print(f"Available columns: {available_columns}")
            # Use only available columns
            output_columns = [col for col in output_columns if col in available_columns]
        
        detailed_summary = detailed_summary[output_columns]
        
        # Sort by gene, primer index, and number of mismatches (best alignments first)
        detailed_summary = detailed_summary.sort_values(['gene_id', 'primer_index', 'mismatches'], 
                                                      ascending=[True, True, True])
        
    else:
        detailed_summary = pd.DataFrame(columns=[
            'gene_id', 'primer_index', 'primer_type', 'primer_sequence', 'gene_strand',
            'aligned_transcript', 'aligned_gene_name', 'alignment_start', 'alignment_end',
            'alignment_length', 'alignment_strand', 'mismatches',
            'distance_to_end', 'transcript_length'
        ])
    
    # Write detailed summary
    detailed_summary.to_csv(out_summary, sep='\t', index=False)
    print(f"Detailed alignment summary written to: {out_summary}")
    
    # Print summary statistics
    print("\n=== PRIMER ALIGNMENT SUMMARY ===")
    print(f"Total primers analyzed: {len(report)}")
    
    print("\nAlignment count distribution:")
    alignment_bins = pd.cut(report['num_alignments'], 
                           bins=[-1, 0, 1, 5, 20, float('inf')], 
                           labels=['No hits', '1 hit', '2-5 hits', '6-20 hits', '>20 hits'])
    alignment_summary = alignment_bins.value_counts()
    for bin_name in alignment_summary.index:
        print(f"  {bin_name}: {alignment_summary[bin_name]}")
    
    # Identify primers with many alignments (potential cross-reactivity)
    many_alignments = report[report['num_alignments'] > 20]
    if not many_alignments.empty:
        print("\n=== PRIMERS WITH MANY ALIGNMENTS ===")
        print("Primers with >20 hits (potential cross-reactivity):")
        for i in range(min(10, len(many_alignments))):
            p = many_alignments.iloc[i]
            print(f"  {p['gene_id']} primer {p['primer_index']} ({p['primer_type']}) - {p['num_alignments']} hits")
        if len(many_alignments) > 10:
            print(f"  ... and {len(many_alignments) - 10} more (see full report)")
    
    # Identify primers with unique alignments (ideal case)
    unique_primers = report[report['num_alignments'] == 1]
    if not unique_primers.empty:
        print("\n=== PRIMERS WITH UNIQUE ALIGNMENTS ===")
        print("Primers with exactly 1 transcriptome hit:")
        for i in range(min(5, len(unique_primers))):
            p = unique_primers.iloc[i]
            print(f"  {p['gene_id']} primer {p['primer_index']} ({p['primer_type']})")
        if len(unique_primers) > 5:
            print(f"  ... and {len(unique_primers) - 5} more unique primers")
    
    print("\nAlignment analysis complete")

def main():
    parser = argparse.ArgumentParser(description="Analyze primer alignment results from bowtie2")
    parser.add_argument('--alignment_bam', required=True, help='Input BAM file from bowtie2 alignment')
    parser.add_argument('--primers_tsv', required=True, help='Original primers TSV file')
    parser.add_argument('--transcriptome_fasta', help='Transcriptome FASTA file for gene name mapping')
    parser.add_argument('--transcript_mapping', help='Pre-built transcript-to-gene mapping TSV file (faster alternative to FASTA)')
    parser.add_argument('--out_report', default='primer_alignment_report.tsv', help='Output alignment report')
    parser.add_argument('--out_summary', default='primer_alignment_summary.tsv', help='Detailed alignment summary per primer')
    
    args = parser.parse_args()
    
    if not args.alignment_bam or not args.primers_tsv:
        print("Error: Both --alignment_bam and --primers_tsv are required")
        sys.exit(1)
    
    analyze_alignments(
        alignment_bam=args.alignment_bam,
        primers_tsv=args.primers_tsv,
        transcriptome_fasta=args.transcriptome_fasta,
        transcript_mapping=args.transcript_mapping,
        out_report=args.out_report,
        out_summary=args.out_summary
    )

if __name__ == "__main__":
    main()
