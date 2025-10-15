/*
========================================================================================
    MODULE: EXTRACT_GENE_SEQUENCES
========================================================================================
    Extract gene sequences from genome using annotation
----------------------------------------------------------------------------------------
*/

process EXTRACT_GENE_SEQUENCES {
    tag "$gene_id"
    label 'process_medium'
    conda "${moduleDir}/../envs/primer_design.yml"
    
    publishDir "${params.outdir}/sequences", mode: 'copy'
    
    input:
    val gene_id
    path genome_fasta
    path annotation_gtf
    
    output:
    tuple val(gene_id), path("${gene_id}_sequence.fasta"), emit: sequences
    path "${gene_id}_info.json", emit: info
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import sys
    from collections import defaultdict
    
    def parse_gtf_attribute(attr_str):
        """Parse GTF attribute string"""
        attrs = {}
        for item in attr_str.strip().split(';'):
            item = item.strip()
            if item:
                key_val = item.split(' ', 1)
                if len(key_val) == 2:
                    key, val = key_val
                    attrs[key] = val.strip('"')
        return attrs
    
    def extract_gene_info(gene_id, gtf_file):
        """Extract gene information from GTF"""
        gene_info = {
            'gene_id': gene_id,
            'chromosome': None,
            'strand': None,
            'start': None,
            'end': None,
            'exons': []
        }
        
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\\t')
                if len(fields) < 9:
                    continue
                
                chrom, source, feature, start, end, score, strand, frame, attributes = fields
                attrs = parse_gtf_attribute(attributes)
                
                # Check if this is our gene
                if attrs.get('gene_id') == gene_id:
                    if feature == 'gene':
                        gene_info['chromosome'] = chrom
                        gene_info['strand'] = strand
                        gene_info['start'] = int(start)
                        gene_info['end'] = int(end)
                    elif feature == 'exon':
                        gene_info['exons'].append({
                            'start': int(start),
                            'end': int(end)
                        })
        
        # Sort exons by position
        gene_info['exons'].sort(key=lambda x: x['start'])
        
        return gene_info
    
    def extract_sequence(chrom, start, end, fasta_file):
        """Extract sequence from FASTA file"""
        from Bio import SeqIO
        
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id == chrom:
                # Convert to 0-based indexing
                seq = str(record.seq[start-1:end])
                return seq
        return None
    
    def reverse_complement(seq):
        """Return reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join([complement.get(base, 'N') for base in seq[::-1]])
    
    # Main execution
    gene_id = "${gene_id}"
    gtf_file = "${annotation_gtf}"
    fasta_file = "${genome_fasta}"
    
    # Extract gene information
    gene_info = extract_gene_info(gene_id, gtf_file)
    
    if gene_info['chromosome'] is None:
        print(f"Error: Gene {gene_id} not found in GTF file", file=sys.stderr)
        sys.exit(1)
    
    # Extract sequence
    sequence = extract_sequence(
        gene_info['chromosome'],
        gene_info['start'],
        gene_info['end'],
        fasta_file
    )
    
    if sequence is None:
        print(f"Error: Could not extract sequence for {gene_id}", file=sys.stderr)
        sys.exit(1)
    
    # Reverse complement if on minus strand
    if gene_info['strand'] == '-':
        sequence = reverse_complement(sequence)
    
    # Write sequence to file
    with open(f"{gene_id}_sequence.fasta", 'w') as f:
        f.write(f">{gene_id}\\n")
        # Write sequence in 60bp lines
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\\n')
    
    # Write gene info
    gene_info['sequence_length'] = len(sequence)
    with open(f"{gene_id}_info.json", 'w') as f:
        json.dump(gene_info, f, indent=2)
    
    print(f"Extracted {len(sequence)} bp for gene {gene_id}")
    """
}