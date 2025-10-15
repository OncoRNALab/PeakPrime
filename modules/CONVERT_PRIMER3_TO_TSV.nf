process CONVERT_PRIMER3_TO_TSV {
    tag "primer3_to_tsv"
    
    conda "${projectDir}/env/extract_3prime_env.yml"
    
    input:
    path primer3_output
    path gene_mapping
    
    output:
    path "primers.tsv", emit: primers_tsv
    
    script:
    def gene_map_arg = gene_mapping.name != 'NO_FILE' ? "--gene_mapping ${gene_mapping}" : ''
    """
    python ${projectDir}/bin/convert_primer3_to_tsv.py \\
        --primer3_output ${primer3_output} \\
        --output primers.tsv \\
        ${gene_map_arg}
    """
}
