process PROCESS_MACS2_PEAKS {
    tag "process_peaks"
    publishDir "${params.outdir}/processed_peaks", mode: 'copy'
    
    // Use environment file for cleaner code and better maintainability
    conda "${projectDir}/env/process_macs2_peaks_env.yml"
    // Alternative inline: conda 'bioconda::bioconductor-genomicfeatures bioconda::bioconductor-rtracklayer bioconda::bioconductor-genomicranges bioconda::bioconductor-iranges bioconda::bioconductor-s4vectors bioconda::bioconductor-biostrings bioconda::bioconductor-txdbmaker conda-forge::r-optparse conda-forge::r-data.table'
    
    input:
    tuple val(sample_id), path(narrowpeak_file)
    path gtf_file
    path genes_file
    path fasta_file
    
    output:
    path "selected_peaks.fa", emit: fasta
    path "selected_peaks.bed", emit: bed  
    path "selected_peaks.tsv", emit: tsv
    path "peaks_qc_summary.tsv", emit: qc
    
    script:
    """
    process_macs2_peaks.R \\
        --narrowpeak ${narrowpeak_file} \\
        --gtf ${gtf_file} \\
        --genes ${genes_file} \\
        --fasta ${fasta_file} \\
        --qvalue_threshold ${params.macs2_qvalue_threshold ?: (params.macs2_pvalue_threshold ?: 0.05)} \\
        --min_peak_score ${params.macs2_min_peak_score ?: 0} \\
        --peak_selection_metric ${params.peak_selection_metric ?: 'score'} \\
        --peak_rank ${params.peak_rank ?: 1} \\
        --min_exonic_fraction ${params.min_exonic_fraction ?: 'NA'} \\
        --trim_to_exon ${params.trim_to_exon ? 'true' : 'false'} \\
        --force_exonic_trimming ${params.force_exonic_trimming ? 'true' : 'false'} \\
        --min_trimmed_length ${params.min_trimmed_length ?: 150} \\
        --out_fa selected_peaks.fa \\
        --out_bed selected_peaks.bed \\
        --out_peaks selected_peaks.tsv \\
        --out_qc peaks_qc_summary.tsv
    """
}
