process OPTIMIZE_PRIMERS_MULTIPEAK {
    tag "Multi-peak optimization"
    publishDir "${params.outdir}/optimized_primers", mode: 'copy'
    
    conda "${projectDir}/env/optimize_primers_multipeak_env.yml"
    
    input:
    path best_primers
    path alignment_summary
    
    output:
    path "optimized_primers_multipeak.tsv", emit: optimized_primers
    path "optimization_stats.txt", emit: stats
    
    script:
    def alignment_arg = alignment_summary.name != 'NO_FILE' ? "--alignment_summary ${alignment_summary}" : ""
    def mismatch_arg = params.max_mismatches ? "--max_mismatches ${params.max_mismatches}" : ""
    """
    python ${projectDir}/bin/optimize_primers_multipeak.py \\
        --best_primers ${best_primers} \\
        ${alignment_arg} \\
        --out optimized_primers_multipeak.tsv \\
        --primers_per_gene ${params.primers_per_gene} \\
        --distance_weight ${params.distance_weight} \\
        --isoform_weight ${params.isoform_weight} \\
        --peak_rank_weight ${params.peak_rank_weight} \\
        ${mismatch_arg} \\
        2>&1 | tee optimization_stats.txt
    
    # Add summary info
    echo "" >> optimization_stats.txt
    echo "=== PARAMETER SUMMARY ===" >> optimization_stats.txt
    echo "Primers per gene: ${params.primers_per_gene}" >> optimization_stats.txt
    echo "Distance weight: ${params.distance_weight}" >> optimization_stats.txt
    echo "Isoform weight: ${params.isoform_weight}" >> optimization_stats.txt
    echo "Peak rank weight: ${params.peak_rank_weight}" >> optimization_stats.txt
    echo "" >> optimization_stats.txt
    echo "Total optimized primers: \$(tail -n +2 optimized_primers_multipeak.tsv | wc -l)" >> optimization_stats.txt
    """
}
