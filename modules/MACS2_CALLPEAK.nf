process MACS2_CALLPEAK {
    tag "$sample_id"
    publishDir "${params.outdir}/macs2_peaks", mode: 'copy'
    
    // Use environment file for better reproducibility across different conda versions
    conda "${projectDir}/env/macs2_env.yml"
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: narrowpeak
    path("${sample_id}_*"), emit: all_outputs
    
    script:
    // Build MACS2 command - use simple default or add fragment modeling if specified
    def use_nomodel = params.macs2_extsize ? true : false
    def fragment_args = use_nomodel ? "--nomodel --extsize ${params.macs2_extsize}" : ""
    def shift_arg = (use_nomodel && params.macs2_shift) ? "--shift ${params.macs2_shift}" : ""
    def summits_arg = use_nomodel ? "--call-summits" : ""
    def qvalue_threshold = params.macs2_qvalue_threshold ?: params.macs2_pvalue_threshold
    def qvalue_arg = qvalue_threshold ? "-q ${qvalue_threshold}" : ""
    
    """
    # First command: Call peaks with MACS2
    macs2 callpeak \\
      -t ${bam_file} \\
      --outdir . \\
      -n ${sample_id} \\
      -g hs \\
      --bdg \\
      --keep-dup auto \\
      ${qvalue_arg} \\
      ${fragment_args} \\
      ${shift_arg} \\
      ${summits_arg}
    
    # Check if the narrowPeak file was created
    if [ ! -f "${sample_id}_peaks.narrowPeak" ]; then
        echo "Error: MACS2 did not produce expected narrowPeak file"
        exit 1
    fi
    """
}
