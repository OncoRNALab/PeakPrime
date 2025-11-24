process BAM_INDEX {
  tag "$bam.baseName"
  publishDir params.outdir, mode: 'copy', pattern: '*.bai'
  conda "${projectDir}/env/bam_index_env.yml"

  input:
  path bam

  output:
  path "${bam}.bai"

  script:
  """
  # Verify samtools is available
  if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH"
    echo "Conda environment may not be activated properly"
    echo "PATH: \$PATH"
    exit 1
  fi
  
  samtools index -@ 4 ${bam}
  """
}
