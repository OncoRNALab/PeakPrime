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
  samtools index -@ 4 ${bam}
  """
}
