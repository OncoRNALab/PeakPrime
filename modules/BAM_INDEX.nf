process BAM_INDEX {
  tag "$bam.baseName"
  publishDir params.outdir, mode: 'copy', pattern: '*.bai'
  conda "samtools=1.19"

  input:
  path bam

  output:
  path "${bam}.bai"

  script:
  """
  samtools index -@ 4 ${bam}
  """
}
