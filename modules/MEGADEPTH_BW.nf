process MEGADEPTH_BW {
  tag "$bam.baseName"
  publishDir params.outdir, mode: 'copy', pattern: '*.bw'
  conda "${projectDir}/env/megadepth_env.yml"

  input:
  path bam
  path bai

  output:
  path "${bam.baseName}.bw"

  script:
  """
  megadepth ${bam} --bigwig --prefix ${bam.baseName}
  # List files to see what was actually created
  ls -la *.bw
  # Find the correct output file and rename it
  if [ -f "${bam.baseName}.all.bw" ]; then
    mv ${bam.baseName}.all.bw ${bam.baseName}.bw
  elif [ -f "${bam.baseName}.coverage.bw" ]; then
    mv ${bam.baseName}.coverage.bw ${bam.baseName}.bw
  elif [ -f "coverage.bw" ]; then
    mv coverage.bw ${bam.baseName}.bw
  else
    echo "Available files:"
    ls -la
    exit 1
  fi
  """
}
