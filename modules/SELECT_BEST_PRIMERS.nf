process SELECT_BEST_PRIMERS {
  tag 'select_best_primers'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.9 pandas=1.5.*"

  input:
  path report
  path summary

  output:
  path 'best_primers.tsv'

  script:
  """
  python ${projectDir}/bin/select_best_primer.py \
    --report ${report} \
    --summary ${summary} \
    --distance_threshold ${params.distance_threshold} \
    --out best_primers.tsv
  """
}
