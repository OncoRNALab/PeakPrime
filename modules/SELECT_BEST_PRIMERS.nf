process SELECT_BEST_PRIMERS {
  tag 'select_best_primers'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/select_best_primers_env.yml"

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
    --max_mismatches ${params.max_mismatches} \
    --out best_primers.tsv
  """
}
