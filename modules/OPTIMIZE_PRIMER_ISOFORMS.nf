process OPTIMIZE_PRIMER_ISOFORMS {
  tag 'optimize_isoforms'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/optimize_primer_isoforms_env.yml"

  input:
  path best_primers_file
  path alignment_summary_file

  output:
  path 'best_primers_optimal.tsv', emit: optimized_primers
  
  script:
  def distance_arg = params.distance_threshold ? "--distance_threshold ${params.distance_threshold}" : ""
  def mismatch_arg = params.max_mismatches ? "--max_mismatches ${params.max_mismatches}" : ""
  """
  optimize_primer_isoforms.py \\
    --best_primers ${best_primers_file} \\
    --alignment_summary ${alignment_summary_file} \\
    --output best_primers_optimal.tsv \\
    ${distance_arg} \\
    ${mismatch_arg}
  """
}
