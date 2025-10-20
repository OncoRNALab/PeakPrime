process OPTIMIZE_PRIMER_ISOFORMS {
  tag 'optimize_isoforms'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.10 pandas=2.0"

  input:
  path best_primers_file
  path alignment_summary_file

  output:
  path 'best_primers_optimal.tsv', emit: optimized_primers
  
  script:
  def distance_arg = params.distance_threshold ? "--distance_threshold ${params.distance_threshold}" : ""
  """
  optimize_primer_isoforms.py \\
    --best_primers ${best_primers_file} \\
    --alignment_summary ${alignment_summary_file} \\
    --output best_primers_optimal.tsv \\
    ${distance_arg}
  """
}
