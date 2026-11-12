process SUMMARIZE_RESULTS {
  tag 'pipeline_summary'
  publishDir "${params.outdir}/summary", mode: 'copy'
  conda "${projectDir}/env/summarize_results_env.yml"

  input:
  path genes
  path qc_summary
  path peaks_tsv
  path cdna_primers
  path best_primers
  path optimized_primers
  path alignment_summary

  output:
  path 'run_summary.txt', emit: text
  path 'gene_funnel.tsv', emit: gene_funnel
  path 'primer_specificity.tsv', emit: specificity
  path 'final_specific_primers.tsv', emit: final_primers

  script:
  def mm = params.max_mismatches != null ? params.max_mismatches : 3
  def dist = params.distance_threshold != null ? params.distance_threshold : 1000
  """
  # Normalize optional missing inputs to NO_FILE for the summary script
  qc_arg=${qc_summary}
  peaks_arg=${peaks_tsv}
  cdna_arg=${cdna_primers}
  best_arg=${best_primers}
  opt_arg=${optimized_primers}
  aln_arg=${alignment_summary}

  [ -s "\$qc_arg" ] || qc_arg=NO_FILE
  [ -s "\$peaks_arg" ] || peaks_arg=NO_FILE
  [ -s "\$cdna_arg" ] || cdna_arg=NO_FILE
  [ -s "\$best_arg" ] || best_arg=NO_FILE
  [ -s "\$opt_arg" ] || opt_arg=NO_FILE
  [ -s "\$aln_arg" ] || aln_arg=NO_FILE

  python ${projectDir}/bin/summarize_primer_run.py \\
    --genes ${genes} \\
    --qc \$qc_arg \\
    --peaks \$peaks_arg \\
    --cdna_primers \$cdna_arg \\
    --best_primers \$best_arg \\
    --optimized_primers \$opt_arg \\
    --alignment_summary \$aln_arg \\
    --max_mismatches ${mm} \\
    --distance_threshold ${dist} \\
    --outdir .
  """
}
