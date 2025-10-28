process SUMMARIZE_RESULTS {
  tag 'pipeline_summary'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/summarize_results_env.yml"

  input:
  path best_primers_file, stageAs: 'best_primers.tsv'
  path cdna_primers_file, stageAs: 'cdna_primers.tsv'
  path qc_summary_file, stageAs: 'peaks_qc_summary.tsv'
  path alignment_summary_file, stageAs: 'primer_alignment_summary.tsv'

  output:
  path 'pipeline_summary.txt'

  script:
  """
  # Create the expected directory structure in the working directory
  mkdir -p processed_peaks
  
  # Only copy if the file is not already in the target location
  if [ "peaks_qc_summary.tsv" -nt "processed_peaks/peaks_qc_summary.tsv" ] || [ ! -f "processed_peaks/peaks_qc_summary.tsv" ]; then
    cp peaks_qc_summary.tsv processed_peaks/
  fi
  
  # The other files are already staged with the correct names in the current directory
  # Run the summary script pointing to current directory
  python ${projectDir}/bin/summarize_pipeline_results.py . > pipeline_summary.txt
  """
}
