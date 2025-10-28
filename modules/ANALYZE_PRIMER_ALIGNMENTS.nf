process ANALYZE_PRIMER_ALIGNMENTS {
  tag 'alignment_analysis'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/analyze_primer_alignments_env.yml"

  input:
  path alignment_bam
  path alignment_bai
  path primers_tsv
  path transcript_mapping

  output:
  path 'primer_alignment_report.tsv'
  path 'primer_alignment_summary.tsv'

  script:
  def mapping_arg = transcript_mapping.name != 'NO_FILE' ? "--transcript_mapping ${transcript_mapping}" : ""
  """
  python ${projectDir}/bin/analyze_primer_alignments.py \
    --alignment_bam ${alignment_bam} \
    --primers_tsv ${primers_tsv} \
    ${mapping_arg} \
    --out_report primer_alignment_report.tsv \
    --out_summary primer_alignment_summary.tsv
  """
}
