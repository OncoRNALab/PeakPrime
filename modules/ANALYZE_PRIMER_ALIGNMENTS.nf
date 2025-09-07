process ANALYZE_PRIMER_ALIGNMENTS {
  tag 'alignment_analysis'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.9 biopython=1.79 pandas=1.5.* pysam=0.21.*"

  input:
  path alignment_bam
  path alignment_bai
  path primers_tsv
  path transcriptome_fasta

  output:
  path 'primer_alignment_report.tsv'
  path 'primer_alignment_summary.tsv'

  script:
  def transcriptome_arg = transcriptome_fasta.name != 'NO_FILE' ? "--transcriptome_fasta ${transcriptome_fasta}" : ""
  """
  python ${projectDir}/bin/analyze_primer_alignments.py \
    --alignment_bam ${alignment_bam} \
    --primers_tsv ${primers_tsv} \
    ${transcriptome_arg} \
    --out_report primer_alignment_report.tsv \
    --out_summary primer_alignment_summary.tsv
  """
}
