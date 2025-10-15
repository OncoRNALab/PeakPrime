process SIMPLE_PRIMERS_TO_FASTA {
  tag 'simple_primers_to_fasta'
  publishDir params.outdir, mode: 'copy'
  conda "python>=3.10"

  input:
  path primer3_output

  output:
  path 'primers_for_alignment.fa'

  script:
  """
  python ${projectDir}/bin/simple_primers_to_fasta.py \
    --primer3-output ${primer3_output} \
    --output-fasta primers_for_alignment.fa
  """
}
