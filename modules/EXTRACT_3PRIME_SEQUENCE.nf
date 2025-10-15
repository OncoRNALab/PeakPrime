process EXTRACT_3PRIME_SEQUENCE {
  tag "$fasta.baseName"
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/extract_3prime_env.yml"

  input:
  path fasta
  val length

  output:
  path '3prime_sequences.fasta'

  script:
  """
  python ${projectDir}/bin/extract_3prime_sequence.py \
    --input ${fasta} \
    --output 3prime_sequences.fasta \
    --length ${length}
  """
}
