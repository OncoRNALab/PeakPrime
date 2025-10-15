process FETCH_MANE_TRANSCRIPTS {
  tag 'fetch_mane'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/fetch_mane_env.yml"

  input:
  path gene_ids_file

  output:
  path 'mane_transcripts.fasta'
  path 'mane_mapping.tsv'

  script:
  """
  ${projectDir}/bin/fetch_mane_transcripts.R \
    --gene-ids ${gene_ids_file} \
    --output-fasta mane_transcripts.fasta \
    --output-mapping mane_mapping.tsv \
    --allow-fallback
  """
}
