process PRIMERS_TO_FASTA {
  tag 'primers_to_fasta'
  publishDir params.outdir, mode: 'copy'
  // Use environment file to avoid r-essentials metapackage (100+ deps) which causes slow/hanging conda resolution
  conda "${projectDir}/env/primers_to_fasta_env.yml"
  // Alternative inline (if env file doesn't work): conda "r-base=4.3.* bioconductor-biostrings r-optparse"

  input:
  path primers_tsv

  output:
  path 'primers_for_alignment.fa'

  script:
  def max_primers_arg = params.max_primers_per_gene ? "--max_primers_per_gene ${params.max_primers_per_gene}" : ""
  """
  Rscript ${projectDir}/bin/primers_to_fasta.R \
    --primers_tsv ${primers_tsv} \
    --out_fasta primers_for_alignment.fa \
    ${max_primers_arg}
  """
}
