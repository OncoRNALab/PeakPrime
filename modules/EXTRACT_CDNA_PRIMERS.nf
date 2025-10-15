process EXTRACT_CDNA_PRIMERS {
  tag 'extract_cdna_primers'
  publishDir params.outdir, mode: 'copy'
  // Use environment file to avoid r-essentials metapackage (100+ deps)
  conda "${projectDir}/env/extract_cdna_primers_env.yml"
  // Alternative: conda "r-base=4.3.* r-optparse r-data.table"

  input:
  path primer3_output
  path peaks_tsv

  output:
  path 'cdna_primers.tsv'

  script:
  """
  Rscript ${projectDir}/bin/extract_cdna_primers.R \
    --primer3_output ${primer3_output} \
    --peaks_tsv ${peaks_tsv} \
    --out_tsv cdna_primers.tsv
  """
}
