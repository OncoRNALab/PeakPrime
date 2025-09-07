process MAKE_PRIMER3_INPUT {
  tag "$targets.baseName"
  publishDir params.outdir, mode: 'copy'
  conda "python>=3.10 primer3=2.6.*"

  input:
  path targets
  path p3_settings

  output:
  path "primer3_input.txt"

  script:
  """
  python ${projectDir}/bin/fasta_to_primer3.py --fasta ${targets} --settings ${p3_settings} --out primer3_input.txt
  """
}
