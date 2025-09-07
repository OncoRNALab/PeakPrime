process RUN_PRIMER3 {
  tag 'primer3'
  publishDir params.outdir, mode: 'copy'
  conda "python>=3.10 primer3=2.6.*"

  input:
  path p3in

  output:
  path 'primer3_output.txt'

  script:
  """
  primer3_core < ${p3in} > primer3_output.txt
  """
}
