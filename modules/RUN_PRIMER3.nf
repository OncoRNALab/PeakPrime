process RUN_PRIMER3 {
  tag 'primer3'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/primer3_env.yml"

  input:
  path p3in

  output:
  path 'primer3_output.txt'

  script:
  def mispriming_lib = file("${projectDir}/resources/humrep_and_simple.txt")
  def stage_lib = mispriming_lib.exists() ? "ln -s ${mispriming_lib} humrep_and_simple.txt" : "echo 'Warning: mispriming library not found'"
  """
  ${stage_lib}
  primer3_core < ${p3in} > primer3_output.txt
  """
}
