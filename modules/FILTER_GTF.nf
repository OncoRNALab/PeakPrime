process FILTER_GTF {
  tag 'filter_gtf'
  conda "${projectDir}/env/filter_gtf_env.yml"
  
  input:
  path gtf
  path genes
  
  output:
  path 'filtered.gtf'
  
  script:
  """
  Rscript ${projectDir}/bin/filter_gtf_by_genes.R \\
    --gtf ${gtf} \\
    --genes ${genes} \\
    --out filtered.gtf
  """
}
