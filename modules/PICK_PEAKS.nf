process PICK_PEAKS {
  tag "$bw.baseName"
  publishDir params.outdir, mode: 'copy'
  conda "r-base>=4.3 r-essentials bioconductor-derfinder bioconductor-genomicfeatures bioconductor-genomicranges bioconductor-rtracklayer bioconductor-iranges bioconductor-s4vectors bioconductor-bsgenome bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-rsamtools bioconductor-txdbmaker r-optparse"

  input:
  path bw
  path gtf
  path genes
  path fasta

  output:
  path "primer_targets.fa"
  path "primer_targets.bed"
  path "peaks.tsv"
  path "qc_coverage_summary.tsv"

  script:
  def fasta_arg = fasta.name != 'NO_FILE' ? "--fasta ${fasta}" : ""
  def sw_flag = params.sliding_window ? "--sliding_window" : ""
  def minmean_arg = params.min_window_mean != null ? "--min_window_mean ${params.min_window_mean}" : ""
  def minmean_pct_arg = params.min_window_mean_pct != null ? "--min_window_mean_pct ${params.min_window_mean_pct}" : ""
  def maxgap_arg = params.max_gap != null ? "--max_gap ${params.max_gap}" : ""
  def slop_arg = params.search_slop != null ? "--search_slop ${params.search_slop}" : ""
  def trim_arg = params.trim_to_exon ? "--trim_to_exon" : ""
  def trim_cov_arg = params.trim_low_coverage_pct != null ? "--trim_low_coverage_pct ${params.trim_low_coverage_pct}" : ""
  def min_exonic_arg = params.min_exonic_fraction != null ? "--min_exonic_fraction ${params.min_exonic_fraction}" : ""
  """
  Rscript ${projectDir}/bin/pick_peaks_span.R \
    --bw ${bw} \
    --gtf ${gtf} \
    --genes ${genes} \
    --genome ${params.genome_package} \
    ${fasta_arg} \
    --pad ${params.pad} \
    --smooth_k ${params.smooth_k} \
    ${sw_flag} \
    ${minmean_arg} \
    ${minmean_pct_arg} \
    ${maxgap_arg} \
    ${slop_arg} \
    ${trim_arg} \
    ${trim_cov_arg} \
    ${min_exonic_arg} \
    --out_fa primer_targets.fa \
    --out_bed primer_targets.bed \
    --out_peaks peaks.tsv \
    --out_qc qc_coverage_summary.tsv
  """
}
