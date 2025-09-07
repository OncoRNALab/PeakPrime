process MAKEPLOTS_NEW {
  conda "r-base>=4.3 r-essentials r-ggplot2 r-data.table r-optparse r-patchwork r-cowplot bioconductor-genomicranges bioconductor-rtracklayer bioconductor-iranges bioconductor-s4vectors bioconductor-genomicfeatures bioconductor-bsgenome bioconductor-bsgenome.hsapiens.ucsc.hg38"
  tag "$gene_id"
  publishDir params.outdir, mode: 'copy'

  input:
  path bw
  path gtf
  path peaks_tsv
  path primer_bed
  path qc_tsv
  val gene_id
  val out_name
  path narrowpeak_file

  output:
  path out_name

  script:
  def narrowpeak_arg = narrowpeak_file.name != 'NO_FILE' ? "--narrowpeak ${narrowpeak_file}" : ""
  """
  Rscript ${projectDir}/bin/MakePlots_new.R \
    --gene ${gene_id} \
    --bw ${bw} \
    --gtf ${gtf} \
    --peaks ${peaks_tsv} \
    --primer ${primer_bed} \
    --qc ${qc_tsv} \
    ${narrowpeak_arg} \
    --out ${out_name}
  """
}
