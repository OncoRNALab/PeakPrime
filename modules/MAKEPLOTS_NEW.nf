process MAKEPLOTS_NEW {
  // Use environment file to avoid r-essentials metapackage (100+ deps) which causes slow/hanging conda resolution
  conda "${projectDir}/env/makeplots_env.yml"
  // Alternative inline: conda "r-base=4.3.* r-ggplot2 r-data.table r-optparse r-patchwork r-cowplot bioconductor-genomicranges bioconductor-rtracklayer bioconductor-iranges bioconductor-s4vectors bioconductor-genomicfeatures"
  tag "$gene_id"
  publishDir "${params.outdir}/plots", mode: 'copy'
  errorStrategy 'ignore'  // Continue pipeline even if individual plots fail

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
  path out_name optional true  // Make output optional in case of failure

  script:
  def narrowpeak_arg = narrowpeak_file.name != 'NO_FILE' ? "--narrowpeak ${narrowpeak_file}" : ""
  """
  set +e  # Don't exit on error
  
  Rscript ${projectDir}/bin/MakePlots_new.R \
    --gene ${gene_id} \
    --bw ${bw} \
    --gtf ${gtf} \
    --peaks ${peaks_tsv} \
    --primer ${primer_bed} \
    --qc ${qc_tsv} \
    ${narrowpeak_arg} \
    --out ${out_name}
  
  exit_code=\$?
  if [ \$exit_code -ne 0 ]; then
    echo "Warning: Plot generation failed for gene ${gene_id} (exit code: \$exit_code)"
    echo "Creating empty placeholder file to continue pipeline"
    touch ${out_name}
  fi
  
  exit 0  # Always exit successfully to continue pipeline
  """
}
