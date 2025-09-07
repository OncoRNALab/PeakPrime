#!/usr/bin/env nextflow
// workflows/makeplots.nf
nextflow.enable.dsl=2

include { MAKEPLOTS_NEW } from '../modules/MAKEPLOTS_NEW.nf'

workflow makeplots {
  def bw_path = params.bw
  def gtf_path = params.gtf
  def peaks_tsv_path = params.peaks_tsv
  def primer_bed_path = params.primer ?: params.primer_targets_bed  // Support both parameter names
  def qc_tsv_path = params.qc_tsv
  def narrowpeak_path = params.narrowpeak ?: 'NO_FILE'  // Default to NO_FILE if not provided

  // Use the same channel logic as primer_design.nf
  genes_ch = Channel.fromPath(params.genes, checkIfExists: true)

  // Read gene IDs from file, one per line
  gene_plot_ch = genes_ch.splitText()
    .map { it.trim() }
    .filter { it }
    .map { gene_id ->
      def out_name = "plot_${gene_id}.png"
      tuple(gene_id, out_name)
    }

  MAKEPLOTS_NEW(
    file(bw_path),
    file(gtf_path),
    file(peaks_tsv_path),
    file(primer_bed_path),
    file(qc_tsv_path),
    gene_plot_ch.map{ it[0] },
    gene_plot_ch.map{ it[1] },
    file(narrowpeak_path)
  )
}
