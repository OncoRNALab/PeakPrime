#!/usr/bin/env nextflow
// workflows/makeplots.nf
nextflow.enable.dsl=2

include { MAKEPLOTS_NEW } from '../modules/MAKEPLOTS_NEW.nf'
include { FILTER_GTF } from '../modules/FILTER_GTF.nf'

workflow makeplots {
  // Resolve paths relative to launch directory
  def bw_path = file(params.bw)
  def gtf_path = file(params.gtf)
  def peaks_tsv_path = file(params.peaks_tsv)
  def primer_bed_path = file(params.primer ?: params.primer_targets_bed)  // Support both parameter names
  def qc_tsv_path = file(params.qc_tsv)
  def narrowpeak_path = params.narrowpeak ? file(params.narrowpeak) : file('NO_FILE')  // Default to NO_FILE if not provided

  // Create channels from input files
  genes_ch = Channel.fromPath(params.genes, checkIfExists: true)
  gtf_ch = Channel.value(gtf_path)
  
  // Pre-filter GTF to target genes only for faster loading (100x speedup)
  filtered_gtf = FILTER_GTF(gtf_ch, genes_ch)
  
  // Read gene IDs from file again for the plotting channel
  genes_list_ch = Channel.fromPath(params.genes, checkIfExists: true)
  gene_plot_ch = genes_list_ch.splitText()
    .map { it.trim() }
    .filter { it }
    .map { gene_id ->
      def out_name = "plot_${gene_id}.pdf"
      tuple(gene_id, out_name)
    }

  // Create singleton channels for input files
  bw = Channel.value(bw_path)
  peaks_tsv = Channel.value(peaks_tsv_path)
  primer_bed = Channel.value(primer_bed_path)
  qc_tsv = Channel.value(qc_tsv_path)
  narrowpeak = Channel.value(narrowpeak_path)

  // Use combine() to create Cartesian product: each gene paired with all singleton inputs
  // This enables parallel execution of plotting for all genes simultaneously
  plot_inputs = gene_plot_ch
    .combine(bw)
    .combine(filtered_gtf)
    .combine(peaks_tsv)
    .combine(primer_bed)
    .combine(qc_tsv)
    .combine(narrowpeak)
    .map { gene_id, out_name, bw_file, gtf_file, peaks_file, bed_file, qc_file, np_file ->
      tuple(gene_id, out_name, bw_file, gtf_file, peaks_file, bed_file, qc_file, np_file)
    }

  // Execute plotting with optimized parallel inputs
  MAKEPLOTS_NEW(
    plot_inputs.map{ it[2] },  // bw
    plot_inputs.map{ it[3] },  // filtered_gtf
    plot_inputs.map{ it[4] },  // peaks_tsv
    plot_inputs.map{ it[5] },  // primer_targets_bed
    plot_inputs.map{ it[6] },  // qc_summary
    plot_inputs.map{ it[0] },  // gene_id
    plot_inputs.map{ it[1] },  // out_name
    plot_inputs.map{ it[7] }   // narrowpeak
  )
}
