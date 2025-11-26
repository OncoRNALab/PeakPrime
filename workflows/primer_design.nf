// workflows/primer_design.nf
nextflow.enable.dsl=2

include { MACS2_CALLPEAK } from '../modules/MACS2_CALLPEAK.nf'
include { PROCESS_MACS2_PEAKS } from '../modules/PROCESS_MACS2_PEAKS.nf'
include { BAM_INDEX } from '../modules/BAM_INDEX.nf'
include { MEGADEPTH_BW } from '../modules/MEGADEPTH_BW.nf'
include { MAKE_PRIMER3_INPUT } from '../modules/MAKE_PRIMER3_INPUT.nf'
include { RUN_PRIMER3 } from '../modules/RUN_PRIMER3.nf'
include { EXTRACT_CDNA_PRIMERS } from '../modules/EXTRACT_CDNA_PRIMERS.nf'
include { PRIMERS_TO_FASTA } from '../modules/PRIMERS_TO_FASTA.nf'
include { ALIGN_PRIMERS_TRANSCRIPTOME } from '../modules/ALIGN_PRIMERS_TRANSCRIPTOME.nf'
include { ANALYZE_PRIMER_ALIGNMENTS } from '../modules/ANALYZE_PRIMER_ALIGNMENTS.nf'
include { SELECT_BEST_PRIMERS } from '../modules/SELECT_BEST_PRIMERS.nf'
include { OPTIMIZE_PRIMER_ISOFORMS } from '../modules/OPTIMIZE_PRIMER_ISOFORMS.nf'
include { OPTIMIZE_PRIMERS_MULTIPEAK } from '../modules/OPTIMIZE_PRIMERS_MULTIPEAK.nf'
include { FILTER_GTF } from '../modules/FILTER_GTF.nf'
include { MAKEPLOTS_NEW } from '../modules/MAKEPLOTS_NEW.nf'
include { SUMMARIZE_RESULTS } from '../modules/SUMMARIZE_RESULTS.nf'


workflow primer_design {
  // Validate required parameters
  if (!params.bam) error "Missing required parameter: --bam"
  if (!params.gtf) error "Missing required parameter: --gtf" 
  if (!params.genes) error "Missing required parameter: --genes"
  
  // Validate optional parameters
  if (params.distance_threshold && params.distance_threshold <= 0) {
    error "Parameter --distance_threshold must be positive (current: ${params.distance_threshold})"
  }
  
  // Create input channels
  bam_ch = Channel.fromPath(params.bam, checkIfExists: true)
  gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
  genes_ch = Channel.fromPath(params.genes, checkIfExists: true)
  
  // Handle optional FASTA input
  if (params.fasta) {
    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true)
  } else {
    // Create a dummy file for when no FASTA is provided
    fasta_ch = Channel.value(file('NO_FILE'))
  }

  // Handle optional transcriptome index for primer alignment QC
  if (params.transcriptome_index) {
    // Validate that index files exist
    transcriptome_index_files = Channel.fromPath("${params.transcriptome_index}.*.bt2", checkIfExists: true).collect()
    transcriptome_index_prefix = params.transcriptome_index
  } else {
    transcriptome_index_prefix = 'NO_INDEX'
  }

  // Handle optional transcript mapping file for gene name mapping
  if (params.transcript_mapping && file(params.transcript_mapping).exists()) {
    transcript_mapping_ch = Channel.fromPath(params.transcript_mapping, checkIfExists: true)
  } else {
    transcript_mapping_ch = Channel.value(file('NO_FILE'))
  }

  main:
    // Extract sample ID from BAM filename for MACS2
    sample_id = params.bam.split('/').last().split('\\.')[0]
    bam_with_id = Channel.of([sample_id, file(params.bam)])
    
    // Ensure BAM is indexed (needed for both MACS2 and BigWig)
    index_bam = BAM_INDEX(bam_ch)
    
    // Run MACS2 peak calling
    MACS2_CALLPEAK(bam_with_id)
    
    // Generate BigWig for plotting compatibility
    bw = MEGADEPTH_BW(bam_ch, index_bam)
    
    // Process MACS2 peaks: filter, select best per gene, extract sequences
    PROCESS_MACS2_PEAKS(
      MACS2_CALLPEAK.out.narrowpeak,
      gtf_ch,
      genes_ch,
      fasta_ch
    )
    
    // Use processed peaks as input for primer design
    primer_targets_fa = PROCESS_MACS2_PEAKS.out.fasta
    primer_targets_bed = PROCESS_MACS2_PEAKS.out.bed
    peaks_tsv = PROCESS_MACS2_PEAKS.out.tsv
    qc_summary = PROCESS_MACS2_PEAKS.out.qc

    // Primer3 input + run
    p3in = MAKE_PRIMER3_INPUT(primer_targets_fa, file(params.primer3_settings))
    primer3_output = RUN_PRIMER3(p3in)

    // Extract cDNA-complementary primers
    cdna_primers = EXTRACT_CDNA_PRIMERS(primer3_output, peaks_tsv)

    // Convert primers to FASTA for alignment
    primers_fasta = PRIMERS_TO_FASTA(cdna_primers)

    // Optional: Align primers to transcriptome for specificity check
    if (transcriptome_index_prefix != 'NO_INDEX') {
      alignment_results = ALIGN_PRIMERS_TRANSCRIPTOME(primers_fasta, transcriptome_index_prefix)
      ANALYZE_PRIMER_ALIGNMENTS(
        alignment_results[0], // BAM file
        alignment_results[1], // BAI file  
        cdna_primers,
        transcript_mapping_ch
      )
      // Select best primers based on alignment summary
      best_primers = SELECT_BEST_PRIMERS(ANALYZE_PRIMER_ALIGNMENTS.out[0], ANALYZE_PRIMER_ALIGNMENTS.out[1])
      
      // Choose optimization strategy based on multi-peak mode
      if (params.optimize_multipeak && params.select_all_peaks) {
        // Multi-peak mode: use distance-based scoring across all peaks
        log.info "Multi-peak optimization enabled: scoring primers across all peaks per gene"
        OPTIMIZE_PRIMERS_MULTIPEAK(
          best_primers,
          ANALYZE_PRIMER_ALIGNMENTS.out[1] // primer_alignment_summary.tsv
        )
        optimized_primers = OPTIMIZE_PRIMERS_MULTIPEAK.out.optimized_primers
      } else {
        // Single-peak mode: optimize for isoform coverage (original behavior)
        OPTIMIZE_PRIMER_ISOFORMS(
          best_primers,
          ANALYZE_PRIMER_ALIGNMENTS.out[1] // primer_alignment_summary.tsv (detailed alignments)
        )
        optimized_primers = OPTIMIZE_PRIMER_ISOFORMS.out.optimized_primers
      }
    } else {
      // No transcriptome alignment - create empty channels for emit section
      best_primers = Channel.empty()
      optimized_primers = Channel.empty()
    } // End of transcriptome alignment block

    // Optional: Generate plots if --makeplots is enabled
    if (params.makeplots) {
      // Create fresh channels for plotting (genes_ch was already consumed)
      genes_plot_ch = Channel.fromPath(params.genes, checkIfExists: true)
      gtf_plot_ch = Channel.fromPath(params.gtf, checkIfExists: true)
      
      // Pre-filter GTF to target genes only for faster loading
      filtered_gtf = FILTER_GTF(gtf_plot_ch, genes_plot_ch)
      
      // Read gene IDs from file again for the plotting channel
      genes_list_ch = Channel.fromPath(params.genes, checkIfExists: true)
      gene_plot_ch = genes_list_ch.splitText()
        .map { it.trim() }
        .filter { it }
        .map { gene_id ->
          def out_name = "plot_${gene_id}.png"
          tuple(gene_id, out_name)
        }

      // Run plotting for each gene using filtered GTF
      // Use combine() to pair each gene with the singleton inputs
      plot_inputs = gene_plot_ch
        .combine(bw)
        .combine(filtered_gtf)
        .combine(peaks_tsv)
        .combine(primer_targets_bed)
        .combine(qc_summary)
        .combine(MACS2_CALLPEAK.out.narrowpeak.map{ it[1] })
        .map { gene_id, out_name, bw_file, gtf_file, peaks_file, bed_file, qc_file, np_file ->
          tuple(gene_id, out_name, bw_file, gtf_file, peaks_file, bed_file, qc_file, np_file)
        }

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

  emit:
    macs2_narrowpeak = MACS2_CALLPEAK.out.narrowpeak
    processed_peaks = PROCESS_MACS2_PEAKS.out.tsv
    bw
    primer_targets_fa
    primer_targets_bed
    peaks_tsv
    qc_summary
    primer3_output
    cdna_primers
    primers_fasta
    best_primers
    optimized_primers
    // Note: optimized_primers will contain either single-peak or multi-peak optimization results
    // depending on params.optimize_multipeak and params.select_all_peaks flags
}
