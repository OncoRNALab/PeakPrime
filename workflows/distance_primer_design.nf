// workflows/distance_primer_design.nf
// Distance-based primer design workflow for 3' end targeting
// Two modes: 1) Fetch MANE transcripts from gene IDs, or 2) Use provided transcript FASTA
nextflow.enable.dsl=2

include { FETCH_MANE_TRANSCRIPTS } from '../modules/FETCH_MANE_TRANSCRIPTS.nf'
include { EXTRACT_3PRIME_SEQUENCE } from '../modules/EXTRACT_3PRIME_SEQUENCE.nf'
include { MAKE_PRIMER3_INPUT } from '../modules/MAKE_PRIMER3_INPUT.nf'
include { RUN_PRIMER3 } from '../modules/RUN_PRIMER3.nf'
include { SIMPLE_PRIMERS_TO_FASTA } from '../modules/SIMPLE_PRIMERS_TO_FASTA.nf'
include { CONVERT_PRIMER3_TO_TSV } from '../modules/CONVERT_PRIMER3_TO_TSV.nf'
include { ALIGN_PRIMERS_TRANSCRIPTOME } from '../modules/ALIGN_PRIMERS_TRANSCRIPTOME.nf'
include { ANALYZE_PRIMER_ALIGNMENTS } from '../modules/ANALYZE_PRIMER_ALIGNMENTS.nf'
include { SELECT_BEST_PRIMERS } from '../modules/SELECT_BEST_PRIMERS.nf'
include { OPTIMIZE_PRIMER_ISOFORMS } from '../modules/OPTIMIZE_PRIMER_ISOFORMS.nf'


workflow distance_primer_design {
  // Validate parameters
  if (!params.genes && !params.transcript_fasta) {
    error "Missing required parameter: either --genes (gene ID list) or --transcript_fasta (transcript sequences)"
  }
  
  if (params.genes && params.transcript_fasta) {
    error "Please provide either --genes OR --transcript_fasta, not both"
  }
  
  if (!params.template_length) {
    error "Missing required parameter: --template_length (number of bases to extract from 3' end)"
  }
  
  if (params.template_length <= 0) {
    error "Parameter --template_length must be positive (current: ${params.template_length})"
  }

  main:
    // Branch 1: Fetch MANE transcripts from gene IDs
    if (params.genes) {
      genes_ch = Channel.fromPath(params.genes, checkIfExists: true)
      
      // Fetch MANE Select transcripts from Ensembl REST API
      fetch_results = FETCH_MANE_TRANSCRIPTS(genes_ch)
      transcript_fasta = fetch_results[0]  // mane_transcripts.fasta
      gene_mapping = fetch_results[1]      // mane_mapping.tsv
      
      println "Mode: Fetching MANE Select transcripts from Ensembl REST API"
    }
    // Branch 2: Use provided transcript FASTA
    else {
      transcript_fasta = Channel.fromPath(params.transcript_fasta, checkIfExists: true)
      gene_mapping = Channel.value(file('NO_FILE'))
      
      println "Mode: Using provided transcript FASTA file"
    }
    
    // Extract N bases from 3' end of transcripts
    threeprime_sequences = EXTRACT_3PRIME_SEQUENCE(
      transcript_fasta,
      params.template_length
    )
    
    // Prepare Primer3 input
    p3in = MAKE_PRIMER3_INPUT(
      threeprime_sequences,
      file(params.primer3_settings)
    )

    // Run Primer3 for primer design
    primer3_output = RUN_PRIMER3(p3in)

    // Convert primers to FASTA for alignment (simple version without gene grouping)
    primers_fasta = SIMPLE_PRIMERS_TO_FASTA(primer3_output)

    // Handle optional transcriptome FASTA for gene name mapping
    if (params.transcriptome_fasta) {
      transcriptome_fasta_ch = Channel.fromPath(params.transcriptome_fasta, checkIfExists: true)
    } else {
      transcriptome_fasta_ch = Channel.value(file('NO_FILE'))
    }

    // Optional: Align primers to transcriptome for specificity checking
    if (params.transcriptome_index) {
      // Validate that index files exist
      transcriptome_index_files = Channel.fromPath("${params.transcriptome_index}.*.bt2", checkIfExists: true).collect()
      
      // Align primers to transcriptome
      alignment_results = ALIGN_PRIMERS_TRANSCRIPTOME(
        primers_fasta,
        params.transcriptome_index
      )
      
      // Convert Primer3 output to TSV format for alignment analysis
      // Note: gene_mapping not needed since primer3 output already has correct gene IDs
      primers_tsv = CONVERT_PRIMER3_TO_TSV(
        primer3_output,
        Channel.value(file('NO_FILE'))
      )
      
      // Analyze primer alignments and generate alignment summary
      ANALYZE_PRIMER_ALIGNMENTS(
        alignment_results[0], // BAM file
        alignment_results[1], // BAI file  
        primers_tsv,          // Converted TSV format
        transcriptome_fasta_ch
      )
      
      // Select best primers based on alignment summary
      best_primers = SELECT_BEST_PRIMERS(
        ANALYZE_PRIMER_ALIGNMENTS.out[0], // alignment summary
        ANALYZE_PRIMER_ALIGNMENTS.out[1]  // filtered primers
      )
      
      // Optimize primer selection to maximize distinct isoform coverage
      optimized_primers = OPTIMIZE_PRIMER_ISOFORMS(
        best_primers,
        ANALYZE_PRIMER_ALIGNMENTS.out[1] // primer_alignment_summary.tsv (detailed alignments)
      )
      
      println "Transcriptome alignment enabled for specificity checking"
    } else {
      println "Transcriptome alignment skipped (no --transcriptome_index provided)"
      best_primers = Channel.empty()
      optimized_primers = Channel.empty()
    }

  emit:
    transcript_fasta
    threeprime_sequences
    primer3_output
    primers_fasta
    gene_mapping
    best_primers
    optimized_primers
}
