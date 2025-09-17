// PeakPrime — RNA‑seq primer design pipeline
// A Nextflow workflow that detects coverage peaks, selects target windows,
// designs cDNA‑appropriate primers with Primer3, and optionally performs
// transcriptome alignment QC.
nextflow.enable.dsl=2

include { makeplots } from './workflows/makeplots.nf'
include { primer_design } from './workflows/primer_design.nf'

workflow {
  if (params.makeplots && !params.bam) {
    // Standalone plotting mode using pre-existing files
    makeplots ()
  } else {
    // Run primer design (which will include plotting if --makeplots is enabled)
    primer_design()
  }
}