// PeakPrime — RNA‑seq primer design pipeline
// A Nextflow workflow that detects coverage peaks, selects target windows,
// designs cDNA‑appropriate primers with Primer3, and optionally performs
// transcriptome alignment QC.
nextflow.enable.dsl=2

include { makeplots } from './workflows/makeplots.nf'
include { primer_design } from './workflows/primer_design.nf'
include { distance_primer_design } from './workflows/distance_primer_design.nf'

workflow {
  if (params.distance_mode) {
    // Distance-based primer design workflow
    // Designs primers at a fixed distance from 3' end
    distance_primer_design()
  } else if (params.makeplots && !params.bam) {
    // Standalone plotting mode using pre-existing files
    makeplots ()
  } else {
    // Run primer design (which will include plotting if --makeplots is enabled)
    primer_design()
  }
}