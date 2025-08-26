// PeakPrime — RNA‑seq primer design pipeline
// A Nextflow workflow that detects coverage peaks, selects target windows,
// designs cDNA‑appropriate primers with Primer3, and optionally performs
// transcriptome alignment QC.
nextflow.enable.dsl=2

// Parameters are now loaded from params.config

process BAM_INDEX {
  tag "$bam.baseName"
  publishDir params.outdir, mode: 'copy', pattern: '*.bai'
  conda "samtools=1.19"

  input:
  path bam

  output:
  path "${bam}.bai"

  script:
  """
  samtools index -@ 4 ${bam}
  """
}

process MEGADEPTH_BW {
  tag "$bam.baseName"
  publishDir params.outdir, mode: 'copy', pattern: '*.bw'
  conda "megadepth=1.2.* samtools=1.19"

  input:
  path bam
  path bai

  output:
  path "${bam.baseName}.bw"

  script:
  """
  megadepth ${bam} --bigwig --prefix ${bam.baseName}
  # List files to see what was actually created
  ls -la *.bw
  # Find the correct output file and rename it
  if [ -f "${bam.baseName}.all.bw" ]; then
    mv ${bam.baseName}.all.bw ${bam.baseName}.bw
  elif [ -f "${bam.baseName}.coverage.bw" ]; then
    mv ${bam.baseName}.coverage.bw ${bam.baseName}.bw
  elif [ -f "coverage.bw" ]; then
    mv coverage.bw ${bam.baseName}.bw
  else
    echo "Available files:"
    ls -la
    exit 1
  fi
  """
}

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

process MAKE_PRIMER3_INPUT {
  tag "$targets.baseName"
  publishDir params.outdir, mode: 'copy'
  conda "python>=3.10 primer3=2.6.*"

  input:
  path targets
  path p3_settings

  output:
  path "primer3_input.txt"

  script:
  """
  python ${projectDir}/bin/fasta_to_primer3.py --fasta ${targets} --settings ${p3_settings} --out primer3_input.txt
  """
}

process RUN_PRIMER3 {
  tag 'primer3'
  publishDir params.outdir, mode: 'copy'
  conda "python>=3.10 primer3=2.6.*"

  input:
  path p3in

  output:
  path 'primer3_output.txt'

  script:
  """
  primer3_core < ${p3in} > primer3_output.txt
  """
}

process EXTRACT_CDNA_PRIMERS {
  tag 'extract_cdna_primers'
  publishDir params.outdir, mode: 'copy'
  conda "r-base>=4.3 r-essentials r-optparse"

  input:
  path primer3_output
  path peaks_tsv

  output:
  path 'cdna_primers.tsv'

  script:
  """
  Rscript ${projectDir}/bin/extract_cdna_primers.R \
    --primer3_output ${primer3_output} \
    --peaks_tsv ${peaks_tsv} \
    --out_tsv cdna_primers.tsv
  """
}

process PRIMERS_TO_FASTA {
  tag 'primers_to_fasta'
  publishDir params.outdir, mode: 'copy'
  conda "r-base>=4.3 r-essentials bioconductor-biostrings r-optparse"

  input:
  path primers_tsv

  output:
  path 'primers_for_alignment.fa'

  script:
  def max_primers_arg = params.max_primers_per_gene ? "--max_primers_per_gene ${params.max_primers_per_gene}" : ""
  """
  Rscript ${projectDir}/bin/primers_to_fasta.R \
    --primers_tsv ${primers_tsv} \
    --out_fasta primers_for_alignment.fa \
    ${max_primers_arg}
  """
}

process ALIGN_PRIMERS_TRANSCRIPTOME {
  tag 'bowtie2_alignment'
  publishDir params.outdir, mode: 'copy'
  conda "bowtie2=2.5.* samtools=1.19"

  input:
  path primers_fasta
  val transcriptome_index_prefix

  output:
  path 'primers_alignment.bam'
  path 'primers_alignment.bam.bai'
  path 'alignment_stats.txt'

  when:
  transcriptome_index_prefix != 'NO_INDEX'

  script:
  """
  # Align primers to transcriptome - report all alignments
  bowtie2 -f -x ${transcriptome_index_prefix} -U ${primers_fasta} -S primers_alignment.sam -a 2> alignment_stats.txt
  
  # Convert to BAM and sort
  samtools view -bS primers_alignment.sam | samtools sort -o primers_alignment.bam
  
  # Index BAM file
  samtools index primers_alignment.bam
  
  # Clean up intermediate files
  rm primers_alignment.sam
  """
}

process ANALYZE_PRIMER_ALIGNMENTS {
  tag 'alignment_analysis'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.9 biopython=1.79 pandas=1.5.* pysam=0.21.*"

  input:
  path alignment_bam
  path alignment_bai
  path primers_tsv
  path transcriptome_fasta

  output:
  path 'primer_alignment_report.tsv'
  path 'primer_alignment_summary.tsv'

  script:
  def transcriptome_arg = transcriptome_fasta.name != 'NO_FILE' ? "--transcriptome_fasta ${transcriptome_fasta}" : ""
  """
  python ${projectDir}/bin/analyze_primer_alignments.py \
    --alignment_bam ${alignment_bam} \
    --primers_tsv ${primers_tsv} \
    ${transcriptome_arg} \
    --out_report primer_alignment_report.tsv \
    --out_summary primer_alignment_summary.tsv
  """
}

process SELECT_BEST_PRIMERS {
  tag 'select_best_primers'
  publishDir params.outdir, mode: 'copy'
  conda "python=3.9 pandas=1.5.*"

  input:
  path report
  path summary

  output:
  path 'best_primers.tsv'

  script:
  """
  python ${projectDir}/bin/select_best_primer.py \
    --report ${report} \
    --summary ${summary} \
    --out best_primers.tsv
  """
}

workflow {
  // Validate required parameters
  if (!params.bam) error "Missing required parameter: --bam"
  if (!params.gtf) error "Missing required parameter: --gtf" 
  if (!params.genes) error "Missing required parameter: --genes"
  
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

  // Handle optional transcriptome FASTA for gene name mapping
  if (params.transcriptome_fasta) {
    transcriptome_fasta_ch = Channel.fromPath(params.transcriptome_fasta, checkIfExists: true)
  } else {
    transcriptome_fasta_ch = Channel.value(file('NO_FILE'))
  }

  main:
    // Ensure BAM is indexed
    index_bam = BAM_INDEX(bam_ch)

    // BigWig from BAM
    bw = MEGADEPTH_BW(bam_ch, index_bam)

    // Isoform‑agnostic peak picking (union of exons)
    PICK_PEAKS(bw, gtf_ch, genes_ch, fasta_ch)
    primer_targets_fa = PICK_PEAKS.out[0]
    primer_targets_bed = PICK_PEAKS.out[1]
    peaks_tsv = PICK_PEAKS.out[2]
    qc_summary = PICK_PEAKS.out[3]

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
        transcriptome_fasta_ch
      )
      // Select best primers based on alignment summary
      best_primers = SELECT_BEST_PRIMERS(ANALYZE_PRIMER_ALIGNMENTS.out[0], ANALYZE_PRIMER_ALIGNMENTS.out[1])
    }

  emit:
    bw
    primer_targets_fa
    primer_targets_bed
    peaks_tsv
    qc_summary
    primer3_output
    cdna_primers
    primers_fasta
  best_primers
}