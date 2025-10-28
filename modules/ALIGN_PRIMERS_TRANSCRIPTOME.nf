process ALIGN_PRIMERS_TRANSCRIPTOME {
  tag 'bowtie2_alignment'
  publishDir params.outdir, mode: 'copy'
  conda "${projectDir}/env/align_primers_env.yml"

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
