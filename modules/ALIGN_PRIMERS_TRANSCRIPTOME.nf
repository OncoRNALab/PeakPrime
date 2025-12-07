process ALIGN_PRIMERS_TRANSCRIPTOME {
  tag 'bowtie_alignment'
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
  # Optimized for ~20nt primers to detect 1-3 mismatches:
  # -l: Seed length (default 15bp allows mismatches outside seed region)
  # -n 3: Allow up to 3 mismatches in the seed region (bowtie1 advantage)
  # -v 3: Allow up to 3 mismatches total across entire primer
  # -a: Report all valid alignments (not just best)
  # --best: Report best alignments when using -a
  bowtie -f -x ${transcriptome_index_prefix} ${primers_fasta} -S primers_alignment.sam -a -l ${params.bowtie_seed_length} -n 3 -v 3 --best 2> alignment_stats.txt
  
  # Convert to BAM and sort
  samtools view -bS primers_alignment.sam | samtools sort -o primers_alignment.bam
  
  # Index BAM file
  samtools index primers_alignment.bam
  
  # Clean up intermediate files
  rm primers_alignment.sam
  """
}
