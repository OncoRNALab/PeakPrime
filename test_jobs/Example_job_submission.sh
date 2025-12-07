#!/bin/bash
#PBS -N PeakPrime_job_
#PBS -l nodes=1:ppn=3
#PBS -l walltime=04:00:00
#PBS -l mem=16gb

#set the enviromental paths, make sure the folders exist beforehand
export CONDA_PKGS_DIRS=/user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs
export NXF_CONDA_CACHEDIR=/user/gent/446/vsc44685/ScratchVO_dir/conda_cache

# go to the pipeline directory
cd path/to/PeakPrime
#load nextflow module
ml Nextflow/25.04.4

# run the complete pipeline with multipeak mode and primer optimization

nextflow run main.nf --makeplots \
  --bam /path/to/bamfile.bam \
  --fasta /path/to/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gtf /path/to/GTF/Homo_sapiens.GRCh38.109.gtf \
  --genes ./path/to/geneidlist/TargetGeneIDs.txt \
  --macs2_qvalue_threshold 0.1 \
  --select_all_peaks \
  --optimize_multipeak \
  --transcriptome_index /path/to/pre-built/bowtie2_index \
  --outdir ./results \
  -profile pbs # to use the slurm scheduler of the HPC

 nextflow run main.nf \
  --bam testdata/RNA034671_S13_L001_dedup.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes testdata/geneids.txt \
  --macs2_qvalue_threshold 0.1 \
  --select_all_peaks --optimize_multipeak \
  --macs2_extsize 150 \
  --macs2_shift 0 \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir ./resultsC3/ \
  -profile local


#if you have to only make plots after running the pipeline use the following 
# some of these files are produced by the pipeline and should be found in the output directory
nextflow run main.nf --makeplots \
  --bw ./results/File.bw \
  --gtf /path/to/GTF/Homo_sapiens.GRCh38.109.gtf \
  --genes ./path/to/geneidlist/TargetGeneIDs.txt \
  --peaks_tsv ./results/processed_peaks/selected_peaks.tsv \
  --qc_tsv ./results/peaks_qc_summary.tsv \
  --primer_targets_bed ./results/processed_peaks/primer_targets.bed \
  --narrowpeak ./results/macs2_peaks/File.narrowPeak \
  --outdir ./results/ \
  -profile pbs

