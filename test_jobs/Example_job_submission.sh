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
 --bam /user/gent/446/vsc44685/ScratchVO_dir/IMR32_RPtrimmed/trimmed_4M/3_Star/RNA033258_S1_Aligned.sortedByCoord.out.bam \
  --fasta /path/to/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gtf /path/to/GTF/Homo_sapiens.GRCh38.109.gtf \
  --genes ./path/to/geneidlist/TargetGeneIDs.txt \
  --macs2_qvalue_threshold 0.1 \
  --select_all_peaks \
  --optimize_multipeak \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie_index/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.transcriptome.chrIS_spikes_45S \
  --outdir ./results \
  -profile pbs # to use the slurm scheduler of the HPC

 nextflow run main.nf \
  --bam /user/gent/446/vsc44685/ScratchVO_dir/IMR32_RPtrimmed/trimmed_4M/3_Star/RNA033258_S1_Aligned.sortedByCoord.out.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes testdata/geneids.txt \
  --macs2_qvalue_threshold 0.1 \
  --select_all_peaks --optimize_multipeak \
  --macs2_extsize 150 \
  --macs2_shift 0 \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie_index/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.transcriptome.chrIS_spikes_45S \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir ./results/ \
  -profile local 
