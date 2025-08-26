#!/bin/bash
#PBS -N PrimerDesign
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -l mem=32gb

# setting up environment variables
#these are the paths where the Nextflow and Singularity caches will be stored
export APPTAINER_TMPDIR=/user/gent/446/vsc44685/ScratchVO_dir/apptainer_cache/tmp
export APPTAINER_CACHEDIR=/user/gent/446/vsc44685/ScratchVO_dir/singularity
export SINGULARITY_CACHEDIR=/user/gent/446/vsc44685/ScratchVO_dir/singularity
export NXF_SINGULARITY_CACHEDIR=/user/gent/446/vsc44685/ScratchVO_dir/apptainer_cache
export CONDA_PKGS_DIRS=/user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs
export NXF_CONDA_CACHEDIR=/user/gent/446/vsc44685/ScratchVO_dir/conda_cache

#activate a conda enviroment
#conda activate /user/gent/446/vsc44685/DataVO_dir/miniconda2/envs/nf-core 
cd /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR
#load nextflow
ml Nextflow/25.04.4
# 1) Launch with conda profile
# nextflow run main.nf -profile local \
#   --bam ./testdata/Merged_S7_S12.unique.bam \
#   --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
#   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
#   --genes ./testdata/Class1_filt.txt \
#   --pad 100 \
#   --smooth_k 31 \
#   --sliding_window true \
#   --min_window_mean_pct 30 \
#   --max_gap 10 \
#   --trim_low_coverage_pct 20 \
#   --min_exonic_fraction 0.9 \
#   --trim_to_exon true \
#   --primer3_settings config/primer3_settings.txt \
#   --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
#   --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
#   --outdir results/Class1_W30_pad100

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class2_filt.txt \
  --pad 100 \
  --smooth_k 31 \
  --sliding_window true \
  --min_window_mean_pct 30 \
  --max_gap 10 \
  --trim_low_coverage_pct 20 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/Class2_W30_pad100

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class3_filt.txt \
  --pad 100 \
  --smooth_k 31 \
  --sliding_window true \
  --min_window_mean_pct 30 \
  --max_gap 10 \
  --trim_low_coverage_pct 20 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/Class3_W30_pad100

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class4_filt.txt \
  --pad 100 \
  --smooth_k 31 \
  --sliding_window true \
  --min_window_mean_pct 30 \
  --max_gap 10 \
  --trim_low_coverage_pct 20 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/Class4_W30_pad100

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class5_Cov.txt \
  --pad 100 \
  --smooth_k 31 \
  --sliding_window true \
  --min_window_mean_pct 30 \
  --max_gap 10 \
  --trim_low_coverage_pct 20 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/Class5_W30_pad100