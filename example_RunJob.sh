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
#go to the directory where the Nextflow script is located
cd /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR
#load nextflow
ml Nextflow/25.04.4
  
  # Test new MACS2-based approach (generates BigWig + MACS2 results)
  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class1_filt.txt \
  --macs2_pvalue_threshold 0.01 \
  --macs2_min_peak_score 10 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/macs2_class1

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class2_filt.txt \
  --macs2_pvalue_threshold 0.01 \
  --macs2_min_peak_score 10 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/macs2_class2v2

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class4_filt.txt \
  --macs2_pvalue_threshold 0.01 \
  --macs2_min_peak_score 10 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/macs2_class4v2

  nextflow run main.nf -profile local \
  --bam ./testdata/Merged_S7_S12.unique.bam \
  --fasta /data/gent/vo/000/gvo00027/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa \
  --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --genes ./testdata/Class5_filt.txt \
  --macs2_pvalue_threshold 0.01 \
  --macs2_min_peak_score 10 \
  --min_exonic_fraction 0.9 \
  --trim_to_exon true \
  --primer3_settings config/primer3_settings.txt \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --transcriptome_fasta /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.fa \
  --outdir results/macs2_class5

  ##make plots for MACS2 results
  nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class1_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class1/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class1/processed_peaks/selected_peaks.tsv \
   --primer_targets_bed /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class1/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class1/processed_peaks/peaks_qc_summary.tsv \
   --outdir results/macs2_class1/plots
  
  nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class2_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2/processed_peaks/selected_peaks.tsv \
   --primer_targets_bed /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2/processed_peaks/peaks_qc_summary.tsv \
   --outdir results/macs2_class2/plots -resume

  nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class3_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class3/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class3/processed_peaks/selected_peaks.tsv \
   --primer_targets_bed /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class3/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class3/processed_peaks/peaks_qc_summary.tsv \
   --outdir results/macs2_class3/plots

nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class4_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class4/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class4/processed_peaks/selected_peaks.tsv \
   --primer_targets_bed /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class4/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class4/processed_peaks/peaks_qc_summary.tsv \
   --outdir results/macs2_class4/plots -resume

nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class5_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class5/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class5/processed_peaks/selected_peaks.tsv \
   --primer_targets_bed /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class5/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class5/processed_peaks/peaks_qc_summary.tsv \
   --outdir results/macs2_class5/plots

nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class2_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/selected_peaks.tsv \
   --primer /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/peaks_qc_summary.tsv \
   --narrowpeak /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/macs2_peaks/Merged_S7_S12_peaks.narrowPeak \
   --outdir results/macs2_class2v2/plots

nextflow run main.nf -profile local --makeplots \
   --gtf /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
   --genes ./testdata/Class2_filt.txt \
   --bw /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/Merged_S7_S12.unique.bw \
   --peaks_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/selected_peaks.tsv \
   --primer /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/selected_peaks.bed \
   --qc_tsv /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/processed_peaks/peaks_qc_summary.tsv \
   --narrowpeak /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/results/macs2_class2v2/macs2_peaks/Merged_S7_S12_peaks.narrowPeak \
   --outdir results/macs2_class2v2/plots