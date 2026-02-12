# PeakPrime manuscript companion instructions

> **Publication:** _[Title Placeholder]_  \
> **Citation:** _[Citation Placeholder]_  \
> **Contact:** _[Corresponding author / email placeholder]_

## License

This project is licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE.md). This means you are free to use, modify, and distribute this software for noncommercial purposes only. See the [LICENSE.md](LICENSE.md) file for full terms.

## Purpose

This branch provides the reference data, and intermediate outputs used to reproduce the analyses described in the publication above.

## Repository Layout

```
PeakPrime/
├── data/            Primer panels (C3 and C4), IMR32/UHRR random-primed counts, targeted counts, target annotations.
├── results/
│   ├── performance_tests/   Nextflow wall-time and memory usage reports for 50-, 100-, and 200-gene primer design runs.
│   ├── figures/             Manuscript-ready plots reproduced from the scripts directory.
│   └── selected_primers/    Selected primer panels.
└── scripts/         R Markdown notebooks for figures: 01 (Figure 2 barplots), 02 (Figure 3 rlog correlations), 03 (Figure 4 fold-changes), plus utils.R helpers.
```

## Reproducing the Analysis

1. **Data Retrieval**
   - Download the raw FASTQ archives from Zenodo (**"PeakPrime manuscript data"**, DOI: [10.5281/zenodo.17821004](https://doi.org/10.5281/zenodo.17821004))
   - Unpack the FASTQ files into `data/raw_fastq/` (or another documented location).


2. **Reference Resources**
   - Retrieve the GRCh38 primary genome assembly and matching GTF annotation:
     ```bash
     mkdir -p data/reference
     cd data/reference
     wget -c https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
     wget -c https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
     gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
     gunzip -k Homo_sapiens.GRCh38.109.gtf.gz
     ```
   - Download the transcriptome FASTA and build the Bowtie2 index used for primer specificity checks:
     ```bash
     wget -c https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
     gunzip -k Homo_sapiens.GRCh38.cdna.all.fa.gz
     bowtie2-build Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.109.transcriptome
     ```
   - Adjust paths and release numbers if different Ensembl snapshots are required for the manuscript.

3. **Install PeakPrime Pipeline**
   - Follow the installation instructions from the repository’s `main` branch (environment setup, dependencies, and optional container caches).
   - Validate the installation by running the bundled smoke test or `nextflow run main.nf --help`.

4. **Generate Alignment BAMs for coverage-based primer design**
   - Process the FASTQs with any QuantSeq-style pipeline to obtain strand-specific, coordinate-sorted BAM files.
   - Accepted options include commercial solutions such as the Lexogen QuantSeqPool pipeline or the open-source implementation maintained by OncoRNALab (`https://github.com/OncoRNALab/QSP_nextflow.git`).
   - Store the resulting BAM files under `data/bam/` and note the exact pipeline version used.
   - The primary primer-design analyses in this branch used the RNA033258 sample for alignment and peak calling.

5. **Input Preparation**
   - Place the prepared BAM file, references files (genome fasta and gtf, and bowite index) and gene lists in the `data/` directory.

6. **Workflow Execution**
    - Launch the main pipeline to reproduce the manuscript primer sets (ensure the file paths below exist exactly as shown or update them accordingly):
       ```bash
       nextflow run main.nf \
          --bam ./data/RNA033258_S1.bam \
          --fasta ./data/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
          --gtf ./data/Homo_sapiens.GRCh38.109.gtf \
          --genes data/genes_list.txt \
          --macs2_qvalue_threshold 0.1 \
          --select_all_peaks --optimize_multipeak \
          --macs2_extsize 150 \
          --macs2_shift 0 \
          --transcriptome_index ./data/bowtie2_index \
          --outdir ./results/ \
          -profile local
       ```
    - Record the Nextflow log and execution report for reproducibility.

7. **Post-processing**

   - Run the R scripts in ./scripts/ to reproduce the figures in the manuscript


### Sample Overview


| RNA_ID    | Filename                          | Class          |
| --------- | --------------------------------- | -------------- |
| RNA028684 | RNA028684_S8_L001_R1_001.fastq.gz | UHRR_random    |
| RNA028684 | RNA028684_S8_L001_R2_001.fastq.gz | UHRR_random    |
| RNA028685 | RNA028685_S9_L001_R1_001.fastq.gz | UHRR_random    |
| RNA028685 | RNA028685_S9_L001_R2_001.fastq.gz | UHRR_random    |
| RNA033258 | RNA033258_S1_R1_001.fastq.gz      | IMR32_random   |
| RNA033258 | RNA033258_S1_R2_001.fastq.gz      | IMR32_random   |
| RNA033259 | RNA033259_S2_R1_001.fastq.gz      | IMR32_random   |
| RNA033259 | RNA033259_S2_R2_001.fastq.gz      | IMR32_random   |
| RNA033944 | RNA033944_S9_R1_001.fastq.gz      | UHRR_targeted  |
| RNA033944 | RNA033944_S9_R2_001.fastq.gz      | UHRR_targeted  |
| RNA033945 | RNA033945_S10_R1_001.fastq.gz     | UHRR_targeted  |
| RNA033945 | RNA033945_S10_R2_001.fastq.gz     | UHRR_targeted  |
| RNA033948 | RNA033948_S13_R1_001.fastq.gz     | UHRR_targeted  |
| RNA033948 | RNA033948_S13_R2_001.fastq.gz     | UHRR_targeted  |
| RNA033949 | RNA033949_S14_R1_001.fastq.gz     | UHRR_targeted  |
| RNA033949 | RNA033949_S14_R2_001.fastq.gz     | UHRR_targeted  |
| RNA033964 | RNA033964_S29_R1_001.fastq.gz     | IMR32_targeted |
| RNA033964 | RNA033964_S29_R2_001.fastq.gz     | IMR32_targeted |
| RNA033965 | RNA033965_S30_R1_001.fastq.gz     | IMR32_targeted |
| RNA033965 | RNA033965_S30_R2_001.fastq.gz     | IMR32_targeted |
| RNA033968 | RNA033968_S33_R1_001.fastq.gz     | IMR32_targeted |
| RNA033968 | RNA033968_S33_R2_001.fastq.gz     | IMR32_targeted |
| RNA033969 | RNA033969_S34_R1_001.fastq.gz     | IMR32_targeted |
| RNA033969 | RNA033969_S34_R2_001.fastq.gz     | IMR32_targeted |

