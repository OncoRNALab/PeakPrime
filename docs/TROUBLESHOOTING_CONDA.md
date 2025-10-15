# Troubleshooting: Conda Environment Issues

## Problem: "Rscript: command not found" Error

If you see this error:
```
Command error:
  .command.sh: line 2: Rscript: command not found
```

This means the conda environment wasn't created or activated properly.

## Solution: Clean Conda Cache and Retry

### Step 1: Clean Nextflow Work Directory
```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR

# Remove work directory
rm -rf work/

# Optional: Clean Nextflow cache
rm -rf .nextflow*
```

### Step 2: Clean Conda Cache
```bash
# Find and remove cached conda environments
nextflow clean -f

# Or manually remove conda cache directory
rm -rf $HOME/.nextflow/conda/*

# Or if using a specific conda cache location
rm -rf work/conda/*
```

### Step 3: Test Conda Environment Creation
```bash
# Test if the environment file is valid
conda env create -f env/fetch_mane_env.yml -n test_fetch_mane

# If successful, test R and biomaRt
conda activate test_fetch_mane
Rscript -e "library(biomaRt); print('Success!')"

# Clean up test environment
conda deactivate
conda env remove -n test_fetch_mane
```

### Step 4: Rerun Workflow
```bash
nextflow run main.nf \
  --distance_mode \
  --genes failed_class4.txt \
  --template_length 300 \
  --outdir results_fixed \
  -resume
```

## Alternative: Use Existing R Installation

If conda continues to have issues, you can modify the module to use the system R installation.

### Option 1: Remove conda directive
Edit `modules/FETCH_MANE_TRANSCRIPTS.nf`:

```groovy
process FETCH_MANE_TRANSCRIPTS {
  tag 'fetch_mane'
  publishDir params.outdir, mode: 'copy'
  // conda "${projectDir}/env/fetch_mane_env.yml"  // Commented out

  input:
  path gene_ids_file

  output:
  path 'mane_transcripts.fasta'
  path 'mane_mapping.tsv'

  script:
  """
  # Load R module (adjust module name for your HPC)
  module load R/4.3.0-foss-2023a  # Adjust version as needed
  
  Rscript ${projectDir}/bin/fetch_mane_transcripts.R \
    --gene-ids ${gene_ids_file} \
    --output-fasta mane_transcripts.fasta \
    --output-mapping mane_mapping.tsv \
    --allow-fallback
  """
}
```

### Option 2: Pre-install R packages
If using system R, make sure biomaRt is installed:

```r
# In R console
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("optparse")
```

## Checking Your HPC R Installation

```bash
# Check if R is available
which R
R --version

# Check if Rscript is available
which Rscript

# Load R module if needed
module avail R
module load R/4.3.0-foss-2023a  # Adjust as needed

# Test R packages
Rscript -e "library(biomaRt)"
Rscript -e "library(optparse)"
```

## Expected Conda Environment Contents

The `fetch_mane_env.yml` should install:
- R base (4.3.x)
- r-optparse (for command-line parsing)
- bioconductor-biomart (for Ensembl queries)

## Common Issues

### Issue 1: Bioconda channel not accessible
```bash
# Add bioconda channel
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Issue 2: Conda environment caching issue
```bash
# Disable conda caching
export NXF_CONDA_CACHEDIR=""

# Or set a fresh cache directory
export NXF_CONDA_CACHEDIR="$PWD/conda_cache"
```

### Issue 3: Package conflicts
```bash
# Update conda
conda update -n base conda

# Clean conda cache
conda clean --all
```

## Testing the R Script Directly

You can test the R script outside of Nextflow:

```bash
# Load R module
module load R/4.3.0-foss-2023a

# Test script
Rscript bin/fetch_mane_transcripts.R \
  --gene-ids failed_class4.txt \
  --output-fasta test_output.fasta \
  --output-mapping test_mapping.tsv \
  --allow-fallback
```

## Nextflow Configuration Options

### Enable conda debugging
Add to `nextflow.config`:

```groovy
conda {
    enabled = true
    createTimeout = '30 min'
    cacheDir = "$PWD/conda_cache"
    useMamba = false
}
```

### Use mamba instead of conda (faster)
```groovy
conda {
    enabled = true
    useMamba = true
}
```

## For HPC/PBS Systems

If running on PBS/Torque system, make sure conda is loaded in your job:

```bash
#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1

# Load conda/miniconda module
module load Miniconda3

# Run Nextflow
cd $PBS_O_WORKDIR
nextflow run main.nf \
  --distance_mode \
  --genes failed_class4.txt \
  --template_length 300 \
  -profile pbs
```

## Contact Information

If issues persist, please provide:
1. Full error message from `.nextflow.log`
2. Output of `conda --version`
3. Output of `which R` and `R --version`
4. HPC system information
