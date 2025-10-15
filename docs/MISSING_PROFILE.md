# CRITICAL: Missing Nextflow Profile

## The Problem

When running:
```bash
nextflow run main.nf --distance_mode --genes ... --outdir ...
```

**No conda environment is created** because you didn't specify a profile!

## The Solution

You MUST use `-profile local` or `-profile pbs` to enable conda:

```bash
nextflow run main.nf \
  --distance_mode \
  --genes testdata/failed_class4.txt \
  --template_length 300 \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --outdir results/results_distance \
  -profile local    # <-- THIS IS REQUIRED!
```

## Why This Matters

Looking at your `nextflow.config`:

```groovy
profiles { 
  local {
    conda.enabled = true       # <-- Conda only enabled in profile
    process.conda = true
    conda.useMamba = false
    process.executor = 'local'
  }
  
  pbs {
    conda.enabled = true
    process.conda = true
    process.executor = 'pbs'
  }
}
```

**Without a profile**: No conda, no R, no Rscript → Error 127
**With `-profile local`**: Conda creates environment with R and biomaRt → Success!

## Correct Commands

### For Local Execution:
```bash
nextflow run main.nf \
  --distance_mode \
  --genes testdata/failed_class4.txt \
  --template_length 300 \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --outdir results/results_distance \
  -profile local
```

### For PBS/HPC Submission:
```bash
nextflow run main.nf \
  --distance_mode \
  --genes testdata/failed_class4.txt \
  --template_length 300 \
  --transcriptome_index /data/gent/vo/000/gvo00027/resources/Bowtie2_index/Homo_sapiens/Transcriptome_Homo_sapiens.GRCh38.109.chrIS_spikes_45S/bowtie2_index \
  --outdir results/results_distance \
  -profile pbs
```

## What Happens With Profile

When you use `-profile local`:

1. **Conda is enabled**: `conda.enabled = true`
2. **Environment is created**: From `env/fetch_mane_env.yml`
3. **R is installed**: r-base=4.3.* 
4. **biomaRt is installed**: bioconductor-biomart
5. **Rscript is available**: In the conda environment's PATH
6. **Script runs successfully**: ✅

## Quick Test

After you run with `-profile local`, check the work directory:
```bash
# Find the conda environment path
ls -la work/conda/
```

You should see a conda environment directory created!

## Summary

**Problem**: Running without profile → No conda → No R → Error 127
**Solution**: Always use `-profile local` or `-profile pbs`

This is why it worked before with other workflows - you were probably using the profile then!
