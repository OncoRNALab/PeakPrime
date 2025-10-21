# Troubleshooting Guide

## Quick Diagnosis

| Error Message | Issue | Section |
|---------------|-------|---------|
| `Rscript: command not found` | Conda environment not created | [Conda Environment Issues](#conda-environment-issues) |
| `numpy.dtype size changed` | MACS2/NumPy incompatibility | [MACS2 NumPy Error](#macs2-numpy-compatibility) |
| `No such file or directory` | File path resolution issue | [File Path Issues](#file-path-issues) |
| `Out of memory` | Insufficient RAM | [Memory Issues](#memory-issues) |
| Environment hangs during creation | Network/cache issues | [Conda Cache Cleanup](#conda-cache-cleanup) |

---

## Conda Environment Issues

### Problem: "Rscript: command not found"

**Error:**
```
Command error:
  .command.sh: line 2: Rscript: command not found
```

**Cause:** Conda environment wasn't created or activated properly.

### Solution 1: Clean and Retry

```bash
# 1. Clean Nextflow work directory
rm -rf work/
rm -rf .nextflow*

# 2. Clean conda cache
nextflow clean -f
rm -rf $HOME/.nextflow/conda/*

# 3. Rerun workflow
nextflow run main.nf [your parameters] -resume
```

### Solution 2: Test Environment Manually

```bash
# Test if environment file is valid
conda env create -f env/filter_gtf_env.yml -n test_env

# If successful, test R
conda activate test_env
Rscript -e "library(rtracklayer); print('Success!')"

# Clean up
conda deactivate
conda env remove -n test_env
```

### Solution 3: Use System R Installation

If conda issues persist, modify the process to use system R:

**Edit module (e.g., `modules/FILTER_GTF.nf`):**
```groovy
process FILTER_GTF {
  // conda "${projectDir}/env/filter_gtf_env.yml"  // Comment out
  
  script:
  """
  module load R/4.3.0  # Load system R instead
  
  Rscript ${projectDir}/bin/filter_gtf.R ...
  """
}
```

---

## Conda Cache Cleanup

### The Problem

Nextflow uses custom conda cache directories that aren't cleaned by standard `conda clean` commands.

**Custom cache locations:**
- `NXF_CONDA_CACHEDIR="/user/gent/446/vsc44685/ScratchVO_dir/conda_cache"`
- `CONDA_PKGS_DIRS="/user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs"`

### Solution 1: Use Cleanup Script (Recommended)

```bash
# Interactive cleanup (asks for confirmation)
bash bin/clean_conda_cache.sh

# Force cleanup without confirmation
bash bin/clean_conda_cache.sh --force
```

**What it cleans:**
- ✅ Standard conda cache (`conda clean --all`)
- ✅ Nextflow conda cache (`$NXF_CONDA_CACHEDIR`)
- ✅ Conda packages cache (`$CONDA_PKGS_DIRS`)
- ✅ Test environments (optional)

### Solution 2: Manual Cleanup

```bash
# 1. Clean standard conda cache
conda clean --all -y

# 2. Clean Nextflow's conda cache
rm -rf /user/gent/446/vsc44685/ScratchVO_dir/conda_cache/*

# 3. Clean conda packages cache
rm -rf /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs/*
```

### Verification

```bash
# Check cache sizes (should be 0 or minimal)
du -sh /user/gent/446/vsc44685/ScratchVO_dir/conda_cache
du -sh /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs

# List contents (should be empty)
ls -la /user/gent/446/vsc44685/ScratchVO_dir/conda_cache/
```

### When to Clean Caches

**Clean when:**
- ❌ Getting conda/environment errors
- 🔄 After updating environment files
- 🐛 Pipeline fails with "environment activation" errors
- 💾 Running out of disk space
- 🧪 Testing reproducibility

**Don't clean when:**
- ✅ Pipeline is working fine (caches save time!)
- ⚡ You need fast subsequent runs

---

## MACS2 NumPy Compatibility

### Problem: "numpy.dtype size changed"

**Error:**
```python
ValueError: numpy.dtype size changed, may indicate binary incompatibility.
Expected 96 from C header, got 88 from PyObject
```

**Cause:** MACS2 compiled against different NumPy version.

### Solution: Use Locked Environment

The pipeline includes a locked environment specification with compatible versions:

**File: `env/macs2_env.yml`**
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - numpy=1.23.5
  - macs2=2.2.9.1
```

### Testing

```bash
# Test MACS2 environment
bash bin/test_macs2_environment.sh

# Test with cache cleaning
bash bin/test_macs2_environment.sh --clean
```

### Manual Test

```bash
# Create test environment
conda env create -f env/macs2_env.yml -n test_macs2

# Activate and test
conda activate test_macs2
python -c "import numpy; print(numpy.__version__)"
macs2 --version

# Clean up
conda deactivate
conda env remove -n test_macs2
```

---

## File Path Issues

### Problem: "No such file or directory"

**Error:**
```
ERROR ~ No such file or directory: /kyukon/scratch/.../file.txt
```

**Cause:** Nextflow resolving relative paths incorrectly.

### Solution 1: Use Absolute Paths

```bash
# Instead of relative paths
--genes ./testdata/genes.txt

# Use absolute paths
--genes /full/path/to/testdata/genes.txt
```

### Solution 2: Check File Exists

```bash
# Verify file exists before running
ls -la ./testdata/genes.txt

# Check current directory
pwd

# Verify all required files
ls -la results/test/processed_peaks/*.tsv
```

### Solution 3: Fix Workflow (if issue persists)

For workflows using `Channel.fromPath()`, ensure files are resolved properly:

```groovy
// Use file() helper to resolve paths
def my_file = file(params.input_file)

// Create channel from resolved file
my_channel = Channel.value(my_file)
```

---

## Memory Issues

### Problem: Out of Memory (OOM)

**Error:**
```
ERROR ~ Process `MAKEPLOTS_NEW (ENSG00000123456)` terminated with an error exit status (137)
```

**Exit code 137 = Out of Memory (killed by system)**

### Solution: Reduce Parallelization

**Edit `nextflow.config`:**
```groovy
process {
  withName: MAKEPLOTS_NEW {
    maxForks = 2  // Reduce from 4 to 2
    memory = '6 GB'
  }
}
```

### Calculate Optimal maxForks

```
maxForks = floor(Available RAM / Process Memory)

Examples:
- 16 GB system → maxForks = 2 (16 / 6 = 2.6)
- 32 GB system → maxForks = 5 (32 / 6 = 5.3)
- 64 GB system → maxForks = 10 (64 / 6 = 10.6)
```

### Alternative: Request More Memory

```groovy
process {
  withName: MAKEPLOTS_NEW {
    memory = { 6.GB * task.attempt }
    maxRetries = 3
    errorStrategy = 'retry'
  }
}
```

---

## Network/Connection Issues

### Problem: Conda Environment Creation Hangs

**Symptoms:**
- Environment creation takes >10 minutes
- No progress in terminal
- Network error messages

### Solution 1: Use Mamba (Faster)

**Edit `nextflow.config`:**
```groovy
conda {
  enabled = true
  useMamba = true  // Add this line
  createTimeout = '30 min'
}
```

**Install mamba:**
```bash
conda install -n base mamba
```

### Solution 2: Check Network

```bash
# Test conda network connection
conda search macs2

# Test bioconda channel
conda search -c bioconda macs2

# Add channels if missing
conda config --add channels conda-forge
conda config --add channels bioconda
```

### Solution 3: Pre-create Environments

```bash
# Create environments manually before running pipeline
conda env create -f env/macs2_env.yml -p ./conda_envs/macs2
conda env create -f env/filter_gtf_env.yml -p ./conda_envs/filter_gtf

# Update nextflow.config to use these
conda.cacheDir = "$PWD/conda_envs"
```

---

## HPC-Specific Issues

### PBS/SLURM Job Issues

**Problem:** Jobs not starting or failing immediately

**Check job status:**
```bash
# SLURM
squeue -u $USER

# PBS
qstat -u $USER

# Check job details
scontrol show job JOBID  # SLURM
qstat -f JOBID          # PBS
```

**Solution: Load required modules in profile**

**Edit `nextflow.config`:**
```groovy
profiles {
  pbs {
    process {
      executor = 'pbs'
      beforeScript = 'module load Miniconda3'  // Add this
    }
  }
}
```

### Quota/Disk Space Issues

**Check disk usage:**
```bash
# Check quota
quota -s

# Check disk usage in work directory
du -sh work/

# Find large files
find work/ -size +1G -exec ls -lh {} \;
```

**Solution:**
```bash
# Clean work directory
rm -rf work/

# Clean conda caches
bash bin/clean_conda_cache.sh --force

# Clean Nextflow cache
nextflow clean -f
```

---

## Common Workflow Errors

### Error: "No peaks found"

**Cause:** No significant peaks called by MACS2.

**Solution:**
```bash
# Check MACS2 output
less results/macs2_peaks/*.narrowPeak

# Adjust q-value threshold (more lenient)
nextflow run main.nf --qvalue_threshold 0.1

# Or use p-value instead
nextflow run main.nf --qvalue_threshold 1.0 --min_peak_score 10
```

### Error: "No primers generated"

**Cause:** Primer3 couldn't design primers for target regions.

**Solution:**
```bash
# Check Primer3 output
less results/primer3_output/*.txt

# Relax Primer3 parameters in config/primer3_settings.txt
# Adjust:
PRIMER_MIN_TM=55  # Lower from 58
PRIMER_MAX_TM=65  # Raise from 62
```

### Error: Gene not in GTF

**Cause:** Gene ID not found in annotation file.

**Solution:**
```bash
# Verify gene IDs in GTF
grep "ENSG00000123456" genome.gtf

# Check GTF format
head -20 genome.gtf

# Ensure gene IDs match (with/without version)
# Bad:  ENSG00000123456.5
# Good: ENSG00000123456
```

---

## Debugging Tools

### Enable Debug Mode

```bash
# Run with debug trace
nextflow run main.nf -with-trace -with-report -with-timeline

# Check generated reports
ls -la trace.txt report.html timeline.html
```

### Check Nextflow Log

```bash
# View recent errors
tail -100 .nextflow.log

# Search for specific errors
grep -i "error" .nextflow.log

# Follow log in real-time
tail -f .nextflow.log
```

### Inspect Failed Process

```bash
# Find failed work directory
find work/ -name ".exitcode" -exec grep -l "^[^0]" {} \; | \
  xargs -I {} dirname {} | head -1

# Navigate to failed work directory
cd work/xx/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Check outputs
cat .command.sh    # Script that ran
cat .command.out   # Standard output
cat .command.err   # Standard error
cat .command.log   # Combined log
```

---

## Getting Help

### Information to Provide

When reporting issues, include:

1. **Full command used:**
   ```bash
   nextflow run main.nf --bam sample.bam --genes genes.txt ...
   ```

2. **Error message:**
   ```bash
   tail -50 .nextflow.log
   ```

3. **System information:**
   ```bash
   nextflow -version
   conda --version
   uname -a
   ```

4. **Process details:**
   ```bash
   # From failed work directory
   cat .command.sh
   cat .command.err
   ```

### Useful Commands

```bash
# Check Nextflow version
nextflow -version

# List available profiles
grep "profiles {" nextflow.config -A 20

# Check conda environments
conda env list

# Test conda environment manually
conda env create -f env/macs2_env.yml -n test
```

---

## Complete Troubleshooting Workflow

```bash
# 1. Clean everything
rm -rf work/ .nextflow*
bash bin/clean_conda_cache.sh --force

# 2. Verify input files exist
ls -la [your input files]

# 3. Test environments
bash bin/test_macs2_environment.sh

# 4. Run with debugging
nextflow run main.nf \
  [your parameters] \
  -with-trace \
  -with-report \
  -with-timeline \
  -resume

# 5. If fails, check logs
tail -100 .nextflow.log
cat trace.txt
```

---

## Related Documentation

- **Setup:** `REPRODUCIBILITY_GUIDE.md`
- **Pipeline:** `pipeline_steps.md`
- **Features:** `features/PEAK_RANKING.md`, `features/DISTANCE_MODE.md`
- **Technical:** `TECHNICAL_APPENDIX.md`

---

*Last updated: October 2025*
