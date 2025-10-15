# Reproducibility Guide for PeakPrime Pipeline

## 🐛 Known Issue: MACS2 NumPy Compatibility Error

### Problem Description
When running the pipeline on different systems or with different Conda versions, you may encounter this error:

```
ValueError: numpy.ndarray size changed, may indicate binary incompatibility. 
Expected 96 from C header, got 88 from PyObject
```

### Root Cause
This is a **binary incompatibility** between MACS2 2.2.7.1 (compiled with older NumPy) and newer NumPy versions (1.24+). Different Conda versions resolve dependencies differently, leading to incompatible package combinations even when version numbers are specified.

---

## ✅ Solution: Environment File Approach

As of the latest update, this pipeline uses a dedicated conda environment file (`env/macs2_env.yml`) to ensure reproducible package versions across different systems and Conda versions.

### What Changed
- **Before**: Inline conda specification in `MACS2_CALLPEAK.nf` process
- **After**: Dedicated environment YAML file with strict version pins

### How It Works
The pipeline now uses `env/macs2_env.yml` which locks:
- Python 3.9.*
- NumPy 1.21.6 (compatible with MACS2 2.2.7.1)
- SciPy 1.7.3
- MACS2 2.2.7.1

---

## 🔧 Troubleshooting

### If You Still Get the Error

#### Option 1: Clear Conda Cache (Recommended)
```bash
# Clear conda cache to force fresh environment resolution
conda clean --all -y

# Remove any existing MACS2 environments
conda env list | grep macs2
conda env remove -n <env_name_if_exists>

# Re-run the pipeline
nextflow run main.nf [your parameters]
```

#### Option 2: Pre-create the Environment
```bash
# Manually create and test the environment
conda env create -f env/macs2_env.yml -n test_macs2

# Activate and test MACS2
conda activate test_macs2
macs2 --version
python -c "import numpy; print(numpy.__version__)"

# Should output:
# macs2 2.2.7.1
# 1.21.6
```

#### Option 3: Use Mamba Instead of Conda
```bash
# Install mamba (faster and more reliable dependency solver)
conda install -n base -c conda-forge mamba

# Update nextflow.config to use mamba
# Change: conda.useMamba = false
# To:     conda.useMamba = true
```

#### Option 4: Switch to Inline Specification
If the YAML file approach doesn't work, you can switch back to inline specification:

Edit `modules/MACS2_CALLPEAK.nf`:
```groovy
# Comment out the YAML line:
// conda "${projectDir}/env/macs2_env.yml"

# Uncomment the inline specification:
conda 'python=3.9.* bioconda::macs2=2.2.7.1 bioconda::homer conda-forge::numpy=1.21.6 conda-forge::scipy=1.7.3'
```

---

## 🔍 Verifying Your Environment

### Check Conda Version
```bash
conda --version
# Should be: conda 4.10.0 or higher
```

### Check Active Conda Channels
```bash
conda config --show channels
# Should include: conda-forge, bioconda, defaults
```

### Inspect Nextflow Conda Environment
```bash
# Find where Nextflow creates conda environments
ls -la $HOME/.conda/envs/
# Or check your work directory
ls -la work/*/*/.conda/

# Activate one and check packages
conda activate <path_to_env>
conda list | grep -E "numpy|macs2|scipy"
```

---

## 📊 Expected Package Versions

When the pipeline runs correctly, you should see these versions:

| Package | Version | Build | Channel |
|---------|---------|-------|---------|
| python | 3.9.* | * | conda-forge |
| numpy | 1.21.6 | * | conda-forge |
| scipy | 1.7.3 | * | conda-forge |
| macs2 | 2.2.7.1 | * | bioconda |

---

## 🚀 Future-Proofing

### Consider Upgrading to MACS2 3.x
If you continue to have issues, consider using MACS2 3.0+, which is compatible with newer NumPy:

```yaml
# env/macs2_env_v3.yml
name: macs2_env_v3
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.9
  - numpy>=1.23
  - macs2>=3.0.0
  - homer
```

**Note**: Test thoroughly as MACS2 3.x may have different default parameters or behavior.

---

## 📝 Reporting Issues

If you encounter problems:

1. **Collect environment information**:
   ```bash
   conda info
   conda list
   python --version
   ```

2. **Check Nextflow work directory**:
   ```bash
   ls -R work/
   cat work/<failed_task_hash>/.command.err
   ```

3. **Share error logs** with the full stack trace

---

## 🔗 Related Resources

- [Conda Documentation](https://docs.conda.io/)
- [MACS2 GitHub](https://github.com/macs3-project/MACS)
- [NumPy Version Compatibility Guide](https://numpy.org/doc/stable/dev/depending_on_numpy.html)
- [Nextflow Conda Integration](https://www.nextflow.io/docs/latest/conda.html)

---

**Last Updated**: October 2025
**Pipeline Version**: Compatible with PeakPrime v1.0+
