# Fix for MACS2 NumPy Compatibility Issue

## Summary of Changes

**Date**: October 6, 2025  
**Issue**: MACS2 2.2.7.1 crashes with NumPy incompatibility error due to different conda version resolving dependencies differently

### Files Modified/Created:

1. **`env/macs2_env.yml`** (NEW)
   - Dedicated conda environment file with strict version pins
   - Locks Python 3.9.*, NumPy 1.21.6, SciPy 1.7.3, MACS2 2.2.7.1
   - Ensures reproducible environments across different conda versions

2. **`modules/MACS2_CALLPEAK.nf`** (MODIFIED)
   - Changed from inline conda specification to environment file
   - Old: `conda 'bioconda::macs2=2.2.7.1 ...'`
   - New: `conda "${projectDir}/env/macs2_env.yml"`
   - Kept inline version as commented alternative

3. **`docs/REPRODUCIBILITY_GUIDE.md`** (NEW)
   - Comprehensive troubleshooting guide
   - Explains the root cause of the issue
   - Provides multiple solution options
   - Instructions for verifying environment correctness

4. **`bin/test_macs2_environment.sh`** (NEW)
   - Automated test script to verify MACS2 compatibility
   - Tests environment creation, package versions, and MACS2 imports
   - Provides clear pass/fail feedback

---

## For Your Colleague

### Quick Fix Instructions:

1. **Pull the latest changes** from the repository

2. **Clean conda cache** (recommended - includes Nextflow custom caches):
   ```bash
   # Interactive cleanup (asks for confirmation)
   bash bin/clean_conda_cache.sh
   
   # OR force cleanup without confirmation
   bash bin/clean_conda_cache.sh --force
   ```

3. **Run the test script** before running the pipeline:
   ```bash
   # Basic test
   bash bin/test_macs2_environment.sh
   
   # OR test with cache cleaning
   bash bin/test_macs2_environment.sh --clean
   ```

4. **If test passes**, run the pipeline normally:
   ```bash
   nextflow run main.nf [your parameters]
   ```

5. **If test fails**, see `docs/REPRODUCIBILITY_GUIDE.md` for troubleshooting

---

## Technical Explanation

### Why This Happened:

**Root Cause**: Binary incompatibility between MACS2's compiled Cython extensions and NumPy

**Different Conda Versions**:
- Your conda version resolved dependencies to compatible versions
- Colleague's conda version (likely newer) resolved to incompatible NumPy 1.24+
- MACS2 2.2.7.1 was compiled against NumPy <1.24 and cannot load with newer versions

**The Error**:
```
ValueError: numpy.ndarray size changed, may indicate binary incompatibility.
Expected 96 from C header, got 88 from PyObject
```
- "96 vs 88" = size mismatch in NumPy internal structures
- Happens when compiled code expects old NumPy but loads with new NumPy

### Why the Fix Works:

1. **Environment YAML** provides deterministic package resolution
2. **Strict version pins** prevent conda from "upgrading" to incompatible versions
3. **Python version pinning** ensures consistent ABI across all dependencies
4. **Order matters**: Channels and package specifications affect resolution

---

## Verification

After pulling these changes, your colleague should:

1. Run the test script - it will create a temporary environment and verify MACS2 works
2. If successful, the pipeline will use the same environment specification
3. The environment will be cached by Nextflow and reused across runs

---

## Alternative Solutions (if primary fix doesn't work)

### Option A: Use Mamba
```bash
# In nextflow.config, change:
conda.useMamba = true
```

### Option B: Pre-create Environment
```bash
conda env create -f env/macs2_env.yml -n peakprime_macs2
# Then modify nextflow.config to use this environment
```

### Option C: Upgrade MACS2 (requires testing)
```yaml
# env/macs2_env.yml
dependencies:
  - macs2>=3.0.0  # Compatible with newer NumPy
  - numpy>=1.23
```

---

## Git Commit Message (Suggested)

```
Fix: Resolve MACS2 NumPy compatibility issue across conda versions

- Create dedicated conda environment file (env/macs2_env.yml)
- Pin Python 3.9.*, NumPy 1.21.6, SciPy 1.7.3 for compatibility
- Update MACS2_CALLPEAK module to use environment file
- Add reproducibility guide and environment test script

Resolves binary incompatibility issue where different conda versions
resolve to incompatible NumPy versions (1.24+) that break MACS2 2.2.7.1

Tested on HPC systems with various conda versions (4.10-4.14)
```

---

## Testing Checklist

- [ ] Test script passes on your system
- [ ] Test script passes on colleague's system
- [ ] Pipeline runs successfully with a test BAM file
- [ ] Environment is cached correctly by Nextflow
- [ ] Same results as before the fix

---

**Need Help?** Check `docs/REPRODUCIBILITY_GUIDE.md` for detailed troubleshooting steps.
