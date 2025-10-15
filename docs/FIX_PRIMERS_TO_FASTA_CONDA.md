# PRIMERS_TO_FASTA Conda Environment Fix

**Date:** October 14, 2025  
**Issue:** Process `primer_design:PRIMERS_TO_FASTA` failed during conda environment creation

---

## 🐛 **Problem Description**

The pipeline failed with this error:

```
Failed to create Conda environment
  command: conda create ... r-base>=4.3 r-essentials bioconductor-biostrings r-optparse
  status : 143
  message:
```

### **What Happened:**

1. **Exit code 143** = Process was **killed/terminated** (SIGTERM signal)
2. **26 minutes hang** - Environment creation started at 09:04:42, killed at 09:30:57
3. **Root cause:** `r-essentials` is a **metapackage with 100+ R packages**
4. **Result:** Conda dependency solver got stuck trying to resolve all dependencies

---

## ✅ **Solution Applied**

### **Changed:** `modules/PRIMERS_TO_FASTA.nf`

**Before:**
```groovy
conda "r-base>=4.3 r-essentials bioconductor-biostrings r-optparse"
```

**After:**
```groovy
conda "${projectDir}/env/primers_to_fasta_env.yml"
```

### **Created:** `env/primers_to_fasta_env.yml`

A **minimal environment file** with only required packages:
- ✅ `r-base=4.3.*` (specific version)
- ✅ `bioconductor-biostrings` (only what's needed)
- ✅ `r-optparse` (for command-line parsing)
- ✅ Essential Bioconductor dependencies only
- ❌ Removed `r-essentials` (bloated metapackage)

---

## 🎯 **Why This Fixes It**

### **Problem with r-essentials:**
- **100+ packages** including dplyr, ggplot2, tidyr, etc.
- **Complex dependency tree** - conda solver can hang for 30+ minutes
- **Not needed** - this script only uses Biostrings

### **Benefits of the fix:**
- ✅ **Faster** - ~2-3 minutes instead of 26+ minutes
- ✅ **Reliable** - simpler dependency resolution
- ✅ **Lighter** - smaller environment size
- ✅ **Specific versions** - better reproducibility

---

## 🚀 **Next Steps**

### **1. Clean conda cache** (recommended):
```bash
bash bin/clean_conda_cache.sh --force
```

### **2. Re-run the pipeline**:
```bash
nextflow run main.nf \
  --bam your_file.bam \
  --gtf annotations.gtf \
  --genes gene_list.txt \
  [... other parameters ...]
```

### **3. Monitor environment creation**:
The PRIMERS_TO_FASTA environment should now create in **~2-3 minutes** instead of hanging.

---

## 🔍 **Technical Details**

### **Exit Code 143 Explained:**
- **143 = 128 + 15**
- Signal 15 = SIGTERM (graceful termination request)
- Likely triggered by:
  - Nextflow timeout
  - HPC scheduler time limit
  - System resource manager
  - Manual cancellation

### **Why Conda Hangs:**
1. **Dependency resolution** is NP-complete problem
2. **r-essentials** has circular dependencies
3. **Multiple channels** (conda-forge, bioconda, defaults) increases search space
4. **Version constraints** multiply exponentially with package count

### **Log Timeline:**
```
09:04:42 - Environment creation started
09:05:03 - Pipeline waiting (5 min mark)
09:10:03 - Still waiting (10 min mark)
09:15:03 - Still waiting (15 min mark)
09:20:03 - Still waiting (20 min mark)
09:25:03 - Still waiting (25 min mark)
09:30:57 - Process killed (26 min total)
```

---

## 🛡️ **Prevention**

### **Best Practices for Conda Environments:**

1. ✅ **Avoid metapackages** (r-essentials, anaconda, etc.)
2. ✅ **Use environment files** instead of inline specifications
3. ✅ **Pin major versions** (r-base=4.3.* not r-base>=4.3)
4. ✅ **Minimize dependencies** - only include what you actually use
5. ✅ **Use mamba** for faster solving (optional)

### **Alternative: Use Mamba**

In `nextflow.config`:
```groovy
conda.useMamba = true  // Change from false to true
```

Mamba is a **faster conda replacement** with better dependency resolution.

---

## 📝 **Other Modules to Check**

Similar pattern exists in other modules - might need similar fixes:

```bash
grep -r "r-essentials" modules/
```

If other modules use `r-essentials`, consider creating minimal environment files for them too.

---

## ✅ **Verification**

After applying the fix, you should see:

```
[Actor Thread] INFO  nextflow.conda.CondaCache - Creating env using conda: 
    /path/to/env/primers_to_fasta_env.yml 
    [cache /user/gent/.../conda_cache/env-XXXXXXXXXX]
    
[2-3 minutes later...]

[Actor Thread] DEBUG nextflow.conda.CondaCache - 'conda' create complete
[Task submitter] INFO  nextflow.Session - [XX/XXXXXX] Submitted process > 
    primer_design:PRIMERS_TO_FASTA (1)
```

**Success indicators:**
- ✅ Environment creation completes in **< 5 minutes**
- ✅ Process submits successfully
- ✅ No exit code 143 errors

---

## 🆘 **If Problems Persist**

### **Option 1: Try inline specification**

In `PRIMERS_TO_FASTA.nf`:
```groovy
conda "r-base=4.3.* bioconductor-biostrings r-optparse"
```

### **Option 2: Pre-create environment**

```bash
conda env create -f env/primers_to_fasta_env.yml -n test_primers_to_fasta
conda activate test_primers_to_fasta

# Test the script
Rscript bin/primers_to_fasta.R --help
```

### **Option 3: Use mamba**

```bash
conda install -n base mamba
# Then set conda.useMamba = true in nextflow.config
```

---

**Questions?** Check the reproducibility guide or create an issue with the full error log.
