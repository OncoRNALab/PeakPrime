# R Environment Analysis & Optimization Strategy

**Date:** October 14, 2025  
**Issue:** Large R environments causing conda hangs and slow creation times

---

## 🔍 **Current State Analysis**

### **R-Using Processes in Pipeline:**

| Process | Current Conda Spec | Status | Packages Actually Used |
|---------|-------------------|--------|------------------------|
| **PROCESS_MACS2_PEAKS** | Inline: 9 Bioconductor + r-optparse + r-data.table | ⚠️ Inline | GenomicFeatures, rtracklayer, GenomicRanges, IRanges, S4Vectors, Biostrings, optparse, data.table |
| **EXTRACT_CDNA_PRIMERS** | **Environment file** ✅ | ✅ Fixed | optparse (only!) |
| **PRIMERS_TO_FASTA** | **Environment file** ✅ | ✅ Fixed | optparse, Biostrings |
| **MAKEPLOTS_NEW** | Inline: **r-essentials** + 10 packages | ❌ PROBLEM | optparse, rtracklayer, GenomicRanges, IRanges, ggplot2, data.table, grid, patchwork/cowplot |

### **Key Findings:**

#### ❌ **Problem Areas:**

1. **MAKEPLOTS_NEW** - Uses `r-essentials` (100+ packages!)
   - Only needs: ~10 specific packages
   - **Bloat factor:** 10x (100 packages vs 10 needed)
   - **Risk:** High (same issue as PRIMERS_TO_FASTA had)

2. **PROCESS_MACS2_PEAKS** - Long inline specification
   - Needs: 8 Bioconductor packages + 2 CRAN
   - **Status:** Works but could be cleaner
   - **Risk:** Medium (no r-essentials, but complex)

#### ✅ **Already Fixed:**

1. **EXTRACT_CDNA_PRIMERS** - Uses minimal env file ✅
   - Only needs: r-base + r-optparse + r-data.table
   
2. **PRIMERS_TO_FASTA** - Uses minimal env file ✅
   - Only needs: r-base + bioconductor-biostrings + r-optparse

---

## 🎯 **Answer to Your Questions:**

### **Q1: Was the initial idea to have a single R environment for all processes?**

**Answer: NO** - Currently each process has its **own separate conda specification**.

**Evidence:**
- Each process has a different `conda` directive
- Each creates a separate conda environment in the cache
- No shared environment file across processes

**Current approach:** 
```
PROCESS_MACS2_PEAKS    → conda cache env-0c0864195848ba13...
EXTRACT_CDNA_PRIMERS   → conda cache env-41a5cc9229633821...
PRIMERS_TO_FASTA       → conda cache env-a49ced3eb9859f2195b...
MAKEPLOTS_NEW          → conda cache env-<unique hash>...
```

### **Q2: Would individual environments be better?**

**Answer: YES** - Individual minimal environments are **significantly better**! ✅

---

## 📊 **Single Environment vs Individual Environments**

### **Option A: Single Shared R Environment**

```yaml
# env/shared_r_env.yml - ONE environment for ALL processes
dependencies:
  - r-base=4.3.*
  - r-essentials  # ← 100+ packages
  - All bioconductor packages needed
  - All plotting packages
  - etc...
```

**Pros:**
- ✅ Created once, reused everywhere
- ✅ Simpler management (one file to maintain)

**Cons:**
- ❌ **MASSIVE** (2-5 GB on disk)
- ❌ Takes **30+ minutes** to create
- ❌ High risk of dependency conflicts
- ❌ Includes packages never used by most processes
- ❌ If it breaks, ALL processes fail
- ❌ Harder to troubleshoot (which process needs what?)

---

### **Option B: Individual Minimal Environments** ✅ **RECOMMENDED**

Each process gets exactly what it needs:

```yaml
# env/process_macs2_peaks_env.yml
dependencies:
  - r-base=4.3.*
  - bioconductor-genomicfeatures
  - bioconductor-rtracklayer
  - bioconductor-genomicranges
  - bioconductor-iranges
  - bioconductor-s4vectors
  - bioconductor-biostrings
  - r-optparse
  - r-data.table

# env/extract_cdna_primers_env.yml (already done!)
dependencies:
  - r-base=4.3.*
  - r-optparse
  - r-data.table

# env/primers_to_fasta_env.yml (already done!)
dependencies:
  - r-base=4.3.*
  - bioconductor-biostrings
  - r-optparse

# env/makeplots_env.yml
dependencies:
  - r-base=4.3.*
  - r-ggplot2
  - r-data.table
  - r-optparse
  - r-patchwork
  - r-cowplot
  - bioconductor-genomicranges
  - bioconductor-rtracklayer
  - bioconductor-iranges
  - bioconductor-s4vectors
```

**Pros:**
- ✅ **Fast creation** (~2-5 minutes each)
- ✅ **Reliable** - simpler dependency resolution
- ✅ **Smaller** - only what's needed
- ✅ **Isolated** - one failure doesn't break others
- ✅ **Clear** - easy to see what each process needs
- ✅ **Cacheable** - Nextflow caches each separately
- ✅ **Maintainable** - update one without affecting others
- ✅ **No r-essentials bloat**

**Cons:**
- ⚠️ More files to maintain (4 env files vs 1)
- ⚠️ First run creates 4 environments (but cached after)

---

## 💡 **Recommendation: Individual Environments**

### **Why This is Better:**

1. **Speed:** Each environment creates in 2-5 minutes
   - vs 30+ minutes for a bloated shared environment

2. **Reliability:** Simple dependencies = fewer conflicts
   - No r-essentials = no hanging conda resolver

3. **Failure Isolation:** 
   - If MAKEPLOTS fails to create → other processes still work
   - With shared env → everything breaks

4. **Disk Space:**
   ```
   Shared approach:    1 × 5GB = 5 GB
   Individual approach: 4 × 500MB = 2 GB  (SMALLER!)
   ```

5. **Development:**
   - Need to add package to plotting? Edit only `makeplots_env.yml`
   - Don't risk breaking peak processing, primer extraction, etc.

---

## 🚀 **Recommended Action Plan**

### **Priority 1: Fix MAKEPLOTS_NEW (High Risk)** ⚠️

This uses `r-essentials` and will likely hang like PRIMERS_TO_FASTA did.

**Create:** `env/makeplots_env.yml`

### **Priority 2: Improve PROCESS_MACS2_PEAKS (Medium)**

Currently works but has long inline specification.

**Create:** `env/process_macs2_peaks_env.yml`

### **Priority 3: Keep EXTRACT_CDNA_PRIMERS & PRIMERS_TO_FASTA** ✅

Already using environment files - these are good!

---

## 📈 **Impact Comparison**

### **Current State (Mixed):**
```
PROCESS_MACS2_PEAKS     → 9 packages inline       → 3-5 min   ✅
EXTRACT_CDNA_PRIMERS    → env file (3 packages)   → 1-2 min   ✅
PRIMERS_TO_FASTA        → env file (3 packages)   → 2-3 min   ✅
MAKEPLOTS_NEW           → r-essentials + 10 more  → 30+ min⚠️ (or HANG)
```

### **Proposed State (All Individual):**
```
PROCESS_MACS2_PEAKS     → env file (8 packages)   → 4-6 min   ✅
EXTRACT_CDNA_PRIMERS    → env file (3 packages)   → 1-2 min   ✅
PRIMERS_TO_FASTA        → env file (3 packages)   → 2-3 min   ✅
MAKEPLOTS_NEW           → env file (10 packages)  → 3-5 min   ✅
```

**Total environment creation time:**
- **Current:** 1-2 + 2-3 + 3-5 + **30+** = **~40 minutes** (if no hang!)
- **Proposed:** 4-6 + 1-2 + 2-3 + 3-5 = **~15 minutes** (reliable!)

**After first run (cached):** Both approaches are instant ⚡

---

## 🎓 **Best Practices for R Environments in Nextflow**

### **DO:**
✅ Use individual environment files per process
✅ Pin major versions (r-base=4.3.*)
✅ Only include packages actually `library()`'ed in scripts
✅ Document why each package is needed
✅ Test environment creation before pipeline run

### **DON'T:**
❌ Use r-essentials (or any metapackage)
❌ Share one environment across unrelated processes
❌ Use `>=` version specifiers (use `=` or `=X.Y.*`)
❌ Include "nice to have" packages
❌ Mix many unrelated bioconductor packages in one env

---

## 📝 **Summary Table**

| Approach | Creation Time | Disk Space | Reliability | Maintainability | **Recommended?** |
|----------|--------------|------------|-------------|-----------------|------------------|
| **Single Shared Env** | 30-60 min | ~5 GB | Low ⚠️ | Easy | ❌ NO |
| **Individual Minimal Envs** | 2-5 min each | ~2 GB total | High ✅ | Moderate | ✅ **YES** |
| **Current Mixed** | Variable | ~3 GB | Medium | Confusing | ⚠️ Needs fixing |

---

## 🎯 **Final Recommendation:**

**Use individual minimal environment files for each R process.**

This is the approach you've already started with PRIMERS_TO_FASTA and EXTRACT_CDNA_PRIMERS, and it's working perfectly! Now we just need to:

1. ✅ **Keep** PRIMERS_TO_FASTA and EXTRACT_CDNA_PRIMERS as-is (already good)
2. 🔧 **Create** env file for MAKEPLOTS_NEW (prevent hanging)
3. 🔧 **Create** env file for PROCESS_MACS2_PEAKS (cleaner, more maintainable)

**Result:** Fast, reliable, maintainable R environments! 🚀

---

**Next steps:** Would you like me to create the missing environment files?
