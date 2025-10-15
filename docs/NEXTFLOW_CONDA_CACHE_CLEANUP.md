# Cleaning Nextflow's Custom Conda Cache

## The Problem

When Nextflow uses custom conda cache directories (configured via environment variables or `nextflow.config`), the standard `conda clean --all` command **will not** clean them. This is because:

1. `conda clean` only operates on conda's default cache locations
2. Nextflow uses custom directories to avoid conflicts and manage environments separately
3. Your pipeline uses these custom paths:
   - `NXF_CONDA_CACHEDIR="/user/gent/446/vsc44685/ScratchVO_dir/conda_cache"`
   - `CONDA_PKGS_DIRS="/user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs"`

---

## Solutions

### ✅ **Solution 1: Use the New Cleanup Script (Recommended)**

We've created a dedicated cleanup script for this:

```bash
# Interactive cleanup (asks for confirmation)
bash bin/clean_conda_cache.sh

# OR force cleanup without confirmation
bash bin/clean_conda_cache.sh --force
```

**What it does:**
- ✅ Cleans standard conda cache (`conda clean --all`)
- ✅ Cleans Nextflow conda cache (`$NXF_CONDA_CACHEDIR`)
- ✅ Cleans conda packages cache (`$CONDA_PKGS_DIRS`)
- ✅ Shows cache sizes before cleaning
- ✅ Offers to remove test environments
- ✅ Safe with confirmation prompts (unless `--force`)

---

### ✅ **Solution 2: Manual Cleanup**

If you prefer to do it manually:

```bash
# 1. Clean standard conda cache
conda clean --all -y

# 2. Clean Nextflow's conda cache
rm -rf /user/gent/446/vsc44685/ScratchVO_dir/conda_cache/*

# 3. Clean conda packages cache
rm -rf /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs/*
```

---

### ✅ **Solution 3: Enhanced Test Script**

The test script now supports cache cleaning:

```bash
# Test environment with cache cleaning
bash bin/test_macs2_environment.sh --clean
```

This will:
1. Clean all conda caches (including Nextflow's custom ones)
2. Create a fresh test environment
3. Verify MACS2 compatibility

---

## Understanding Nextflow's Conda Cache

### **How Nextflow Uses Conda**

When you run a Nextflow pipeline with conda enabled:

1. **First run**: Nextflow creates conda environments in `$NXF_CONDA_CACHEDIR`
2. **Subsequent runs**: Nextflow reuses cached environments (much faster!)
3. **Cache location**: By default in `.nextflow/` but your setup uses custom paths

### **Why Custom Cache Directories?**

From your `clean_caches_force.sh`, you're using custom cache directories likely because:
- **Shared HPC system**: Default home directory has limited space
- **Better performance**: ScratchVO_dir is usually on faster storage
- **Avoid conflicts**: Separate caches for different projects

### **When to Clean These Caches**

Clean the Nextflow conda cache when:
- ❌ Getting conda/environment errors (like the NumPy issue)
- 🔄 After updating environment specifications
- 🐛 Pipeline fails with "environment activation" errors
- 💾 Running out of disk space
- 🧪 Testing environment reproducibility

**Don't clean if:**
- ✅ Pipeline is working fine (cached environments save time!)
- ⚡ You want fast subsequent runs

---

## Verification

After cleaning, verify the caches are empty:

```bash
# Check sizes
du -sh /user/gent/446/vsc44685/ScratchVO_dir/conda_cache
du -sh /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs

# Check contents (should be empty or minimal)
ls -la /user/gent/446/vsc44685/ScratchVO_dir/conda_cache/
ls -la /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs/
```

---

## Complete Workflow for Your Colleague

**Step-by-step troubleshooting process:**

```bash
# 1. Pull latest changes with the fix
git pull

# 2. Clean all conda caches
bash bin/clean_conda_cache.sh --force

# 3. Test environment creation
bash bin/test_macs2_environment.sh

# 4. If test passes, run your pipeline
nextflow run main.nf \
  --bam your_file.bam \
  --gtf annotations.gtf \
  --genes genes.txt \
  --outdir results/ \
  -profile local
```

---

## Advanced: Configure Custom Cache Locations

If your colleague needs different cache paths, they can:

### **Option A: Environment Variables**
```bash
export NXF_CONDA_CACHEDIR="/your/custom/path/conda_cache"
export CONDA_PKGS_DIRS="/your/custom/path/conda_pkgs"
nextflow run main.nf [parameters]
```

### **Option B: Nextflow Config**
Add to `nextflow.config`:
```groovy
conda.cacheDir = "/your/custom/path/conda_cache"
```

### **Option C: Command Line**
```bash
nextflow run main.nf -with-conda -conda-cache-dir /your/custom/path/conda_cache
```

---

## Troubleshooting

### **Issue: Permission Denied**
```bash
# Fix permissions
chmod -R u+w /user/gent/446/vsc44685/ScratchVO_dir/conda_cache
chmod -R u+w /user/gent/446/vsc44685/ScratchVO_dir/conda_pkgs

# Then clean
bash bin/clean_conda_cache.sh --force
```

### **Issue: Disk Space Still Full**
```bash
# Check what's using space
du -sh /user/gent/446/vsc44685/ScratchVO_dir/*

# Also clean Nextflow work directory
rm -rf work/

# Clean singularity/apptainer caches if using containers
bash bin/clean_caches_force.sh
```

### **Issue: Conda Command Not Found**
```bash
# Initialize conda for your shell
conda init bash
source ~/.bashrc

# Or use absolute path
/path/to/conda/bin/conda clean --all -y
```

---

## Summary: Key Points

1. **Standard `conda clean` doesn't touch Nextflow's custom caches**
2. **Use the new `bin/clean_conda_cache.sh` script** for comprehensive cleaning
3. **Clean caches when troubleshooting environment issues**
4. **Keep caches if pipeline is working** (saves time on subsequent runs)
5. **After cleaning, test with `bin/test_macs2_environment.sh`**

---

## Related Files

- `bin/clean_conda_cache.sh` - Interactive conda cache cleanup (NEW)
- `bin/clean_caches_force.sh` - Clean ALL caches including Singularity (EXISTING)
- `bin/test_macs2_environment.sh` - Test MACS2 environment (supports `--clean` flag)
- `env/macs2_env.yml` - Locked environment specification
- `docs/REPRODUCIBILITY_GUIDE.md` - Full troubleshooting guide

---

**Questions?** Check `docs/REPRODUCIBILITY_GUIDE.md` for more details!
