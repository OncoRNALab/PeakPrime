# Running PeakPrime Explorer App on UGent HPC RStudio Server

## 🎯 **Overview**

This guide shows you how to run the PeakPrime standalone Shiny app using the UGent HPC web portal's RStudio server. This approach gives you a full RStudio environment with web-based access - perfect for interactive data exploration without local R setup requirements.

---

## 🚀 **Step-by-Step Instructions**

### **Step 1: Access UGent HPC Web Portal**

1. **Navigate to the HPC login portal:**
   - 🌐 **URL**: https://login.hpc.ugent.be/
   - 📖 **Documentation**: https://docs.hpc.ugent.be/macOS/web_portal/

2. **Log in and authorize:**
   - Use your UGent credentials
   - Complete any required 2FA authentication
   - Accept authorization prompts

### **Step 2: Request Interactive RStudio Session**

1. **Go to Interactive Sessions:**
   - In the web portal, navigate to `Interactive Apps` → `RStudio Server`

2. **Configure your session:**
   ```
   Cluster: Choose appropriate cluster (e.g., doduo, donphan)
   Time: 4-8 hours (recommended for data exploration)
   Cores: 2-4 cores (sufficient for Shiny apps)
   Memory: 8-16 GB (depending on dataset size)
   R Version: R/4.3.0 or newer
   ```

3. **Submit and wait:**
   - Click `Launch`
   - Wait for session to start (usually 1-5 minutes)
   - Click `Connect to RStudio Server` when ready

### **Step 3: Set Up Your Environment**

1. **Navigate to your PeakFindR directory:**
   ```r
   # In RStudio console
   setwd("/user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR")
   getwd()  # Verify location
   ```

2. **Install required packages (if needed):**
   ```r
   # Check if packages are available
   required_packages <- c("shiny", "DT", "ggplot2", "GenomicRanges", 
                         "IRanges", "data.table", "grid", "rtracklayer")
   
   # Install any missing packages
   missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
   if (length(missing) > 0) {
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
     }
     BiocManager::install(missing)
   }
   ```

### **Step 4: Preprocess Your Data (If Needed)**

1. **Check if preprocessing is needed:**
   ```r
   # Check for RDS files in your results directory
   results_dir <- "results/testboth"  # Adjust path as needed
   required_files <- c("qc_data.rds", "data_manifest.rds")
   file.exists(file.path(results_dir, required_files))
   ```

2. **Run preprocessing if required:**
   ```r
   # Method 1: Auto-detect GTF file
   source("preprocess_for_standalone.R")
   preprocess_peakprime_standalone("results/testboth")
   
   # Method 2: Specify custom GTF path
   preprocess_peakprime_standalone("results/testboth", "/path/to/your.gtf")
   
   # Method 3: Skip GTF processing entirely
   preprocess_peakprime_standalone("results/testboth", gtf_path = NA)
   ```

### **Step 5: Troubleshoot Missing Data (If Needed)**

If you see "GTF data: ❌ Missing" or "Coverage data: ❌ Missing", check the following:

1. **Verify preprocessed files exist:**
   ```r
   # Check what files are actually in your results directory
   results_dir <- "results/testboth"
   list.files(results_dir, pattern = "\\.rds$")
   
   # Should show: qc_data.rds, data_manifest.rds, gtf_data.rds, coverage_index.rds
   ```

2. **Check file sizes (troubleshoot corruption):**
   ```r
   # Check file info
   results_dir <- "results/testboth"
   rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
   sapply(rds_files, function(f) file.info(f)$size)
   ```

3. **Re-run preprocessing if files are missing:**
   ```r
   # Re-run preprocessing with verbose output
   source("preprocess_for_standalone.R")
   preprocess_peakprime_standalone("results/testboth")
   ```

### **Step 6: Launch the Standalone App**

1. **Method 1: Two-step approach (Recommended for RStudio)**
   ```r
   # Step 1: Load the app (loads data and creates app object)
   source("app_standalone.R")
   # When prompted, enter: results/testboth
   
   # Step 2: Launch the app
   runApp(peakprime_app)
   ```

2. **Method 2: Auto-launch**
   ```r
   # Set launch flag and source to auto-launch
   .standalone_launch <- TRUE
   source("app_standalone.R")
   # When prompted, enter: results/testboth
   ```

3. **Method 3: Command line execution**
   ```r
   # Launch using Rscript (runs automatically)
   system("Rscript app_standalone.R results/testboth", wait = FALSE)
   ```

3. **The app will start and display:**
   ```
   📁 Using results directory: /path/to/results/testboth
   ✅ Found all required preprocessed data files
   🚀 Loading preprocessed data...
   ⚡ Loaded XXX genes in X.X secs
   
   Listening on http://127.0.0.1:XXXX
   ```

### **Step 7: Access Your Shiny App**

1. **Note the app URL:**
   - The app will show something like: `Listening on http://127.0.0.1:3838`
   - This is your local Shiny server address

2. **Open in RStudio viewer or browser:**
   - **In RStudio**: The app may automatically open in the Viewer pane
   - **In browser**: Copy the URL and open in a new browser tab within your RStudio session
   - **Alternative**: Use RStudio's "Show in new window" option

---

## 🔧 **Troubleshooting**

### **Common Issues & Solutions**

1. **"GTF data: ❌ Missing" or "Coverage data: ❌ Missing"**
   ```r
   # Check what preprocessing files exist
   results_dir <- "results/testboth"
   list.files(results_dir, pattern = "\\.rds$")
   
   # If missing files, re-run preprocessing
   source("preprocess_for_standalone.R") 
   preprocess_peakprime_standalone(results_dir)  # Will auto-detect or prompt for GTF
   
   # Or specify GTF path explicitly
   preprocess_peakprime_standalone(results_dir, "/path/to/your.gtf")
   
   # Check file sizes (should be > 0 bytes)
   rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
   sapply(rds_files, function(f) paste(basename(f), ":", file.info(f)$size, "bytes"))
   ```

2. **"Cannot find required packages"**
   ```r
   # Install Bioconductor packages
   if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
   BiocManager::install(c("GenomicRanges", "IRanges", "rtracklayer"))
   ```

3. **"Directory not found"**
   ```r
   # Check your current working directory
   getwd()
   # List available results directories
   list.dirs("results", recursive = FALSE)
   ```

4. **"App won't load in browser"**
   - Try refreshing the browser
   - Check the R console for the correct URL
   - Use RStudio's Viewer pane instead
   - Restart the R session if needed

5. **"Memory issues with large datasets"**
   - Request more memory in your HPC session (16-32 GB)
   - Reduce `max_coverage_points` in preprocessing
   - Filter to fewer genes before preprocessing

6. **"Error: object 'peakprime_app' not found"**
   ```r
   # Re-source the app to create the object
   source("app_standalone.R")
   # Enter directory when prompted, then:
   runApp(peakprime_app)
   ```

### **Performance Tips**

1. **For large datasets:**
   ```r
   # Preprocess with reduced coverage points for speed
   preprocess_peakprime_standalone("results/testboth", max_coverage_points = 500)
   ```

2. **For multiple results directories:**
   ```r
   # Preprocess multiple directories
   dirs <- c("results/Class1", "results/Class2", "results/testboth")
   for (dir in dirs) {
     if (dir.exists(dir)) {
       cat("Processing", dir, "...\n")
       preprocess_peakprime_standalone(dir)
     }
   }
   ```

---

## 💡 **Advanced Usage**

### **Running Multiple Apps Simultaneously**

```r
# Launch apps on different ports for comparison
library(shiny)

# App 1 (port 3838)
runApp(list(
  ui = source("app_standalone.R", local = TRUE)$value$ui,
  server = source("app_standalone.R", local = TRUE)$value$server
), port = 3838, launch.browser = FALSE)

# App 2 (port 3839) - different results
# Modify app_standalone.R to point to different results directory
runApp("app_standalone.R", port = 3839, launch.browser = FALSE)
```

### **Batch Processing Multiple Datasets**

```r
# Process and launch apps for multiple analyses
results_dirs <- c("results/Class1_W20_pad100", "results/Class2_W30_pad100", "results/testboth")

for (dir in results_dirs) {
  cat("=== Processing", dir, "===\n")
  
  # Preprocess if needed
  if (!file.exists(file.path(dir, "qc_data.rds"))) {
    preprocess_peakprime_standalone(dir)
  }
  
  # Launch app (will prompt for manual interaction)
  cat("Ready to launch app for", dir, "\n")
  cat("Run: source('app_standalone.R') and enter directory:", dir, "\n\n")
}
```

---

## ⚡ **Quick Reference**

### **Essential Commands**
```r
# Set working directory
setwd("/user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR")

# Preprocess data (if needed)
source("preprocess_for_standalone.R")
preprocess_peakprime_standalone("results/your_results_dir")

# Launch app (recommended method)
.standalone_launch <- TRUE
source("app_standalone.R")
# Enter: results/your_results_dir

# Alternative: Load then launch manually
# source("app_standalone.R")  # Enter directory when prompted
# runApp(peakprime_app)      # Launch manually
```

### **File Locations**
- **App files**: `/user/gent/446/vsc44685/ScratchVO_dir/OncoRNA_peakprime/Primer_PeakFindR/`
- **Results**: `results/` subdirectories
- **Preprocessed data**: `*.rds` files in each results directory

### **Session Requirements**
- **Time**: 4-8 hours recommended
- **Memory**: 8-16 GB (more for large datasets)  
- **Cores**: 2-4 cores sufficient
- **R Version**: R/4.3.0 or newer

---

## 🎉 **Benefits of HPC RStudio Approach**

✅ **No local R installation needed** - Everything runs on HPC  
✅ **Web-based interface** - Access from any browser  
✅ **High-performance computing** - More memory and CPU than laptops  
✅ **Persistent sessions** - Can run long analyses  
✅ **Data locality** - Direct access to HPC-stored results  
✅ **Collaborative** - Share session URLs with colleagues  

This setup gives you the full power of the PeakPrime standalone app with the convenience and performance of the UGent HPC infrastructure! 🚀