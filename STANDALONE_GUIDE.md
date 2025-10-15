# PeakPrime Explorer - Complete Guide

A fast Shiny application for exploring PeakPrime pipeline outputs with comprehensive visualization capabilities.

## � **Quick Start**

### **Step 1: Preprocessing** (Required once per results directory)

```bash
# From your PeakFindR directory
Rscript preprocess_for_standalone.R /path/to/results/directory

# Examples:
Rscript preprocess_for_standalone.R results/test
```

### **Step 2: Launch App** 

```bash
# Run the standalone app (after preprocessing)
Rscript app_standalone.R /path/to/results/directory

# Examples:
Rscript app_standalone.R results/test
```

### **R Console Usage:**
```r
# Preprocess data (flexible GTF handling)
source("preprocess_for_standalone.R")
preprocess_peakprime_standalone("results/FirstOrder_auto_300distThr")  # Auto-detect GTF
# OR specify custom GTF:
preprocess_peakprime_standalone("results/your_analysis", "/path/to/custom.gtf")

# Launch app
source("app_standalone.R")  # Enter directory when prompted
runApp(peakprime_app)      # Launch the interface
```

## 🎯 **Key Features**

### ⚡ **Fast Performance**
- **Ultra-fast loading**: Core gene and coverage data loaded from optimized RDS files (2-3 seconds)
- **Sub-second gene switching**: Instant visualization updates for fluid exploration
- **Memory efficient**: Smart data structures minimize RAM usage
- **Scalable**: Handles thousands of genes without performance degradation

### 🔬 **Comprehensive Peak Analysis**
- **Multi-peak visualization**: Display all detected MACS2 peaks per gene, not just the selected best peak
- **Peak ranking system**: Visual hierarchy showing peak confidence and selection logic  
- **Conditional display**: Toggle between single-peak and multi-peak modes
- **Scientific insight**: Understand peak detection diversity and selection criteria

### 🧬 **Advanced Gene Visualization**
- **Multi-isoform support**: All transcript variants displayed with visible intron lines
- **Strand-aware primer design**: Single primer arrows pointing toward 3' end
- **Peak highlighting**: Selected peak prominently highlighted across all tracks
- **Coverage integration**: Seamless alignment of coverage data with gene structure

### 🎛️ **Interactive Controls**
- **Gene selection**: Fast dropdown with all available genes
- **Peak display modes**: Switch between selected peak only and all detected peaks
- **Peak count control**: Adjustable slider for maximum peaks displayed (1-10)
- **Primer visualization**: Toggle primer location display on/off
- **Y-axis scaling**: Percentage, absolute, or log10 coverage scales

## 📁 **File Requirements**

### **Raw Nextflow outputs** (automatically discovered):
- `peaks_qc_summary.tsv` - Gene selection results and QC metrics
- `*.bw` - BigWig coverage files for visualization  
- `*_peaks.narrowPeak` - MACS2 peak calls for multi-peak analysis
- `selected_peaks.tsv` - Final selected peak regions (optional)

### **GTF annotation** (flexible options):
- **Auto-detected**: Default UGent HPC path (`/data/gent/vo/000/gvo00027/resources/...`)
- **Local discovery**: Searches for `*.gtf` files in results directory
- **Custom path**: Specify during preprocessing or via interactive prompt
- **Optional**: Can skip GTF processing (limits gene structure visualization)

### **Generated files** (created by preprocessing):
- `qc_data.rds` - Enhanced QC data with coordinate parsing  
- `gtf_data.rds` - Processed gene structure annotations
- `coverage_index.rds` - Optimized coverage data for all genes
- `data_manifest.rds` - Processing metadata and performance info

## 🎮 **User Interface**

### **Sidebar Controls**
1. **Gene Selection**
   - Fast dropdown menu with all available genes
   - Auto-populated from preprocessed gene data

2. **Peak Display Options**
   - ☐ **"Show all detected peaks"**: Toggle multi-peak visualization
   - 🎚️ **"Max peaks to display"**: Slider control (1-10 peaks)
   - Only visible when multi-peak mode is enabled

3. **Visualization Controls**
   - ☐ **"Show primers"**: Toggle primer arrow display
   - 📊 **"Y-axis scale"**: Percentage/Absolute/Log10 options

4. **Instructions Panel**
   - Built-in usage guide and feature explanations

### **Main Panel**
1. **Coverage Plot** (Top)
   - High-resolution coverage data visualization
   - Selected peak highlighting (yellow highlight)
   - Y-axis scaling options for different data ranges

2. **Gene Structure Plot** (Bottom)
   - Multi-isoform transcript display with visible intron lines
   - Primer track with strand-aware single arrows
   - Peak track showing detected peaks with visual hierarchy

## 🔬 **Multi-Peak Visualization System**

When **"Show all detected peaks"** is enabled:

**Selected Peak (Highlighted)**:
- Yellow highlighting across coverage and gene structure plots
- Dark red rectangle on peak track with prominent borders
- Primer always points to selected peak center

**Additional Detected Peaks**:
- Light gray rectangles on peak track only (no cross-plot highlighting)
- Muted borders to distinguish from selected peak
- No peak numbering labels for clean visualization

**Peak Ranking**:
- Peaks ordered by MACS2 score (highest = rank #1)
- Visual hierarchy through color and opacity
- Configurable display count (1-10 peaks maximum)

## 🔧 **Troubleshooting**

### 🚨 **Common Issues**

1. **"GTF data: ❌ Missing" or "Coverage data: ❌ Missing"**
   
   This means preprocessing didn't complete successfully. **Step-by-step fix:**
   
   ```r
   # Step 1: Check what RDS files exist
   results_dir <- "results/testboth"  # Adjust your path
   list.files(results_dir, pattern = "\\.rds$")
   # Expected: qc_data.rds, data_manifest.rds, gtf_data.rds, coverage_index.rds
   
   # Step 2: Load required packages FIRST
   library(data.table)
   library(rtracklayer) 
   library(GenomicRanges)
   library(IRanges)
   
   # Step 3: Re-run preprocessing with packages loaded and ensure the GTf file exists
   source("preprocess_for_standalone.R")
   preprocess_peakprime_standalone(results_dir, "/path/to/your.gtf")
   
   # Step 4: Verify files were created
   list.files(results_dir, pattern = "\\.rds$")
   sapply(list.files(results_dir, pattern = "\\.rds$", full.names = TRUE), 
          function(f) paste(basename(f), ":", file.info(f)$size, "bytes"))
   ```

2. **"Error: could not find function 'fread'"**
   ```r
   # Install required packages
   install.packages("data.table")
   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
   BiocManager::install(c("rtracklayer", "GenomicRanges", "IRanges"))
   ```

3. **"Error: object 'peakprime_app' not found"**
   ```r
   # Re-source the app to create the object
   source("app_standalone.R")  # Enter directory when prompted
   runApp(peakprime_app)       # Then launch
   ```

### **Performance Optimization**
- **Large datasets**: Consider filtering to top genes before preprocessing
- **Memory constraints**: Reduce `max_coverage_points` in preprocessing (default: 1000)
- **Network storage**: Copy RDS files locally for faster access
- **HPC usage**: Use RStudio Server for better performance than SSH+terminal

## 🌐 **Running on UGent HPC**

### **RStudio Server (Recommended)**
For UGent users, you can run the app directly on the HPC web portal:

📖 **See detailed guide**: [HPC_RSTUDIO_GUIDE.md](HPC_RSTUDIO_GUIDE.md)

**Quick summary:**
1. 🌐 Login: https://login.hpc.ugent.be/
2. 🖥️ Request RStudio Server session (8-16GB memory, 4-8 hours)
3. 📁 Navigate to PeakFindR directory in RStudio
4. 🚀 Run preprocessing + app launch from R console

**Benefits:**
- ✅ No local R setup needed
- ✅ High-performance HPC resources  
- ✅ Web-based access from anywhere
- ✅ Direct access to HPC-stored results

## 📊 **Performance Characteristics**

### **Loading Times**
- **Initial startup**: 2-5 seconds (loads all preprocessed data)
- **Gene switching**: <0.1 seconds (instant)
- **Multi-peak mode**: +0.2 seconds (narrowPeak parsing overhead)

### **Memory Usage**
- **Base app**: ~100-500 MB (depends on dataset size)
- **Per additional peak**: ~1-5 MB (depends on peak density)
- **Recommended RAM**: 4-8 GB for large datasets

### **Scalability**
- **Tested with**: 10,000+ genes, 50,000+ peaks
- **Optimal performance**: Up to 1,000 genes per class
- **Peak display limit**: 10 peaks maximum (user-configurable)

## 🎯 **Integration with PeakPrime Pipeline**

### **Complete Workflow**
1. **Run PeakPrime pipeline** → Generate TSV/BED/narrowPeak outputs
2. **Preprocessing step** → Convert to RDS + organize for app
3. **Launch explorer** → Interactive analysis and validation
4. **Export findings** → Publication-ready plots and peak coordinates

---

**📖 Additional guides:**
- `HPC_RSTUDIO_GUIDE.md` - UGent HPC RStudio Server setup  
- `README.md` - Complete PeakPrime pipeline documentation