# PeakPrime Explorer

A fast Shiny application for exploring PeakPrime pipeline outputs with comprehensive visualization capabilities.

## Overview

The PeakPrime Explorer combines **preprocessed RDS files** for lightning-fast loading with **direct narrowPeak file reading** for comprehensive peak analysis. This hybrid approach provides both speed and scientific depth for peak analysis workflows.

## Key Features

### ⚡ Fast Performance
- **Preprocessed data loading**: Core gene and coverage data loaded from optimized RDS files
- **Sub-second gene switching**: Instant visualization updates for fluid exploration
- **Memory efficient**: Smart data structures minimize RAM usage
- **Scalable**: Handles thousands of genes without performance degradation

### 🔬 Comprehensive Peak Analysis
- **Multi-peak visualization**: Display all detected MACS2 peaks per gene, not just the selected best peak
- **Peak ranking system**: Visual hierarchy showing peak confidence and selection logic  
- **Conditional display**: Toggle between single-peak and multi-peak modes
- **Scientific insight**: Understand peak detection diversity and selection criteria

### 🧬 Advanced Gene Visualization
- **Multi-isoform support**: All transcript variants displayed with visible intron lines
- **Strand-aware primer design**: Single primer arrows pointing toward 3' end
- **Peak highlighting**: Selected peak prominently highlighted across all tracks
- **Coverage integration**: Seamless alignment of coverage data with gene structure

### 🎛️ Interactive Controls
- **Gene selection**: Fast dropdown with all available genes
- **Peak display modes**: Switch between selected peak only and all detected peaks
- **Peak count control**: Adjustable slider for maximum peaks displayed (1-10)
- **Primer visualization**: Toggle primer location display on/off
- **Y-axis scaling**: Percentage, absolute, or log10 coverage scales

## Installation & Usage

### Prerequisites
```r
# Required R packages
install.packages(c("shiny", "ggplot2", "data.table", "DT"))

# Optional (for enhanced plot layouts)
install.packages(c("patchwork", "cowplot"))
```

### Integration with Peak Prime
The app requires preprocessing of PeakPrime outputs into RDS format:

```bash
# Example preprocessing command

```

```bash
# 1. Navigate to your PeakPrime results directory
cd /path/to/results/
# 2. Preprocess the ouputs of peakPrimer before launching the app
Rscript preprocess_for_hybrid.R --input /path/to/peakprime/results --output /path/to/ClassX_final/
# 3. Launch the explorer
R -e "shiny::runApp('app_hybrid.R')"
```

### File Requirements

**Core files** (required for basic functionality):
- `genes_data.rds` - Preprocessed gene annotations and metadata
- `coverage_data.rds` - Preprocessed coverage data for all genes  
- `qc_summary.rds` - Quality control metrics and peak information

**Peak analysis files** (required for multi-peak visualization):
- `*_peaks.narrowPeak` - MACS2 peak calls for comprehensive peak analysis
- Pattern examples: `Class1_peaks.narrowPeak`, `treatment_peaks.narrowPeak`

**Auto-detected files** (enhance functionality if present):
- `*.gtf` or `*.gff` - Gene annotation files for coordinate validation
- `primer_*.bed` - Primer target coordinates for enhanced primer display

## User Interface

### Sidebar Controls
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

### Main Panel
1. **Coverage Plot** (Top)
   - High-resolution coverage data visualization
   - Selected peak highlighting (yellow highlight)
   - Y-axis scaling options for different data ranges

2. **Gene Structure Plot** (Bottom)
   - Multi-isoform transcript display with visible intron lines
   - Primer track with strand-aware single arrows
   - Peak track showing detected peaks with visual hierarchy

## Advanced Features

### Multi-Peak Visualization System
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

### Hybrid Data Architecture
**Preprocessed Components**:
- Gene coordinates and metadata → `genes_data.rds`
- Coverage matrices → `coverage_data.rds` 
- QC metrics and selected peaks → `qc_summary.rds`

**Real-time Components**:
- narrowPeak file parsing for comprehensive peak lists
- Dynamic gene-to-peaks mapping
- On-demand peak coordinate retrieval

**Performance Benefits**:
- Core data loads in seconds regardless of dataset size
- Gene switching is instantaneous
- Peak analysis adds minimal overhead
- Memory usage scales linearly with displayed peaks

## Technical Implementation

### Core Functions

**Data Loading**:
```r
# Preprocessed data (loaded once at startup)
genes_data <- readRDS("genes_data.rds")
coverage_data <- readRDS("coverage_data.rds")
qc_summary <- readRDS("qc_summary.rds")

# Dynamic peak analysis
all_peaks_data <- load_all_peaks_data()  # From narrowPeak files
```

**Peak Analysis**:
```r
get_all_peaks <- function(gene_id, all_peaks_data, max_peaks = 5) {
  # Retrieve and rank all detected peaks for a gene
  # Returns top N peaks sorted by score
}
```

**Visualization Engine**:
```r
plot_gene_comprehensive <- function(gene_id, show_all_peaks = FALSE, max_peaks = 5, ...) {
  # Unified plotting function supporting both single and multi-peak modes
  # Uses ggplot2 with immediate evaluation to avoid lazy evaluation bugs
}
```

### Error Handling & Robustness
- **Missing file graceful degradation**: App functions with partial data
- **Peak coordinate validation**: Automatic bounds checking and correction
- **Memory management**: Smart data loading prevents memory overflow
- **ggplot2 evaluation fixes**: Uses `annotate()` to avoid lazy evaluation issues

## Performance Characteristics

### Loading Times
- **Initial startup**: 2-5 seconds (loads all preprocessed data)
- **Gene switching**: <0.1 seconds (instant)
- **Multi-peak mode**: +0.2 seconds (narrowPeak parsing overhead)

### Memory Usage
- **Base app**: ~100-500 MB (depends on dataset size)
- **Per additional peak**: ~1-5 MB (depends on peak density)
- **Recommended RAM**: 4-8 GB for large datasets

### Scalability
- **Tested with**: 10,000+ genes, 50,000+ peaks
- **Optimal performance**: Up to 1,000 genes per class
- **Peak display limit**: 10 peaks maximum (user-configurable)

## Troubleshooting

### Common Issues

1. **"No RDS files found"**
   - Ensure you're in a directory with preprocessed `*.rds` files
   - Check that preprocessing step completed successfully

2. **"No peaks detected for gene"**
   - Verify narrowPeak files exist and are readable
   - Check gene ID matches between RDS and narrowPeak data

3. **"App loading slowly"**
   - RDS files may be corrupted - try reprocessing
   - Check available RAM (large datasets require more memory)

4. **"Selected peak disappears in multi-peak mode"**
   - This was a previous bug - ensure you're using the latest version
   - Selected peak should always be visible and highlighted

### Performance Optimization
- **Large datasets**: Consider filtering to top genes before preprocessing
- **Memory constraints**: Reduce max peaks displayed or limit gene set
- **Network storage**: Copy RDS files locally for faster access

## Integration with PeakPrime Pipeline

### Preprocessing Requirements
The app requires preprocessing of PeakPrime outputs into RDS format:
```bash
# Example preprocessing command
Rscript preprocess_for_hybrid.R --input /path/to/peakprime/results --output /path/to/ClassX_final/
```

### File Compatibility
- **Input formats**: Standard PeakPrime TSV/BED/narrowPeak outputs  
- **Output formats**: Optimized RDS files + original narrowPeak preservation
- **Versioning**: Compatible with PeakPrime v2.0+ outputs

### Workflow Integration
1. **Run PeakPrime pipeline** → Generate TSV/BED/narrowPeak outputs
2. **Preprocessing step** → Convert to RDS + organize for hybrid app
3. **Launch hybrid explorer** → Interactive analysis and validation
4. **Export findings** → Publication-ready plots and peak coordinates

This hybrid approach provides the best of both worlds: the speed of preprocessed data with the scientific depth of comprehensive peak analysis, making it ideal for both exploratory analysis and detailed validation workflows.