# PeakPrime Documentation

Welcome to the PeakPrime documentation! This guide will help you find the information you need.

---

## 📚 Quick Navigation

### Getting Started
- **[Pipeline Overview](pipeline_steps.md)** - Understand the workflow and pipeline steps
- **[Reproducibility Guide](REPRODUCIBILITY_GUIDE.md)** - Setup, installation, and reproducibility
- **[Troubleshooting](troubleshooting/TROUBLESHOOTING.md)** - Common issues and solutions

### Core Documentation
- **[Scientific Manuscript](MANUSCRIPT.md)** - Detailed methodology and validation
- **[Technical Appendix](TECHNICAL_APPENDIX.md)** - File formats and technical specifications

---

## 🎯 Features

### Peak-Based Primer Design (Default Mode)
Design primers targeting RNA-seq coverage peaks:
- Automatic peak calling with MACS2
- Peak quality filtering and exonic overlap checks
- Primer3 optimization for cDNA applications

### Distance-Based Primer Design
Design primers at fixed distance from 3' end:
- **[Distance Mode Guide](features/DISTANCE_MODE.md)** - For 3' RNA-seq protocols (QuantSeq, etc.)
- Automatic MANE transcript fetching
- No BAM file required

### Advanced Features
- **[Isoform Optimization](features/ISOFORM_OPTIMIZATION.md)** - Maximize distinct isoforms targeted
- **[Peak Ranking](features/PEAK_RANKING.md)** - Understand how peaks are ranked and selected
- **[Plotting Optimization](features/PLOTTING_OPTIMIZATION.md)** - Fast parallel visualization (40× speedup)

### Interactive Visualization
- **[Standalone App Guide](STANDALONE_GUIDE.md)** - Interactive Shiny app for exploring results
- **[HPC RStudio Guide](HPC_RSTUDIO_GUIDE.md)** - Running the app on UGent HPC web portal

---

## 🚀 Quick Start Examples

### Example 1: Peak-Based Design with Plotting
```bash
nextflow run main.nf \
  --bam sample.bam \
  --gtf genome.gtf \
  --genes genes.txt \
  --makeplots \
  --outdir results \
  -profile conda
```

### Example 2: Distance-Based Design (3' RNA-seq)
```bash
nextflow run main.nf \
  --distance_mode \
  --genes genes.txt \
  --template_length 300 \
  --outdir results \
  -profile conda
```

### Example 3: Standalone Plotting
```bash
nextflow run main.nf \
  --makeplots \
  --genes genes.txt \
  --bw sample.bw \
  --gtf genome.gtf \
  --peaks_tsv results/peaks.tsv \
  --primer_targets_bed results/primers.bed \
  --qc_tsv results/qc_summary.tsv \
  --outdir plots \
  -profile conda
```

---

## 📖 Documentation by Topic

### Setup & Configuration
| Document | Description |
|----------|-------------|
| [Reproducibility Guide](REPRODUCIBILITY_GUIDE.md) | Installation, setup, and environment configuration |
| [Troubleshooting](troubleshooting/TROUBLESHOOTING.md) | Solutions for common issues (conda, memory, paths, etc.) |

### Pipeline Usage
| Document | Description |
|----------|-------------|
| [Pipeline Steps](pipeline_steps.md) | Detailed workflow and process descriptions |
| [Peak Ranking](features/PEAK_RANKING.md) | How peaks are ranked (score vs q-value) |
| [Distance Mode](features/DISTANCE_MODE.md) | Distance-based primer design workflow |

### Advanced Features
| Document | Description |
|----------|-------------|
| [Isoform Optimization](features/ISOFORM_OPTIMIZATION.md) | Optimize primer selection for maximum isoform coverage |
| [Plotting Optimization](features/PLOTTING_OPTIMIZATION.md) | Fast parallel plotting (40× speedup) |
| [Standalone App Guide](STANDALONE_GUIDE.md) | Interactive Shiny app for exploring and visualizing results |
| [HPC RStudio Guide](HPC_RSTUDIO_GUIDE.md) | Running the standalone app on UGent HPC RStudio Server |

### Technical Reference
| Document | Description |
|----------|-------------|
| [Technical Appendix](TECHNICAL_APPENDIX.md) | File formats, specifications, and parameters |
| [Manuscript](MANUSCRIPT.md) | Scientific methodology and validation |

### Examples
| Document | Description |
|----------|-------------|
| [QC Analysis Example](examples/ENSG00000108960_ANALYSIS.md) | Understanding QC failures and exonic filtering |

---

## 🔧 Workflow Modes

### Mode 1: Peak-Based Primer Design (Default)
**Input:** BAM file + GTF + gene list  
**Process:** MACS2 peak calling → peak filtering → primer design  
**Best for:** Standard RNA-seq data with full transcript coverage  

**Key parameters:**
- `--bam` - Aligned RNA-seq reads
- `--qvalue_threshold` - Peak significance threshold (default: 0.05)
- `--peak_rank` - Which peak to select (1=best, 2=second-best, etc.)

### Mode 2: Distance-Based Primer Design
**Input:** Gene list (or transcript FASTA)  
**Process:** MANE transcript fetch → 3' extraction → primer design  
**Best for:** 3' RNA-seq protocols, no BAM file available  

**Key parameters:**
- `--distance_mode` - Enable distance mode
- `--template_length` - Distance from 3' end (default: 300bp)

### Mode 3: Standalone Plotting
**Input:** Pre-computed results (BigWig, peaks, primers)  
**Process:** Generate plots only  
**Best for:** Re-plotting with different genes or parameters  

**Key parameters:**
- `--makeplots` - Enable plotting
- `--bw` - BigWig coverage file
- `--peaks_tsv`, `--primer_targets_bed`, `--qc_tsv` - Result files

---

## 🆘 Common Issues

### Quick Troubleshooting

| Issue | Solution | Reference |
|-------|----------|-----------|
| `Rscript: command not found` | Clean conda cache and retry | [Conda Issues](troubleshooting/TROUBLESHOOTING.md#conda-environment-issues) |
| `numpy.dtype size changed` | Use locked MACS2 environment | [MACS2 NumPy](troubleshooting/TROUBLESHOOTING.md#macs2-numpy-compatibility) |
| Out of memory (exit 137) | Reduce `maxForks` in config | [Memory Issues](troubleshooting/TROUBLESHOOTING.md#memory-issues) |
| No peaks found | Adjust q-value threshold | [Workflow Errors](troubleshooting/TROUBLESHOOTING.md#common-workflow-errors) |
| File path errors | Use absolute paths | [File Paths](troubleshooting/TROUBLESHOOTING.md#file-path-issues) |

**Full troubleshooting guide:** [troubleshooting/TROUBLESHOOTING.md](troubleshooting/TROUBLESHOOTING.md)

---

## 📊 Output Files

### Main Output Files
```
results/
├── primer_targets.bed         # Primer target regions
├── primer_targets.fasta       # Target sequences
├── peaks_qc_summary.tsv       # QC metrics for all genes
├── best_primers.tsv           # Best primer per gene
├── best_primers_optimal.tsv   # Isoform-optimized primers (optional)
├── plots/                     # Visualization plots (if --makeplots)
│   └── plot_ENSG*.pdf
├── macs2_peaks/               # Peak calling results
│   └── *_peaks.narrowPeak
└── primer3_output/            # Primer3 results
    └── *_primer3.txt
```

**Detailed specifications:** [TECHNICAL_APPENDIX.md](TECHNICAL_APPENDIX.md)

---

## 🎓 Learning Path

### New Users
1. Start with [Pipeline Overview](pipeline_steps.md)
2. Follow [Reproducibility Guide](REPRODUCIBILITY_GUIDE.md) for setup
3. Run your first analysis (see Quick Start above)
4. Check [Troubleshooting](troubleshooting/TROUBLESHOOTING.md) if issues arise

### Advanced Users
1. Explore [Peak Ranking](features/PEAK_RANKING.md) for peak selection
2. Learn [Distance Mode](features/DISTANCE_MODE.md) for 3' RNA-seq
3. Optimize with [Isoform Optimization](features/ISOFORM_OPTIMIZATION.md)
4. Speed up plotting with [Plotting Optimization](features/PLOTTING_OPTIMIZATION.md)

### Developers
1. Review [Technical Appendix](TECHNICAL_APPENDIX.md) for specifications
2. Read [Manuscript](MANUSCRIPT.md) for methodology
3. Check [QC Analysis Example](examples/ENSG00000108960_ANALYSIS.md) for QC logic

---

## 📈 Performance Tips

### For Large Datasets
- **Enable plotting optimization:** Already implemented, provides 40× speedup
- **Increase parallelization:** Adjust `maxForks` in `nextflow.config`
- **Use GTF filtering:** Automatically enabled for plotting
- **Pre-filter genes:** Only include genes of interest

### For HPC Clusters
- Use appropriate profile: `-profile slurm` or `-profile pbs`
- Request sufficient resources: See [Reproducibility Guide](REPRODUCIBILITY_GUIDE.md)
- Use mamba for faster conda: Set `conda.useMamba = true` in config

**Details:** [Plotting Optimization](features/PLOTTING_OPTIMIZATION.md)

---

## 🔗 External Resources

### Nextflow
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

### Tools Used
- [MACS2](https://github.com/macs3-project/MACS) - Peak calling
- [Primer3](https://primer3.org/) - Primer design
- [SAMtools](http://www.htslib.org/) - BAM processing
- [Bioconductor](https://bioconductor.org/) - R analysis tools

---

## 📝 Citation

If you use PeakPrime in your research, please cite:

> [Citation to be added - see MANUSCRIPT.md]

---

## 📞 Support

- **Issues:** Check [Troubleshooting](troubleshooting/TROUBLESHOOTING.md) first
- **Questions:** Open an issue on GitHub
- **Bug reports:** Include `.nextflow.log` and command used

---

## 📅 Last Updated

October 21, 2025

---

## 🗂️ Complete File List

```
docs/
├── README.md                    # This file
├── MANUSCRIPT.md                # Scientific paper
├── pipeline_steps.md            # Pipeline overview
├── TECHNICAL_APPENDIX.md        # Technical specifications
├── REPRODUCIBILITY_GUIDE.md     # Setup guide
├── STANDALONE_GUIDE.md          # Interactive Shiny app guide
├── HPC_RSTUDIO_GUIDE.md         # HPC RStudio Server guide
│
├── features/
│   ├── PEAK_RANKING.md          # Peak ranking explained
│   ├── DISTANCE_MODE.md         # Distance-based workflow
│   ├── ISOFORM_OPTIMIZATION.md  # Isoform optimization
│   └── PLOTTING_OPTIMIZATION.md # Plotting performance
│
├── troubleshooting/
│   └── TROUBLESHOOTING.md       # Comprehensive troubleshooting
│
└── examples/
    └── ENSG00000108960_ANALYSIS.md  # QC analysis example
```
