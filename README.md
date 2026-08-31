# Caucasus Kefir Metagenomics

Analysis scripts for **"Shotgun metagenomic analysis reveals regional and strain-level structuring in traditional Caucasus kefir"** (submitted to *Microbial Genomics*, manuscript # MGEN-D-26-00317).

This repository contains the Unix/HPC pipeline scripts and R scripts used for read preprocessing, taxonomic/functional/strain profiling, metagenomic assembly and MAG recovery, and all downstream statistical analyses and visualizations reported in the manuscript.

## Repository structure

Scripts are organized into two top-level folders, each with numbered subfolders corresponding to the same analysis stages:

```
unix_scripts/
├── 01_preprocessing/         # Read QC, trimming, host/contaminant removal
├── 02_assembly/              # De novo metagenomic assembly
├── 03_bin_classification/    # Contig binning, MAG quality assessment, taxonomic classification
├── 04_taxonomic_profiling/   # MetaPhlAn4 species-level profiling
├── 05_functional_profiling/  # HUMAnN3 pathway profiling
├── 06_strain_profiling/      # StrainPhlAn4 strain-level phylogenies
└── 07_downstream_analysis/   # Supporting HPC/command-line steps for downstream analyses

r_scripts/
├── 04_taxonomic_profiling/   # Taxonomic composition and diversity analyses (MetaPhlAn4 output)
├── 05_functional_profiling/  # Functional pathway analyses (HUMAnN3 output)
├── 06_strain_profiling/      # Strain-level distance and clustering analyses (StrainPhlAn4 output)
└── 07_downstream_analysis/   # Statistical analyses and figure generation
```

Numbering follows the order of the Bioinformatics Analysis and Statistical Analysis sections of the manuscript (Methods, sections 6.1–6.4).

## Software and key package versions

| Step | Tool | Version |
|---|---|---|
| Read QC/trimming | Trim Galore | v0.6.10 |
| Host/contaminant read removal | Bowtie2 | v2.5.4 |
| Taxonomic profiling | MetaPhlAn4 (database: mpa_vOct22_CHOCOPhlAnSGB_202403) | v4.1.1 |
| Functional profiling | HUMAnN3 | v3.8 |
| Strain profiling | StrainPhlAn4 | v4.0.6 |
| Metagenomic assembly | SPAdes (metaSPAdes mode) | v3.15.5 |
| Contig filtering | seqkit | v2.8.0 |
| Contig binning | metaWRAP | v1.3.2 |
| MAG quality assessment | CheckM2 | v1.1.0 |
| MAG taxonomic classification | GTDB-Tk (GTDB release 226) | v2.4.0 |
| MAG functional annotation | DRAM | v1.5.0 |
| Antibiotic resistance gene identification | Abricate (CARD database, release 2025-01-14) | v1.0.1 |
| Statistical analysis | R / RStudio | v4.5.1 |
| Diversity/ordination analyses | vegan | v2.7-2 |
| Statistical tests/plots | ggpubr | v0.6.2 |
| Phylogenetic distance/PCoA | ape | v5.8-1 |
| Post hoc testing | rstatix | 1.1.0 |
| Data visualization | ggplot2 | v4.0.0 |

Exact parameters and command-line options for each step are described in the manuscript Methods (Bioinformatics Analysis; Statistical Analysis and Data Visualization).

## Data availability

Raw sequencing data used in this study are available from the European Nucleotide Archive (ENA) at EMBL-EBI under project accession number **PRJEB111907**.

## Citation

If you use these scripts, please cite the associated manuscript (citation to be updated upon publication):

> Ancla J.B., et al. Shotgun metagenomic analysis reveals regional and strain-level structuring in traditional Caucasus kefir. *Microbial Genomics* (submitted).

## Contact

For questions about these scripts, please open an issue on this repository or contact the corresponding author.
