# LupusNephritis_Neutrophil_scRNA

This repository contains code and analysis from the project **"Single-cell profiling of neutrophil heterogeneity in transcriptomics and immunometabolism underlying murine lupus nephritis progression reveals potential treatments."** The work was conducted as part of Lupus Nephritis multi-omics studies from Xin Sheng's lab, Liangzhu laboratory, Zhejiang University.

## 🔬 Project Overview

Systemic lupus erythematosus (SLE) and its severe renal manifestation, lupus nephritis (LN), involve complex immune dysregulation. This project leverages whole blood single-cell RNA sequencing (scRNA-seq) to uncover transcriptional heterogeneity, ISGs, metabolism, and druggable features of **neutrophils** during LN progression in a murine model.

Key analytical components include:
- Cell clustering and annotation
- Pseudotime and trajectory inference (Monocle3)
- RNA velocity analysis
- Transcription factor inference (SCENIC)
- CUT&Tag integration and correlation
- Metabolism post-GWAS analysis and functional enrichment (MAGMA, KEGG)
- Cell-cell communication (CellChat)
- LINCS-based druggability analysis

## 📁 Repository Structure

```bash
├── README.md
├── environment.yml             # Conda environment setup
├── .gitignore
│
├── neutrophil/                 # All neutrophil-specific analyses
│   ├── LINCS_AUC/              # Drug response via LINCS & AUC
│   ├── trajectory_analysis/    # Monocle3 pseudotime and gene trends
│   ├── rnaseq_velocity/        # RNA velocity analysis
│   ├── cell_communication/     # CellChat signaling analysis
│   ├── scenic_tf/              # SCENIC transcription factor activity
│   ├── cut_tag/                # CUT&Tag integration & correlation
│   ├── functional_enrichment/  # MAGMA, KEGG metabolism
│   ├── deg_hvg/                # HVG and DE gene identification
│   └── clustering/             # UMAP projection and clustering
│
├── wholeblood/                 # Whole blood single-cell analysis
│   ├── clustering/             # UMAP projection, clustering, immune landscape
│   ├── sample_comparison/      # Group-based comparison and proportions
│   ├── gene_modules/           # Global ISG & inflammatory modules
│   └── quality_control/        # QC, doublet removal, filtering steps
│
└── data/                       # Placeholder for processed data or metadata
```

## 🔧 Environments

This project uses different Conda environments for R and Python components.

- R 4.1: Used for Seurat, SCENIC, Monocle3 and CellChat analysis → `environment_r41.yml`
- R 4.2: Used for CUT&Tag ananotation → `environment_r42.yml`

Create and activate environments as needed:

```bash
conda env create -f environment_r41.yml
conda activate r41_env
```
