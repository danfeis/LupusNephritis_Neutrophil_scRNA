# LupusNephritis_Neutrophil_scRNA

This repository contains code and analysis from the project **"Single-cell profiling of neutrophil heterogeneity in transcriptomics and immunometabolism underlying murine lupus nephritis progression reveals potential treatments."** The work was conducted as part of **Lupus Nephritis** studies from Xin Sheng's lab from Liangzhu laboratory, Zhejiang University.

## 🔬 Project Overview

Systemic lupus erythematosus (SLE) and its severe renal manifestation, lupus nephritis (LN), involve complex immune dysregulation. This project leverages single-cell RNA sequencing (scRNA-seq) to uncover transcriptional heterogeneity, metabolic states, and druggable features of neutrophils during LN progression in a murine model.

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
├── LINCS_AUC/                  # Drug response analysis via LINCS & AUC
├── trajectory_analysis/        # Monocle3 pseudotime and gene trends
├── rnaseq_velocity/            # RNA velocity analysis
├── cell_communication/         # CellChat signaling analysis
├── scenic_tf/                  # SCENIC transcription factor activity
├── cut_tag/                    # CUT&Tag integration & correlation
├── functional_enrichment/      # MAGMA, KEGG metabolism
├── dimensionality_reduction/   # UMAP and dimensionality projection
└── data/                       # Placeholder for processed data or metadata
