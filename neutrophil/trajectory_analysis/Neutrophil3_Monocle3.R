## Monocle 3 analysis for Neutrophil version 3

# conda activate R_v3

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(viridis)

load('RData/Neutrophil3.RData')

neutrophil.cds <- as.cell_data_set(neutrophil)
neutrophil.cds <- cluster_cells(cds = neutrophil.cds, reduction_method = "UMAP", k = 5, resolution = 0.7, cores = 10)
neutrophil.cds <- learn_graph(neutrophil.cds, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 8)) # , verbose = TRUE, learn_graph_control = list(minimal_branch_len = 10)
neutrophil.cds <- order_cells(neutrophil.cds, reduction_method = "UMAP")


plot_cells(neutrophil.cds, color_cells_by = "cluster", graph_label_size = 3)


save(neutrophil.cds, file = 'RData/Neutrophil3_Monocle3.RData')