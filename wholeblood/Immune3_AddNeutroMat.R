## add neutrophil subcluster information into immune.combined3@meta.data

library(Seurat)
library(tidyverse)

load('RData/20250218immune.combined3.RData')
load('RData/Neutrophil3.RData')

immune.combined3@meta.data$Neutro_MatureGroups <- as.vector(immune.combined3@meta.data$celltype2)
immune.combined3@meta.data$Neutro_MatureClusters <- as.vector(immune.combined3@meta.data$celltype2)

## if cell barcodes from immune.combiend3 and neutrophil object identical or not
neutrophil_cells <- rownames(immune.combined3@meta.data %>% filter(celltype2 == 'Neutrophil'))
identical(neutrophil_cells, rownames(neutrophil@meta.data))

## add mature groups
immune.combined3@meta.data[neutrophil_cells, 'Neutro_MatureGroups'] <- paste0('Neutrophil_', neutrophil@meta.data$mature_groups)
immune.combined3@meta.data[neutrophil_cells, 'Neutro_MatureClusters'] <- paste0('Neutrophil_', neutrophil@meta.data$mature_clusters)

save(immune.combined3, file = 'RData/Immune.combined3_Neutro.RData')