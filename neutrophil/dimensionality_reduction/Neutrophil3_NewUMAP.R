## subset neutrophils from immune.combined3 and re-clustered

library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)

load('RData/20250218immune.combined3.RData')

## subset neutrophils
neutrophil <- subset(immune.combined3, ident = 'Neutrophil')

## re clustering
neutrophil <- NormalizeData(neutrophil)
neutrophil <- FindVariableFeatures(neutrophil, selection.method = "vst", nfeatures = 2000)
neutrophil <- ScaleData(neutrophil, verbose = FALSE)
neutrophil <- RunPCA(neutrophil, npcs = 30, verbose = FALSE)
neutrophil <- RunHarmony(neutrophil, group.by.vars = "orig.ident", plot_convergence = TRUE)
ElbowPlot(neutrophil, reduction = 'harmony')
neutrophil <- RunUMAP(neutrophil, reduction = "harmony", dims = 1:15, return.model=TRUE)
DimPlot(neutrophil)

for(k in c(5,10,15,20)){
        for(r in c(0.3,0.4,0.5,0.6,0.7)){
        neutrophil <- FindNeighbors(neutrophil, reduction = "harmony", dims = 1:15, k.param = k) #20
        neutrophil <- FindClusters(neutrophil, resolution = r)
        p <- DimPlot(neutrophil, group.by = 'seurat_clusters', label = TRUE)
        f <- paste0('UMAP_Neutrophil/Neutrophil_k',k,'_reso0',10*r,'.png')
        ggsave(f, p, width = 10, height = 10)}
}

neutrophil <- FindNeighbors(neutrophil, reduction = "harmony", dims = 1:15, k.param = 15) # dim = 15, k = 15
neutrophil <- FindClusters(neutrophil, resolution = 0.7) # resolution = 0.7
DimPlot(neutrophil, label =T)


## neutrophil subcluster relabeling
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='0'] <- 'C7'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='1'] <- 'C5'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='2'] <- 'C8'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='3'] <- 'C9'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='4'] <- 'C6'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='5'] <- 'C2'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='6'] <- 'C4'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='7'] <- 'C0'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='8'] <- 'C3'
neutrophil@meta.data$cluster[neutrophil@meta.data$seurat_clusters=='9'] <- 'C1'

neutrophil@meta.data$cluster <- factor(neutrophil@meta.data$cluster, levels = paste0('C',0:9))

## neutrophil maturation labeling
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C0'] <- 'Imm_G1'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C1'] <- 'Imm_G1'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C2'] <- 'Imm_G2'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C3'] <- 'Imm_G2'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C4'] <- 'Mat_G3'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C5'] <- 'Mat_G3'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C6'] <- 'Mat_G4'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C7'] <- 'Mat_G4'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C8'] <- 'Mat_G4'
neutrophil@meta.data$mature_groups[neutrophil@meta.data$cluster=='C9'] <- 'Mat_G4'

neutrophil@meta.data$mature_groups <- factor(neutrophil@meta.data$mature_groups, levels = c('Imm_G1','Imm_G2','Mat_G3','Mat_G4'))

## neutrophil maturation cluster labeling
neutrophil@meta.data$mature_clusters <- paste(neutrophil@meta.data$mature_groups, neutrophil@meta.data$cluster, sep = '_')
neutrophil@meta.data$mature_clusters <- factor(neutrophil@meta.data$mature_clusters, levels = c('Imm_G1_C0','Imm_G1_C1','Imm_G2_C2','Imm_G2_C3','Mat_G3_C4','Mat_G3_C5','Mat_G4_C6','Mat_G4_C7','Mat_G4_C8','Mat_G4_C9'))

## set Idents
Idents(neutrophil) <- 'mature_clusters'

## save object
save(neutrophil, file = 'RData/Neutrophil3.RData')
