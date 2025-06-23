## try to start a new round of clustering of immune
## to check whether CD4 and CD8 T cell clustering correct or not

library(ggplot2)
library(tidyverse)
library(Seurat)

load('RData/20230920immune.combined2.RData')

## clustering
for(d in c(10,15,20,25,30)){
    immune.combined2 <- RunUMAP(immune.combined2, reduction = "harmony", dims = 1:d, return.model=TRUE)
    for(k in c(5,10,15,20)){
            for(r in c(0.3,0.4,0.5,0.6,0.7)){
                immune.combined2 <- FindNeighbors(immune.combined2, reduction = "harmony", dims = 1:d, k.param = k) #20
                immune.combined2 <- FindClusters(immune.combined2, resolution = r, graph.name = 'integrated_snn')
                p <- DimPlot(immune.combined2, group.by = 'seurat_clusters', label = TRUE)
                f <- paste0('UMAP_immune/Immune_d',d,'_k',k,'_reso0',10*r,'.png')
                ggsave(f, p, width = 10, height = 10)}
    }
}



## do not change UMAP, just change cluster with knn and louvain
for(k in c(5,10,15,20)){
    for(r in c(0.3,0.4,0.5,0.6,0.7,0.8)){
        immune.combined2 <- FindNeighbors(immune.combined2, reduction = "harmony", dims = 1:14, k.param = k) #20
        immune.combined2 <- FindClusters(immune.combined2, resolution = r, graph.name = 'integrated_snn')
        p <- DimPlot(immune.combined2, group.by = 'seurat_clusters', label = TRUE)
        f <- paste0('UMAP_immune/Immune2_k',k,'_reso0',10*r,'.png')
        ggsave(f, p, width = 10, height = 10)
    }
}



## parameters ensured
dim = 14
k = 5
reso = 0.7
immune.combined3 <- immune.combined2 %>% 
  #RunUMAP(reduction = "harmony", dims = 1:14) %>% 
  FindNeighbors(reduction = "harmony", k.param = k , dims = 1:dim) %>% 
  FindClusters(resolution = reso,  graph.name = 'integrated_snn')

## labeling
markers <- c("Cd3e", "Cd2", "Cd4", "Cd8a", "Cd79a", "Ms4a1", "Klrb1c", "Gzma",  "Mcpt8", "Ccl3", "Ly6g", "Lcn2", "S100a8","Adgre1", "Csf1r", "Nrp1", "Siglech", "Pf4", "Pdpd", "Hbb-bt")
p <- VlnPlot(immune.combined3, markers, pt.size=0, stack = T)
ggsave('UMAP_immune/VlnPlot_Immune2_k5_reso07.pdf', p, height = 12, width = 18)


immune.combined3 <- RenameIdents(immune.combined3, '0' = 'DNT', '1' = 'DNT', '2' = 'DNT', '8' = 'DNT', '14' = 'DNT', '16' = 'DNT', '20' = 'DNT', 
                                                   '11' = 'CD4+T',
                                                   '5' = 'CD8+T', '18' = 'CD8+T',
                                                   '6' = 'B', '25' = 'B', '26' = 'B', '27' = 'B',
                                                   '17' = 'NK',
                                                   '21' = 'Basophil',
                                                   '3' = 'Neutrophil', '4' = 'Neutrophil', '7' = 'Neutrophil', '10' = 'Neutrophil', '19' = 'Neutrophil', '22' = 'Neutrophil',
                                                   '13' = 'Monocyte','15' = 'Monocyte',
                                                   '23' = 'pDC',
                                                   '9' = 'Megaka',
                                                   '12' = 'Erythroid', '24' = 'Erythroid')

################################# Annotate and Factorize
# Annotation: celltype, celltype2
immune.combined3@meta.data$celltype2 <- Idents(immune.combined3) # cell type
immune.combined3@meta.data$celltype2.type <- paste(immune.combined3@meta.data$celltype2, immune.combined3@meta.data$type, sep = "_") # cell type with disease stage
immune.combined3@meta.data$celltype <- immune.combined3@meta.data %>%
  group_by(celltype2) %>%
  mutate(celltype = paste0(celltype2, match(seurat_clusters, unique(seurat_clusters)))) %>%
  ungroup() %>%
  pull(celltype)  # # cell type with subclusters, extracts the column as a vector to avoid tibble issues

## Factor my_type, samples, celltype.type
my_type <- c('mild','moderate','severe')
my_samples <- c('MNL3181BL','MNL4186BL','MML2182BL','MML41810BL','MSL7182BL','MSL4187BL')
my_celltype2type <- paste(rep(levels(immune.combined3@meta.data$celltype2),each = 3), c('mild','moderate','severe'), sep = '_')

immune.combined3@meta.data$orig.ident <- factor(immune.combined3@meta.data$orig.ident, levels = my_samples) ## factor samples
immune.combined3@meta.data$type <- factor(immune.combined3@meta.data$type, levels = my_type) ## factor samples
immune.combined3@meta.data$celltype2.type <- factor(immune.combined3@meta.data$celltype2.type, levels = my_celltype2type) ## factor celltype with disease stages
  
save(immune.combined3, file = 'RData/20250218immune.combined3.RData')




# ## parameters ensured
# dim = 14
# k = 20
# reso = 0.8
# immune.combined3 <- immune.combined2 %>% 
#   #RunUMAP(reduction = "harmony", dims = 1:14) %>% 
#   FindNeighbors(reduction = "harmony", k.param = k , dims = 1:dim) %>% 
#   FindClusters(resolution = reso,  graph.name = 'integrated_snn')
# DimPlot(immune.combined3, label = T)

# ## labeling
# markers <- c("Cd3e", "Cd2", "Cd4", "Cd8a", "Cd79a", "Ms4a1", "Klrb1c", "Gzma",  "Mcpt8", "Ccl3", "Ly6g", "Lcn2", "S100a8","Adgre1", "Csf1r", "Nrp1", "Siglech", "Pf4", "Pdpd", "Hbb-bt")
# p <- VlnPlot(immune.combined3, markers, pt.size=0, stack = T)
# ggsave('UMAP_immune/VlnPlot_Immune2_k20_reso08.pdf', p, height = 12, width = 18)

# immune.combined3 <- RenameIdents(immune.combined3, '0' = 'DNT', '1' = 'DNT', '2' = 'DNT', '5' = 'DNT', '15' = 'DNT', '16' = 'DNT', '17' = 'DNT', '22' = 'DNT', 
#                                                    '12' = 'CD4+T',
#                                                    '7' = 'CD8+T', '20' = 'CD8+T',
#                                                    '6' = 'B', '27' = 'B',
#                                                    '19' = 'NK',
#                                                    '23' = 'Basophil',
#                                                    '3' = 'Neutrophil', '4' = 'Neutrophil', '8' = 'Neutrophil', '9' = 'Neutrophil', '11' = 'Neutrophil', '21' = 'Neutrophil', '24' = 'Neutrophil',
#                                                    '14' = 'Monocyte','18' = 'Monocyte',
#                                                    '25' = 'pDC',
#                                                    '10' = 'Megaka',
#                                                    '13' = 'Erythroid', '26' = 'Erythroid')
# immune.combined3@meta.data$celltype2 <- Idents(immune.combined3) # cell type
# DimPlot(immune.combined3, label = T)
# ## check neutrophils identical or not
# neutrophil2 <- subset(immune.combined2, ident = c('Neutrophil')) # 14713
# neutrophil3 <- subset(immune.combined3, ident = c('Neutrophil')) # 14706
# diff_cells1 <- setdiff(colnames(neutrophil2),colnames(neutrophil3)) # cells not identified as 'neutrophils' in new version of clustering
# diff_cells2 <- setdiff(colnames(neutrophil3),colnames(neutrophil2)) # cells were not identified as 'neutrophils' in previous version of clustering
# immune.combined3@meta.data[diff_cells1,] # 11, most were re-clustered and labeled into DNT
# immune.combined2@meta.data[diff_cells2,] # 4, these were clustered into DNT,CD4



# ## parameters ensured
# dim = 14
# k = 10
# reso = 0.8
# immune.combined3 <- immune.combined2 %>% 
#   #RunUMAP(reduction = "harmony", dims = 1:14) %>% 
#   FindNeighbors(reduction = "harmony", k.param = k , dims = 1:dim) %>% 
#   FindClusters(resolution = reso,  graph.name = 'integrated_snn')
# DimPlot(immune.combined3, label = T)

# ## labeling
# markers <- c("Cd3e", "Cd2", "Cd4", "Cd8a", "Cd79a", "Ms4a1", "Klrb1c", "Gzma",  "Mcpt8", "Ccl3", "Ly6g", "Lcn2", "S100a8","Adgre1", "Csf1r", "Nrp1", "Siglech", "Pf4", "Pdpd", "Hbb-bt")
# p <- VlnPlot(immune.combined3, markers, pt.size=0, stack = T)
# ggsave('UMAP_immune/VlnPlot_Immune2_k10_reso08.pdf', p, height = 12, width = 18)

# immune.combined3 <- RenameIdents(immune.combined3, '0' = 'DNT', '1' = 'DNT', '2' = 'DNT', '3' = 'DNT', '5' = 'DNT', '13' = 'DNT', '15' = 'DNT', '16' = 'DNT', '17' = 'DNT', '23' = 'DNT', 
#                                                    '11' = 'CD4+T',
#                                                    '8' = 'CD8+T', '19' = 'CD8+T',
#                                                    '7' = 'B', '30' = 'B', '31' = 'B',
#                                                    '18' = 'NK',
#                                                    '24' = 'Basophil',
#                                                    '4' = 'Neutrophil', '6' = 'Neutrophil', '9' = 'Neutrophil', '12' = 'Neutrophil', '20' = 'Neutrophil', '22' = 'Neutrophil', '26' = 'Neutrophil', '27' = 'Neutrophil', 
#                                                    '14' = 'Monocyte','17' = 'Monocyte',
#                                                    '28' = 'pDC',
#                                                    '10' = 'Megaka',
#                                                    '21' = 'Erythroid', '25' = 'Erythroid', '29' = 'Erythroid')
# DimPlot(immune.combined3, label = T)
# ## check neutrophils identical or not
# neutrophil2 <- subset(immune.combined2, ident = c('Neutrophil')) # 14713
# neutrophil3 <- subset(immune.combined3, ident = c('Neutrophil')) # 14715
# diff_cells1 <- setdiff(colnames(neutrophil2),colnames(neutrophil3)) # cells not identified as 'neutrophils' in new version of clustering
# diff_cells2 <- setdiff(colnames(neutrophil3),colnames(neutrophil2)) # cells were not identified as 'neutrophils' in previous version of clustering
# immune.combined3@meta.data[diff_cells1,] # 18, most were re-clustered and labeled into neutrophils
# immune.combined2@meta.data[diff_cells2,] # 20, these were clustered into DNT,CD4,CD8




# ## check cluster 24 markers
# markers_cluster24 <- FindMarkers(immune.combined3, ident.1 = 24, only.pos = TRUE)

# ## check neutrophils identical or not
# neutrophil2 <- subset(immune.combined2, ident = c('Neutrophil')) # 14713
# neutrophil3 <- subset(immune.combined3, ident = c('Neutrophil')) # 14618
# diff_cells1 <- setdiff(colnames(neutrophil2),colnames(neutrophil3)) # cells not identified as 'neutrophils' in new version of clustering
# diff_cells2 <- setdiff(colnames(neutrophil3),colnames(neutrophil2)) # cells were not identified as 'neutrophils' in previous version of clustering

# immune.combined3@meta.data[diff_cells1,] # 110, most were re-clustered and labeled into DNT
# immune.combined2@meta.data[diff_cells2,] # 15, these were clustered into CD8+, CD4+ and DNT

# # Create a new column in metadata for coloring
# immune.combined3$highlight <- "Other"  # Default to grey
# # Assign specific colors to identified cells
# immune.combined3$highlight[colnames(immune.combined3) %in% diff_cells1] <- "Neu_Old"
# immune.combined3$highlight[colnames(immune.combined3) %in% diff_cells2] <- "Neu_New"
# # Convert to factor for correct ordering in the plot
# immune.combined3$highlight <- factor(immune.combined3$highlight, levels = c("Other", "Neu_Old", "Neu_New"))
# DimPlot(immune.combined3, group.by = "highlight", cols = c("Other" = "grey", "Neu_Old" = "red", "Neu_New" = "blue")) +
#   ggtitle("Highlighted Cells in Red and Blue")
