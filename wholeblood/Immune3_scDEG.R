## find scDEG from immune.combined3 without cell types
## keep all protein coding genes (PCGs)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rtracklayer)

load('RData/20250218immune.combined3.RData')

####################################################################################### scDEGs were identified from mild and severe (Severe VS. Mild) disease stages with Seurat
cells_mild <- row.names(immune.combined3@meta.data %>% filter(type=='mild'))
cells_severe <- row.names(immune.combined3@meta.data %>% filter(type=='severe'))

scDEG <- FindMarkers(immune.combined3, ident.1 = cells_severe, ident.2 = cells_mild) # 9262, severe versus mild, same as scDEGs identified from immune.combined2
scDEG_upregulate <- scDEG %>% filter(p_val_adj < 0.05, avg_log2FC > 1) # 1182
scDEG_downregulate <- scDEG %>% filter(p_val_adj < 0.05, avg_log2FC <(-1)) # 315

write.table(scDEG, file = 'Genelist/Genelist3_Immune_scDEGs.txt', quote = FALSE)

####################################################################################### scDEGs in PCGs
## extract pcg from mm10 gtf file
gtf_file <- '/public/home/Shenglab/DataShare/Reference/mouse/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf'
gtf <- import(gtf_file)
pcg <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"] # protein coding genes: 21700
pcg_symbol <- pcg %>% data.frame() %>% select(gene_id, gene_name)

## scDEGs in PCGs
scDGE_PCG_ind <- rownames(scDEG) %in% pcg_symbol$gene_name
scDEG_PCG <- scDEG[scDGE_PCG_ind,] # subset rows for pcgs, 8734
scDEG_PCG_upregulate <- scDEG_PCG %>% filter(p_val_adj < 0.05, avg_log2FC > 1) # 1047
scDEG_PCG_downregulate <- scDEG_PCG %>% filter(p_val_adj < 0.05, avg_log2FC <(-1)) # 293

write.table(scDEG_PCG, file = 'Genelist/Genelist3_Immune_scDEGs_PCG.txt', quote = FALSE)
write.table(rownames(scDEG_PCG_upregulate), file = 'Genelist/Genelist3_Immune_scDEGs_PCG_UpregulateSymbol.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(scDEG_PCG_downregulate), file = 'Genelist/Genelist3_Immune_scDEGs_PCG_DownregulateSymbol.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
