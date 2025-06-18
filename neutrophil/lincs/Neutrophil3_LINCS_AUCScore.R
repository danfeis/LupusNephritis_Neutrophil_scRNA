## Compute AUCScore of leading edge of drugs

library(Seurat)
library(AUCell)
library(viridis)
library(ggplot2)
library(tidyverse)
library(cowplot)

load('RData/Neutrophil3.RData')

## LINCS: Neutrophil scDEGs
LINCS_Neu_up_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_UP_sigAnno.txt', header = 1) %>% arrange(NES) %>% filter(NES < 0) # 360
LINCS_Neu_down_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_DOWN_sigAnno.txt', header = 1) %>% arrange(desc(NES)) %>% filter(NES > 0) # 123

## to upper
LINCS_Neu_up_sig$pert_iname <- toupper(LINCS_Neu_up_sig$pert_iname)
LINCS_Neu_down_sig$pert_iname <- toupper(LINCS_Neu_down_sig$pert_iname)

## leading genes of each drug, up and down
leading_geneset_up <- lapply(LINCS_Neu_up_sig$leadingEdge, function(genes){
    genes_split <- strsplit(genes, ', ')[[1]]
})
names(leading_geneset_up) <- LINCS_Neu_up_sig$pert_iname
leading_geneset_down <- lapply(LINCS_Neu_down_sig$leadingEdge, function(genes){
    genes_split <- strsplit(genes, ', ')[[1]]
})
names(leading_geneset_down) <- LINCS_Neu_down_sig$pert_iname

## leading genes of all drugs targeting up- and down-
LINCS_genests <- list(LINCS_up = unique((LINCS_Neu_up_sig %>% separate_rows(leadingEdge, sep = ",\\s*"))$leadingEdge), # all leading edges in drugs targeting updegulated genes
                      LINCS_down = unique((LINCS_Neu_down_sig %>% separate_rows(leadingEdge, sep = ",\\s*"))$leadingEdge)) # all leading edges in drug targeting downregulated genes

## combine
LINCS_genests <- c(LINCS_genests, leading_geneset_up, leading_geneset_down) # 485

## AUCell score
## up
exprMatrix_RNA <- neutrophil@assays$RNA@data
RNA_cells_rankings <- AUCell_buildRankings(exprMatrix_RNA, plotStats=FALSE)
RNA_cells_AUC_up <- AUCell_calcAUC(LINCS_genests, RNA_cells_rankings)
if(identical(colnames(neutrophil), colnames(getAUC(RNA_cells_AUC_up)))){
   AUC_assays_up <- t(getAUC(RNA_cells_AUC_up))
   colnames(AUC_assays_up) <- paste('LINCS', colnames(AUC_assays_up), sep = '_')
   neutrophil@meta.data <- cbind(neutrophil@meta.data, AUC_assays_up)}

save(neutrophil, file = 'RData/Neutrophil3_LINCS_AUCScore.RData')

