## AUCScore

library(Seurat)
library(AUCell)
library(viridis)
library(ggplot2)
library(tidyverse)
library(cowplot)

load('RData/20250218immune.combined3.RData')

## LINCS immune scDEGs
LINCS_immune_up_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_UP_sigAnno.txt', header=1) # 1794
LINCS_immune_down_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_DOWN_sigAnno.txt', header=1) # 2425

LINCS_immune_up_sig$pert_iname <- toupper(LINCS_immune_up_sig$pert_iname) # to upper
LINCS_immune_down_sig$pert_iname <- toupper(LINCS_immune_down_sig$pert_iname)

## leading genes of each drug, up and down
leading_geneset_up <- lapply(LINCS_immune_up_sig$leadingEdge, function(genes){
    genes_split <- strsplit(genes, ', ')[[1]]
})
names(leading_geneset_up) <- LINCS_immune_up_sig$pert_iname
leading_geneset_down <- lapply(LINCS_immune_down_sig$leadingEdge, function(genes){
    genes_split <- strsplit(genes, ', ')[[1]]
})
names(leading_geneset_down) <- LINCS_immune_down_sig$pert_iname

## leading genes of all drugs targeting up- and down-
LINCS_genests <- list(LINCS_up = unique((LINCS_immune_up_sig %>% separate_rows(leadingEdge, sep = ",\\s*"))$leadingEdge), # all leading edges in drugs targeting updegulated genes
                      LINCS_down = unique((LINCS_immune_down_sig %>% separate_rows(leadingEdge, sep = ",\\s*"))$leadingEdge)) # all leading edges in drug targeting downregulated genes

## combine
LINCS_genests <- c(LINCS_genests, leading_geneset_up, leading_geneset_down) # 4221


## AUCell score
exprMatrix_RNA <- immune.combined3@assays$RNA@data
RNA_cells_rankings <- AUCell_buildRankings(exprMatrix_RNA, plotStats=FALSE)
RNA_cells_AUC_up <- AUCell_calcAUC(LINCS_genests, RNA_cells_rankings)
if(identical(colnames(immune.combined3), colnames(getAUC(RNA_cells_AUC_up)))){
   AUC_assays_up <- t(getAUC(RNA_cells_AUC_up))
   colnames(AUC_assays_up) <- paste('LINCS', colnames(AUC_assays_up), sep = '_')
   immune.combined3@meta.data <- cbind(immune.combined3@meta.data, AUC_assays_up)}


save(immune.combined3, file = 'RData/Immune.combined3_LINCS_AUCScore.RData')