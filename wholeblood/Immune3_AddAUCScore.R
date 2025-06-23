## add AUC score for IFN, Cytokine

library(AUCell)

load('RData/20250218immune.combined3.RData')

IFN_100_list <- read.table('Genelist/Genelist_IFN100.txt')$V1
Cytokine_list <- read.table('Genelist/Genelist_Immune_Cytokines.txt')$V1

genesets <- list(IFN_100 = IFN_100_list,
                 Cytokine = Cytokine_list)

exprMatrix_RNA <- immune.combined3@assays$RNA@data
RNA_cells_rankings <- AUCell_buildRankings(exprMatrix_RNA, plotStats=TRUE)
RNA_cells_AUC <- AUCell_calcAUC(genesets, RNA_cells_rankings)

if(identical(colnames(immune.combined3), colnames(getAUC(RNA_cells_AUC)))){
   AUC_assays <- t(getAUC(RNA_cells_AUC))
   colnames(AUC_assays) <- paste('AUCell', colnames(AUC_assays), sep = '_')
   immune.combined3@meta.data <- cbind(immune.combined3@meta.data, AUC_assays)}

save(immune.combined3, file = 'RData/Immune.combined3_AUCScore.RData')
