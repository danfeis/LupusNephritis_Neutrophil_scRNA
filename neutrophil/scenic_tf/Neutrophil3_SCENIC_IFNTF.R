## 1. Identification of TFs regulating ISGs with SCENIC TF-target inference.
## 2. Compute ISG score for each TF in each group

library(Seurat)
library(tidyverse)
library(AUCell)
library(RColorBrewer)

load('RData/Neutrophil3.RData')

##---------------------------------------------- 1. SCENIC, Neu_TF and IFN
## SCENIC
TF_target <- read_tsv('Genelist/Genelist3_TFTarget_SCENIC_All.tsv')
TF_target_sig <- TF_target %>% filter(Coef > 0.01) # filter weak regulations, 297 TFs with 106,971 TF-target pairs

## Neu_TF and IFN
Neu_TF <- read.table('Genelist/Genelist3_TFlist_SCENIC_All.txt', header = 1)
Neu_TF_All <- Neu_TF$TF
Neu_TF_ImmG1 <- (Neu_TF %>% filter(TFGroup=='Imm_G1'))$TF # 4
Neu_TF_ImmG2 <- (Neu_TF %>% filter(TFGroup=='Imm_G2'))$TF # 22
Neu_TF_Mat <- (Neu_TF %>% filter(TFGroup=='Mat'))$TF # 7
Neu_TF_Common <- (Neu_TF %>% filter(TFGroup=='Common'))$TF # 12
IFN_100_list <- read.table('Genelist/Genelist_IFN100.txt')$V1


##----------------------------------------------- IFN score of each TF
# selecting TFs regulting ISGs in neutrophils
# 1. TF in Neu_TF
# 2. IFN as TargetGene
# 3. Create new IFN list subsets, one regulated by Imm Neu_TFs (ImmG1, Imm_G2), one regulated by Mat Neu_TFs (Mat, Common).
IFN_target_df_TF <- TF_target_sig %>% filter(TF %in% Neu_TF_All, TargetGene %in% IFN_100_list) # 32 TF
IFN_target_df_Imm <- TF_target_sig %>% filter(TF %in% c(Neu_TF_ImmG1,Neu_TF_ImmG2), TargetGene %in% IFN_100_list) # 81, 13 TFs and 15 targets, with 31 TF-target pairs
IFN_target_df_Mat <- TF_target_sig %>% filter(TF %in% c(Neu_TF_Mat, Neu_TF_Common), TargetGene %in% IFN_100_list) # 3829, 7 TFs and 69 ISGs, with 601 TF-target pairs
# IFN_target_df_Mat <- TF_target_sig %>% filter(TF %in% Neu_TF_Mat, TargetGene %in% IFN_100_list) # 3829, 7 TFs and 60 ISGs, with 182 TF-target pairs
# IFN_target_df_Common <- TF_target_sig %>% filter(TF %in% Neu_TF_Common, TargetGene %in% IFN_100_list) # 11595, 12 TFs and 67 ISGs, with 11,595 TF-target pairs

IFN_Imm_list <- unique(IFN_target_df_Imm$TargetGene)
IFN_Mat_list <- unique(IFN_target_df_Mat$TargetGene)

## expr rank
exprMatrix_RNA <- neutrophil@assays$RNA@data
exprMatrix_RNA_ranked <- AUCell_buildRankings(exprMatrix_RNA, plotStats=FALSE)

geneSets1 <- list(IFN_Imm = IFN_Imm_list)
geneSets2 <- list(IFN_Mat = IFN_Mat_list)
geneSets3 <- list(IFN_100 = IFN_100_list)



## IFN score
IFN_score_df <- data.frame()

for(tf in Neu_TF_All){
   TF_target_df <- TF_target_sig %>% filter(TF==tf)
   TF_target_list <- unique(TF_target_df$TargetGene)
   TF_target_expr_ranked <- exprMatrix_RNA_ranked[TF_target_list,]

   print(tf)

   # immature
   TF_target_AUC_immature_assays <- tryCatch(
      {
         print('immature')
         TF_target_AUC_immature <- AUCell_calcAUC(geneSets1, TF_target_expr_ranked,  aucMaxRank  = nrow(TF_target_expr_ranked)) # less then 20% genes in geneSets are included
         TF_target_AUC_immature_assays <- getAUC(TF_target_AUC_immature)
         rownames(TF_target_AUC_immature_assays) <- paste(get('tf'), rownames(TF_target_AUC_immature_assays), sep = '_')
         
         TF_target_AUC_immature_assays
      },
      error = function(e) {
         print('immature error')
         TF_target_AUC_immature_assays <- data.frame(matrix(0, nrow=1, ncol = ncol(TF_target_expr_ranked)))
         rownames(TF_target_AUC_immature_assays) <- paste(get('tf'), 'IFN_Imm', sep = '_')
         
         TF_target_AUC_immature_assays
      })

   # mature
   TF_target_AUC_mature_assays <- tryCatch(
      {
         print('mature')
         TF_target_AUC_mature <- AUCell_calcAUC(geneSets2, TF_target_expr_ranked,  aucMaxRank  = nrow(TF_target_expr_ranked)) # less then 20% genes in geneSets are included
         TF_target_AUC_mature_assays <- getAUC(TF_target_AUC_mature)
         rownames(TF_target_AUC_mature_assays) <- paste(get('tf'), rownames(TF_target_AUC_mature_assays), sep = '_')
         
         TF_target_AUC_mature_assays
      },
      error = function(e) {
         print('mature error')
         TF_target_AUC_mature_assays <- data.frame(matrix(0, nrow=1, ncol = ncol(TF_target_expr_ranked)))
         rownames(TF_target_AUC_mature_assays) <- paste(get('tf'), 'IFN_Mat', sep = '_')
         
         TF_target_AUC_mature_assays
      })
    
   #  # common
   # TF_target_AUC_common_assays <- tryCatch(
   #    {
   #       print('common')
   #       TF_target_AUC_common <- AUCell_calcAUC(geneSets3, TF_target_expr_ranked,  aucMaxRank  = nrow(TF_target_expr_ranked)) # less then 20% genes in geneSets are included
   #       TF_target_AUC_common_assays <- getAUC(TF_target_AUC_common)
   #       rownames(TF_target_AUC_common_assays) <- paste(get('tf'), rownames(TF_target_AUC_common_assays), sep = '_')
         
   #       TF_target_AUC_common_assays
   #    },
   #    error = function(e) {
   #       print('common error')
   #       TF_target_AUC_common_assays <- data.frame(matrix(0, nrow=1, ncol = ncol(TF_target_expr_ranked)))
   #       rownames(TF_target_AUC_common_assays) <- paste(get('tf'), 'IFN_Common', sep = '_')
         
   #       TF_target_AUC_common_assays
   #    })
   
   # IFN 100
   TF_target_AUC_IFN100_assays <- tryCatch(
      {
         print('IFN 100')
         TF_target_AUC_IFN100 <- AUCell_calcAUC(geneSets3, TF_target_expr_ranked,  aucMaxRank  = nrow(TF_target_expr_ranked)) # less then 20% genes in geneSets are included
         TF_target_AUC_IFN100_assays <- getAUC(TF_target_AUC_IFN100)
         rownames(TF_target_AUC_IFN100_assays) <- paste(get('tf'), rownames(TF_target_AUC_IFN100_assays), sep = '_')
         
         TF_target_AUC_IFN100_assays
      },
      error = function(e) {
         print('IFN 100 error')
         TF_target_AUC_IFN100_assays <- data.frame(matrix(0, nrow=1, ncol = ncol(TF_target_expr_ranked)))
         rownames(TF_target_AUC_IFN100_assays) <- paste(get('tf'), 'IFN_100', sep = '_')
         
         TF_target_AUC_IFN100_assays
      })

   colnames(TF_target_AUC_immature_assays) <- colnames(TF_target_expr_ranked)
   colnames(TF_target_AUC_mature_assays) <- colnames(TF_target_expr_ranked)
   colnames(TF_target_AUC_IFN100_assays) <- colnames(TF_target_expr_ranked)


   IFN_score_df <- rbind(IFN_score_df, TF_target_AUC_immature_assays, TF_target_AUC_mature_assays, TF_target_AUC_IFN100_assays)
   rm(TF_target_AUC_immature_assays, TF_target_AUC_mature_assays, TF_target_AUC_IFN100_assays)
}

IFN_score_df <- data.frame(t(IFN_score_df))
neutrophil@meta.data <- cbind(neutrophil@meta.data, IFN_score_df)
TF_IFN_list <- colnames(IFN_score_df)[colSums(IFN_score_df)>0] # 200-->50
TF_IFN_list2 <- unique(gsub('_IFN_100|_IFN_Imm|_IFN_Mat', '', TF_IFN_list)) # 17 TFs significantly reulating ISGs, 1 ImmG2, 5 Mat and 11 Common

save(neutrophil, file = 'RData/Neutrophil3_SCENIC_IFNScore.RData')
write.table(TF_IFN_list2, file = 'Genelist/Genelist3_Neutrophil_TFlist_SCENIC_IFN.txt',quote = FALSE, col.names = FALSE, row.names = FALSE)