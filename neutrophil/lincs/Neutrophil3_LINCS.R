## Drugability analysis with LINCS database
## we use up- and down-regulated scDEGs and bulkDEGs for LINCS analysis
## we do not integrate up- and down- genes as we wang to know how these gene expression reverted by these drugs
## Neutrophil scDEGs cutoff (padj 0.05, logFC 0.5)
## Immune scDEGs cutoff (padj 0.05, logFC 1)
## bulkDEGs cutoff (FDR 0.05, logFC 1)

library(fgsea)
library(cmapR)
library(doParallel)
library(foreach)
library(tidyverse)
library(dplyr)
library(stringr)

## LINCS
PRL_gmt <- parse_gmt('/public/home/zhudf/LINCS/LINCS_Out2/LINCS_PRL_top250.gmt')

############################################################# 1. LINCS on immune scDEGs, up- and down-regulated separately
## scDEG, immune, PCG
scDEG_up_df <- read.table('Genelist/Genelist3_Immune_scDEGs_PCG.txt') %>% filter(p_val_adj < 0.05, avg_log2FC > 1) # 1047
scDEG_down_df <- read.table('Genelist/Genelist3_Immune_scDEGs_PCG.txt') %>% filter(p_val_adj < 0.05, avg_log2FC < (-1)) # 293

## logFC as stat in fgsea and in decreasing order
scDEG_up_df <- scDEG_up_df %>% arrange(desc(avg_log2FC))
scDEG_down_df <- scDEG_down_df %>% arrange(desc(avg_log2FC))

## scDEG genelists
scDEG_up <- scDEG_up_df$avg_log2FC
names(scDEG_up) <- rownames(scDEG_up_df)
scDEG_down <- scDEG_down_df$avg_log2FC
names(scDEG_down) <- rownames(scDEG_down_df)

## gsea
cl_cores <- 60
cl <- makeCluster(cl_cores)
registerDoParallel(cl)
scDEG_fgseaRes_up <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'scDEG_up') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = scDEG_up,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}
scDEG_fgseaRes_down <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'scDEG_down') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = scDEG_down,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}
stopCluster(cl) 

## list as vector
scDEG_fgseaRes_up$leadingEdge <- sapply(scDEG_fgseaRes_up$leadingEdge, paste, collapse = ", ") # 52372
scDEG_fgseaRes_down$leadingEdge <- sapply(scDEG_fgseaRes_down$leadingEdge, paste, collapse = ", ") # 24089

write.table(scDEG_fgseaRes_up, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_UP_All.txt', sep = '\t', row.names = FALSE)
write.table(scDEG_fgseaRes_down, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_DOWN_All.txt', sep = '\t', row.names = FALSE)

## significant
scDEG_fgseaRes_up_sig <- scDEG_fgseaRes_up %>% filter(padj < 0.05) # 1794
scDEG_fgseaRes_down_sig <- scDEG_fgseaRes_down %>% filter(padj < 0.05) # 2425

## pathway annotation
pert_anno_phase1 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
pert_anno_phase2 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz'), sep = '\t')

scDEG_fgseaRes_up_sig_anno <- lapply(scDEG_fgseaRes_up_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }
    
    return(pert_anno)
})
scDEG_fgseaRes_down_sig_anno <- lapply(scDEG_fgseaRes_down_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }
    return(pert_anno)
})

scDEG_fgseaRes_up_sig_anno_df <- Reduce(rbind, scDEG_fgseaRes_up_sig_anno) # 1794
scDEG_fgseaRes_down_sig_anno_df <- Reduce(rbind, scDEG_fgseaRes_down_sig_anno) # 2425

## combine
identical(scDEG_fgseaRes_up_sig$pathway, scDEG_fgseaRes_up_sig_anno_df$pert_id)
identical(scDEG_fgseaRes_down_sig$pathway, scDEG_fgseaRes_down_sig_anno_df$pert_id)
scDEG_fgseaRes_up_sig_anno_df <- cbind(data.frame(scDEG_fgseaRes_up_sig), scDEG_fgseaRes_up_sig_anno_df) %>% arrange(NES)
scDEG_fgseaRes_down_sig_anno_df <- cbind(data.frame(scDEG_fgseaRes_down_sig), scDEG_fgseaRes_down_sig_anno_df) %>% arrange(desc(NES))

write.table(scDEG_fgseaRes_up_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_UP_sigAnno.txt', sep = '\t', row.names = FALSE)
write.table(scDEG_fgseaRes_down_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_DOWN_sigAnno.txt', sep = '\t', row.names = FALSE)


############################################################# 2. LINCS on immune scDEGs, integrating up- and down-regulated in one list
# ## scDEG
# scDEG_all <- read.table('Genelist/Genelist3_Immune_scDEGs_PCG.txt')
# scDEG_all <- scDEG_all %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>% arrange(desc(avg_log2FC)) # 1340
# scDEG_genelist <- scDEG_all$avg_log2FC
# names(scDEG_genelist) <- rownames(scDEG_all)

# ## gsea
# cl_cores <- 60
# cl <- makeCluster(cl_cores)
# registerDoParallel(cl)
# scDEG_fgseaRes_all <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'scDEG_genelist') %dopar% {
#     prl_gsea <- fgsea(pathways = prl, 
#                       stats    = scDEG_genelist,
#                       minSize  = 10,
#                       maxSize  = 500)
#     prl_gsea$pathway <- prl$head
#     return(prl_gsea)
# }
# stopCluster(cl) 

# ## list as vector
# scDEG_fgseaRes_all$leadingEdge <- sapply(scDEG_fgseaRes_all$leadingEdge, paste, collapse = ", ") # 52381
# write.table(scDEG_fgseaRes_all, 'Genelist/Genelist3_LINCS_GSEA_scDEG_All.txt', sep = '\t', row.names = FALSE)
# ## significant
# scDEG_fgseaRes_all_sig <- scDEG_fgseaRes_all %>% filter(padj < 0.05) %>% arrange(desc(NES)) # 3376
# ## pathway annotation
# pert_anno_phase1 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
# pert_anno_phase2 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz'), sep = '\t')

# scDEG_fgseaRes_all_sig_anno <- lapply(scDEG_fgseaRes_all_sig$pathway, function(pert){
#     if(pert %in% pert_anno_phase1$pert_id){
#         pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
#     } else {
#         pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
#     }
    
#     return(pert_anno)
# })

# scDEG_fgseaRes_all_sig_anno_df <- Reduce(rbind, scDEG_fgseaRes_all_sig_anno) # 3374
# scDEG_fgseaRes_all_sig2 <- scDEG_fgseaRes_all_sig %>% filter(pathway %in% scDEG_fgseaRes_all_sig_anno_df$pert_id) # 3374, remove pathways not in pert anno
# identical(scDEG_fgseaRes_all_sig2$pathway, scDEG_fgseaRes_all_sig_anno_df$pert_id)
# scDEG_fgseaRes_all_sig_anno_df <- cbind(data.frame(scDEG_fgseaRes_all_sig2), scDEG_fgseaRes_all_sig_anno_df) # combine anno

# write.table(scDEG_fgseaRes_all_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_All_sigAnno.txt', sep = '\t', row.names = FALSE)


############################################################# 3. LINCS on Neutrophil scDEGs, up- and down-regulated separately
## Neutrophil scDEG
scDEG_Neu_up_df <- read.table('Genelist/Genelist3_Neutrophil_scDEGs_PCG.txt') %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% arrange(desc(avg_log2FC)) # 392 (logFC=0.5), 127 (logFC=1)
scDEG_Neu_down_df <- read.table('Genelist/Genelist3_Neutrophil_scDEGs_PCG.txt') %>% filter(p_val_adj < 0.05, avg_log2FC < (-0.5)) %>% arrange(desc(avg_log2FC)) # 1142 (logFC=0.5), 404 (logFC=1)

## genelists
scDEG_Neu_up <- scDEG_Neu_up_df$avg_log2FC
names(scDEG_Neu_up) <- rownames(scDEG_Neu_up_df)
scDEG_Neu_down <- scDEG_Neu_down_df$avg_log2FC
names(scDEG_Neu_down) <- rownames(scDEG_Neu_down_df)

## gsea
cl_cores <- 60
cl <- makeCluster(cl_cores)
registerDoParallel(cl)
scDEG_Neu_fgseaRes_up <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'scDEG_Neu_up') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = scDEG_Neu_up,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}
scDEG_Neu_fgseaRes_down <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'scDEG_Neu_down') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = scDEG_Neu_down,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}
stopCluster(cl) 


## list as vector
scDEG_Neu_fgseaRes_up$leadingEdge <- sapply(scDEG_Neu_fgseaRes_up$leadingEdge, paste, collapse = ", ") # 34234
scDEG_Neu_fgseaRes_down$leadingEdge <- sapply(scDEG_Neu_fgseaRes_down$leadingEdge, paste, collapse = ", ") # 52381

write.table(scDEG_Neu_fgseaRes_up, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_UP_All.txt', sep = '\t', row.names = FALSE)
write.table(scDEG_Neu_fgseaRes_down, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_DOWN_All.txt', sep = '\t', row.names = FALSE)

## significant
scDEG_Neu_fgseaRes_up_sig <- scDEG_Neu_fgseaRes_up %>% filter(padj < 0.05) %>% arrange(desc(NES)) # 1803
scDEG_Neu_fgseaRes_down_sig <- scDEG_Neu_fgseaRes_down %>% filter(padj < 0.05) %>% arrange(desc(NES)) # 3127


## pathway annotation
pert_anno_phase1 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
pert_anno_phase2 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz'), sep = '\t')

scDEG_Neu_fgseaRes_up_sig_anno <- lapply(scDEG_Neu_fgseaRes_up_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }
    
    return(pert_anno)
})
scDEG_Neu_fgseaRes_down_sig_anno <- lapply(scDEG_Neu_fgseaRes_down_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }
    return(pert_anno)
})

scDEG_Neu_fgseaRes_up_sig_anno_df <- Reduce(rbind, scDEG_Neu_fgseaRes_up_sig_anno) # 1376
scDEG_Neu_fgseaRes_down_sig_anno_df <- Reduce(rbind, scDEG_Neu_fgseaRes_down_sig_anno) # 2952

## identical
identical(scDEG_Neu_fgseaRes_up_sig$pathway, scDEG_Neu_fgseaRes_up_sig_anno_df$pert_id)
identical(scDEG_Neu_fgseaRes_down_sig$pathway, scDEG_Neu_fgseaRes_down_sig_anno_df$pert_id)

## match
scDEG_Neu_fgseaRes_up_sig2 <- scDEG_Neu_fgseaRes_up_sig %>% filter(pathway %in% scDEG_Neu_fgseaRes_up_sig_anno_df$pert_id) # 1803, remove pathways not in pert anno
scDEG_Neu_fgseaRes_down_sig2 <- scDEG_Neu_fgseaRes_down_sig %>% filter(pathway %in% scDEG_Neu_fgseaRes_down_sig_anno_df$pert_id) # 3125
identical(scDEG_Neu_fgseaRes_up_sig2$pathway, scDEG_Neu_fgseaRes_up_sig_anno_df$pert_id)
identical(scDEG_Neu_fgseaRes_down_sig2$pathway, scDEG_Neu_fgseaRes_down_sig_anno_df$pert_id)

## combine
scDEG_Neu_fgseaRes_up_sig_anno_df <- cbind(data.frame(scDEG_Neu_fgseaRes_up_sig2), scDEG_Neu_fgseaRes_up_sig_anno_df) %>% arrange(NES) # neg NES to pos NES
scDEG_Neu_fgseaRes_down_sig_anno_df <- cbind(data.frame(scDEG_Neu_fgseaRes_down_sig2), scDEG_Neu_fgseaRes_down_sig_anno_df) %>% arrange(desc(NES)) # pos NES to neg NES

write.table(scDEG_Neu_fgseaRes_up_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_UP_sigAnno.txt', sep = '\t', row.names = FALSE)
write.table(scDEG_Neu_fgseaRes_down_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_DOWN_sigAnno.txt', sep = '\t', row.names = FALSE)




############################################################# 3. LINCS on bulkDEGs, up- and down-regulated
## RNAseq DEG
bulkDEG <- read.table('SLE_DEG_RNASeq/SLE_Con_DEG_new_all.csv', sep = ',', row.names = 1, header = TRUE)
bulkDEG$gene_name <- str_to_title(bulkDEG$gene_name) # To lower case for mouse genes
bulkDEG_up <- bulkDEG %>% filter(FDR<0.05, logFC>1) %>% arrange(desc(logFC)) # 784
bulkDEG_down <- bulkDEG %>% filter(FDR<0.05, logFC<(-1)) %>% arrange(desc(logFC)) # 81

## up and down genelists
bulkDEG_up_list <- bulkDEG_up$logFC
names(bulkDEG_up_list) <- bulkDEG_up$gene_name
bulkDEG_down_list <- bulkDEG_down$logFC
names(bulkDEG_down_list) <- bulkDEG_down$gene_name

## gsea
cl_cores <- 60
cl <- makeCluster(cl_cores)
registerDoParallel(cl)

bulkDEG_fgseaRes_up <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'bulkDEG_up_list') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = bulkDEG_up_list,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}

bulkDEG_fgseaRes_down <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'bulkDEG_down_list') %dopar% {
    prl_gsea <- fgsea(pathways = prl, 
                      stats    = bulkDEG_down_list,
                      minSize  = 10,
                      maxSize  = 500)
    prl_gsea$pathway <- prl$head
    return(prl_gsea)
}
stopCluster(cl) 

## list as vector
bulkDEG_fgseaRes_up$leadingEdge <- sapply(bulkDEG_fgseaRes_up$leadingEdge, paste, collapse = ", ") # 51757
bulkDEG_fgseaRes_down$leadingEdge <- sapply(bulkDEG_fgseaRes_down$leadingEdge, paste, collapse = ", ") # 52381

write.table(bulkDEG_fgseaRes_up, 'Genelist/Genelist3_LINCS_GSEA_bulkDEG_UP_All.txt', sep = '\t', row.names = FALSE)
write.table(bulkDEG_fgseaRes_down, 'Genelist/Genelist3_LINCS_GSEA_bulkDEG_DOWN_All.txt', sep = '\t', row.names = FALSE)

## significant padj<0.05
bulkDEG_fgseaRes_up_sig <- bulkDEG_fgseaRes_up %>% filter(padj<0.05) %>% arrange(NES) # 6963
bulkDEG_fgseaRes_down_sig <- bulkDEG_fgseaRes_down %>% filter(padj<0.05) %>% arrange(desc(NES)) # 0


## pert annotation
pert_anno_phase1 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
pert_anno_phase2 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz'), sep = '\t')
rownames(pert_anno_phase1) <- pert_anno_phase1$pert_id
rownames(pert_anno_phase2) <- pert_anno_phase2$pert_id

bulkDEG_fgseaRes_up_sig_anno <- lapply(bulkDEG_fgseaRes_up_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }
    
    return(pert_anno)
})

bulkDEG_fgseaRes_down_sig_anno <- lapply(bulkDEG_fgseaRes_down_sig$pathway, function(pert){
    if(pert %in% pert_anno_phase1$pert_id){
        pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    } else {
        pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
    }

    return(pert_anno)
})

bulkDEG_fgseaRes_up_sig_anno_df <- Reduce(rbind, bulkDEG_fgseaRes_up_sig_anno) # 6961
bulkDEG_fgseaRes_down_sig_anno_df <- Reduce(rbind, bulkDEG_fgseaRes_down_sig_anno) # 0, no enrichment found


## identical
identical(bulkDEG_fgseaRes_up_sig$pathway, bulkDEG_fgseaRes_up_sig_anno_df$pert_id)
identical(bulkDEG_fgseaRes_down_sig$pathway, bulkDEG_fgseaRes_down_sig_anno_df$pert_id)

## match, only upregulated
bulkDEG_fgseaRes_up_sig2 <- bulkDEG_fgseaRes_up_sig %>% filter(pathway %in% bulkDEG_fgseaRes_up_sig_anno_df$pert_id) # 6961, remove pathways not in pert anno
identical(bulkDEG_fgseaRes_up_sig2$pathway, bulkDEG_fgseaRes_up_sig_anno_df$pert_id)

## combine
bulkDEG_fgseaRes_up_sig_anno_df <- cbind(data.frame(bulkDEG_fgseaRes_up_sig2), bulkDEG_fgseaRes_up_sig_anno_df) %>% arrange(NES) # neg NES to pos NES

write.table(bulkDEG_fgseaRes_up_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_bulkDEG_sigAnno.txt', sep = '\t', row.names = FALSE)

# ############################################################# 4. LINCS on bulkDEGs integrating up- and down-regulated in one list
# ## RNAseq DEG
# bulkDEG <- read.table('SLE_DEG_RNASeq/SLE_Con_DEG_new_all.csv', sep = ',', row.names = 1, header = TRUE)
# bulkDEG <- bulkDEG %>% filter(FDR<0.05, abs(logFC)>1) %>% arrange(desc(logFC)) # 865
# bulkDEG$gene_name <- str_to_title(bulkDEG$gene_name) # To lower case for mouse genes
# bulkDEG_list <- bulkDEG$logFC
# names(bulkDEG_list) <- bulkDEG$gene_name

# ## gsea
# cl_cores <- 60
# cl <- makeCluster(cl_cores)
# registerDoParallel(cl)
# bulkDEG_fgseaRes_all <- foreach(prl=PRL_gmt, .packages = 'fgsea', .combine = 'rbind', .export = 'bulkDEG_list') %dopar% {
#     prl_gsea <- fgsea(pathways = prl, 
#                       stats    = bulkDEG_list,
#                       minSize  = 10,
#                       maxSize  = 500)
#     prl_gsea$pathway <- prl$head
#     return(prl_gsea)
# }
# stopCluster(cl) 

# ## list as vector
# bulkDEG_fgseaRes_all$leadingEdge <- sapply(bulkDEG_fgseaRes_all$leadingEdge, paste, collapse = ", ") # 47754
# write.table(bulkDEG_fgseaRes_all, 'Genelist/Genelist3_LINCS_GSEA_bulkDEG_All.txt', sep = '\t', row.names = FALSE)

# ## significant
# bulkDEG_fgseaRes_all_sig <- bulkDEG_fgseaRes_all %>% filter(padj < 0.05) %>% arrange(desc(NES)) # 7536

# ## pathway annotation
# pert_anno_phase1 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
# pert_anno_phase2 <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz'), sep = '\t')

# bulkDEG_fgseaRes_all_sig_anno <- lapply(bulkDEG_fgseaRes_all_sig$pathway, function(pert){
#     if(pert %in% pert_anno_phase1$pert_id){
#         pert_anno <- pert_anno_phase1 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
#     } else {
#         pert_anno <- pert_anno_phase2 %>% filter(pert_id==pert) %>% select(pert_id, pert_iname, pert_type)
#     }
    
#     return(pert_anno)
# })

# ## keep matched pathways and annotation
# bulkDEG_fgseaRes_all_sig_anno_df <- Reduce(rbind, bulkDEG_fgseaRes_all_sig_anno) # 7534
# bulkDEG_fgseaRes_all_sig2 <- bulkDEG_fgseaRes_all_sig %>% filter(pathway %in% bulkDEG_fgseaRes_all_sig_anno_df$pert_id) # 7534, remove pathways not in pert anno
# identical(bulkDEG_fgseaRes_all_sig2$pathway, bulkDEG_fgseaRes_all_sig_anno_df$pert_id)
# bulkDEG_fgseaRes_all_sig_anno_df <- cbind(data.frame(bulkDEG_fgseaRes_all_sig2), bulkDEG_fgseaRes_all_sig_anno_df) # combine anno

# write.table(bulkDEG_fgseaRes_all_sig_anno_df, 'Genelist/Genelist3_LINCS_GSEA_bulkDEG_All_sigAnno.txt', sep = '\t', row.names = FALSE)
