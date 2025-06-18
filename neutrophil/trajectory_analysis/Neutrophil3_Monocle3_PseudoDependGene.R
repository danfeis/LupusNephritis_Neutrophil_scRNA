## identify genes that pseudotime dependent

# conda activate R_v3

library(monocle3)
library(ggplot3)

load('RData/Neutrophil3_Monocle3.RData')

TF_list <- read.table('Genelist/Genelist3_TFlist_SCENIC_All.txt', header=1)

## 1. pseudotime-dependent TFs
# morans_I	Moran's I statistic (how correlated gene expression is along pseudotime)
# morans_test_statistic	Test statistic for evaluating significance
# p_value	p-value for significance of Moran's I
cds_subset_up <- neutrophil.cds[rownames(neutrophil.cds) %in% TF_list$TF,]
ciliated_cds_pr_test_res_TF <- graph_test(cds_subset_up, neighbor_graph="principal_graph", cores=20)
ciliated_cds_pr_test_res_TF_pos <- ciliated_cds_pr_test_res_TF %>% filter(q_value<0.05, morans_I > 0.1) %>% arrange(desc(morans_I)) # 45-->19
write.table(ciliated_cds_pr_test_res_TF_pos, file = 'Genelist/Genelist3_Monocle3_PseudoDependTFs.txt', quote = FALSE)

## 2. pseudotime-dependent genes
ciliated_cds_pr_test_res_genes <- graph_test(neutrophil.cds, neighbor_graph="principal_graph", cores=20)
ciliated_cds_pr_test_res_genes_pos <- ciliated_cds_pr_test_res_genes %>% filter(q_value<0.05, morans_I > 0.1) %>% arrange(desc(morans_I)) # 18374-->831
write.table(ciliated_cds_pr_test_res_genes_pos, file = 'Genelist/Genelist3_Monocle3_PseudoDependGenes.txt', quote = FALSE)

## 3. pseudotime-dependent scDEGs in Neu
Neu_scDEGs_all <- read.table('Genelist/Genelist3_Neutrophil_scDEGs_PCG.txt')
Neu_scDEGs_sig <- Neu_scDEGs_all %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>% mutate(Regulate = ifelse(avg_log2FC < (-1), 'Down', 'Up')) # 531
cds_subset_up <- neutrophil.cds[rownames(neutrophil.cds) %in% rownames(Neu_scDEGs_sig),]
ciliated_cds_pr_test_res_scDEG <- graph_test(cds_subset_up, neighbor_graph="principal_graph", cores=20)
ciliated_cds_pr_test_res_scDEG_pos <- ciliated_cds_pr_test_res_scDEG %>% filter(q_value<0.05, morans_I > 0.1) %>% arrange(desc(morans_I)) # 531-->138
write.table(ciliated_cds_pr_test_res_scDEG_pos, file = 'Genelist/Genelist3_Monocle3_PseudoDependNeuscDEGs.txt', quote = FALSE)


## 4. TF targets that are also pseudotime-dependent scDEGs in Neu
## 4.1 TFs are pseudotime-dependent
## 4.2 targets of these TFs are pseudotime-dependent
## 4.3 targets are also differential expressed as LN progressed.
TF_target <- read_tsv('Genelist/Genelist_SCENIC_Targetgene.txt')
pseudo_depend_TF <- rownames(read.table('Genelist/Genelist3_Monocle3_PseudoDependTFs.txt'))
pseudo_depend_scDEGs <- rownames(read.table('Genelist/Genelist3_Monocle3_PseudoDependNeuscDEGs.txt'))

pseudo_depend_TF_target <- TF_target %>% 
        filter(TF %in% pseudo_depend_TF, TargetGene %in% pseudo_depend_scDEGs, TargetGeneCoef > 0.1) %>% 
        select(TF, TargetGene) %>% 
        unique() # 19 TFs with 368 pseudotime-dependent scDEG
## TF modules that pseudotime-dependent TFs lie in
TF_list %>% filter(TF %in% pseudo_depend_TF) # 11 TF_M2 (Imm_G2), 2 TF_M3 (Mat), 6 TF_M4 (Common)
pseudo_depend_TF_target %>% filter(TF=='Stat1') %>% data.frame()
pseudo_depend_TF_target %>% filter(TF=='Ets2') %>% data.frame()
pseudo_depend_TF_target %>% filter(TF=='Runx3') %>% data.frame()