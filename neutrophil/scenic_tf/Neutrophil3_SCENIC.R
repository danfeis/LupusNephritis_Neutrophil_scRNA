## add SCENIC to Neutrophil Seurat Object

library(tidyverse)

load('RData/Immune.combined3_SCENIC.RData')
load('RData/Neutrophil3.RData')

## binarized regulons for immune.combined3
AUC_mtx <- immune.combined3@meta.data %>% select(starts_with("Regulon")) # 69452 * 297

## subset regulons for neutrophils
cells_neu <- rownames(neutrophil@meta.data)
AUC_mtx_neu <- AUC_mtx[cells_neu,] # 14618 * 297

## combined regulons to neutrophil Seurat object
identical(cells_neu, rownames(AUC_mtx_neu))
neutrophil@meta.data <- cbind(neutrophil@meta.data, AUC_mtx_neu)

# save(neutrophil, file = 'RData/Neutrophil3_SCENIC.RData')

############################################## subset TF regulons enriched in neutrophil
library(tidyverse)
library(Seurat)
library(pheatmap)

load('RData/Neutrophil3_SCENIC.RData')

## filter TF regulons activated less than 5%
AUC_mtx_neu <- neutrophil@meta.data %>% select(starts_with("Regulon"))
AUC_mtx_neu_selected <- AUC_mtx_neu[,colSums(AUC_mtx_neu) > nrow(AUC_mtx_neu)*0.05] # 297 --> 145

## cell annotation
anno_df <- neutrophil@meta.data %>% select(mature_clusters, mature_groups, type)
colnames(anno_df) <- c('MatureClusters','MatureGroups','DiseaseStages')
my_color_diseasestages <- c('#eed2cc','#e8998d','#6c9a8b')
my_color_maturegroup <- c('#f8eea4','#a8d1d1','#d66b78','#7b7fc4')
my_color_matureclusters <- c('#f8eea4','#ecd416','#a8d1d1','#4d9393','#d66b78','#982a37','#c1c3e4', '#7b7fc4', '#363a7a','#181a36', '#627a9d')
names(my_color_diseasestages) <- levels(neutrophil$type)
names(my_color_maturegroup) <- levels(neutrophil$mature_groups)
names(my_color_matureclusters) <- levels(neutrophil$mature_clusters)
anno_colors <- list(MatureClusters = my_color_matureclusters,
                    MatureGroups = my_color_maturegroup,
                    DiseaseStages = my_color_diseasestages)

p <- pheatmap(t(AUC_mtx_neu_selected), 
              show_colnames = F,
              treeheight_col = 0,
              cutree_rows = 6,
              fontsize = 12,
              fontsize_row = 8,
              annotation_col = anno_df,
              annotation_colors = anno_colors)
ggsave('test.pdf', p, width = 10, height = 25)


## TF activity percentage in each maturation groups
TFActivity_freq <- neutrophil@meta.data %>% 
                        select(mature_groups, colnames(AUC_mtx_neu_selected)) %>% 
                        pivot_longer(cols = colnames(AUC_mtx_neu_selected), names_to = 'TF', values_to = 'TFActivity') %>%
                        mutate(TF = str_remove(TF, "^Regulon_")) %>% # remove 'Regulon_' in TF strings
                        group_by(mature_groups, TF) %>%
                        summarise(freq = sum(TFActivity)/n()) %>%
                        ungroup()

TFActivity_freq2 <- TFActivity_freq %>% pivot_wider(names_from = TF, values_from = freq, values_fill = 0) %>% select(-1) %>% t() %>% data.frame()
colnames(TFActivity_freq2) <- levels(neutrophil$mature_groups)

p <- pheatmap(TFActivity_freq2,cutree_rows = 6, cluster_cols = FALSE)
ggsave('NeuTF_Freq.pdf', p, width = 10, height = 20)

## shared TFs or TFs specific to maturation groups, cutoff set to 0.5 and 0.9
TFActivity_freq_select_common <- TFActivity_freq2 %>% filter(Imm_G1>0.9, Imm_G2>0.9, Mat_G3>0.9, Mat_G4>0.9) # 12, all > 0.9
TFActivity_freq_select_Mat <- TFActivity_freq2 %>% filter(Imm_G1<0.6, Imm_G2<0.6, Mat_G3>0.9, Mat_G4>0.9) # 7, Mat > 0.9, and Imm < 0.6. do not set this too low as the Imm_G1 may have many of these Tfs activated aroung 0.4 and 0.5
TFActivity_freq_select_ImmG1 <- TFActivity_freq2 %>% filter(Imm_G1>0.9, Imm_G2<0.6, Mat_G3<0.6, Mat_G4<0.6) # 4, Imm_G1 > 0.9, others < 0.6
TFActivity_freq_select_ImmG2 <- TFActivity_freq2 %>% filter(Imm_G1<0.3, Imm_G2>0.5, Mat_G3<0.3, Mat_G4<0.3) # 22, lower activation frequency in Imm_G2 compared with other maturation groups.

write.table(TFActivity_freq_select_ImmG1, file = 'Genelist/Genelist3_TFActivityFraction_SCENIC_ImmG1.txt', quote = FALSE)
write.table(TFActivity_freq_select_ImmG2, file = 'Genelist/Genelist3_TFActivityFraction_SCENIC_ImmG2.txt', quote = FALSE)
write.table(TFActivity_freq_select_Mat, file = 'Genelist/Genelist3_TFActivityFraction_SCENIC_Mat.txt', quote = FALSE)
write.table(TFActivity_freq_select_common, file = 'Genelist/Genelist3_TFActivityFraction_SCENIC_Common.txt', quote = FALSE)

## combined TF
TFlist_all <- rbind(data.frame(TF = rownames(TFActivity_freq_select_ImmG1), TFGroup = 'Imm_G1'),
                    data.frame(TF = rownames(TFActivity_freq_select_ImmG2), TFGroup = 'Imm_G2'),
                    data.frame(TF = rownames(TFActivity_freq_select_Mat), TFGroup = 'Mat'),
                    data.frame(TF = rownames(TFActivity_freq_select_Common), TFGroup = 'Common'))

write.table(TFlist_all, file = 'Genelist/Genelist3_Neutrophil_TFlist_SCENIC_All.txt', quote = FALSE, row.names = FALSE)
