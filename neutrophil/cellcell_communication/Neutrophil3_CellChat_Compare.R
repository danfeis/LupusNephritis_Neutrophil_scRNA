## compare L-R signaling in mild and severe stages

library(Seurat)
library(pheatmap)
library(network)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(CellChat)
library(future)
library(reshape2)
library(cowplot)
library(stringr)


load('RData/Immune.combined3_CellChat.RData')
load('RData/Immune.combined2_CellChat2.RData')

## significant pathways
pathways_netP_mild <- cellchat_mild@netP$pathways # 135 pathways and 720 LR pairs in LR slot, 54 significant pathways in netP
pathways_netP_severe <- cellchat_severe@netP$pathways # 151 pathways, 833 LR pairs in LR slot, 58 significant pathways in netP

setdiff(pathways_netP_mild, pathways_netP_severe) # NULL
setdiff(pathways_netP_severe, pathways_netP_mild) # SPP1, CD23, VISFATIN, PD-L1

## LR specific to mild or severe stages
LR_mild <- cellchat_mild@LR$LRsig %>% filter(pathway_name %in% pathways_netP_mild) %>% select(interaction_name, pathway_name, ligand, receptor) # 338 LR
LR_severe <- cellchat_severe@LR$LRsig %>% filter(pathway_name %in% pathways_netP_severe) %>% select(interaction_name, pathway_name, ligand, receptor) # 368 LR

unique_LR_mild <- setdiff(LR_mild$interaction_name, LR_severe$interaction_name) # LR specific to mild
unique_LR_severe <- setdiff(LR_severe$interaction_name, LR_mild$interaction_name) # LR specific to severe

mild_unique <- LR_mild[unique_LR_mild, , drop = FALSE] # signaling include CCL, MIF, PARS, MHC-I
severe_unique <- LR_severe[unique_LR_severe, , drop = FALSE] # CCL, CXCL, IL2, IL1, SPP1, VISFATIN, LAMININ, THBS, CD23, CD39, CDH, ICAM, PD-L1

## selected signaling
signaling_selected <- c('CCL','CXCL','MIF','IL2','CD23','SPP1','VISFATIN','PD-L1')
LR_selected <- rbind(mild_unique %>% filter(pathway_name %in% signaling_selected) %>% mutate(Condition = 'Mild'),
                     severe_unique %>% filter(pathway_name %in% signaling_selected) %>% mutate(Condition = 'Severe'))
LR_selected_molecules <-  LR_selected %>% 
                            separate(receptor, into = c('receptor1', 'receptor2'), sep = '_') %>% # separate if there are more than 1 receptor like ITGA9_ITGB1
                            mutate(receptor1 = str_to_title(str_to_lower(receptor1)), 
                                   receptor2 = str_to_title(str_to_lower(receptor2))) %>% # capitalize
                            pivot_longer(cols = c(receptor1, receptor2), names_to = NULL, values_to = 'receptor') %>%
                            drop_na(receptor) %>% # remove rows if receptor is NA
                            pivot_longer(cols = c(ligand, receptor), names_to = 'LR', values_to = 'Genes') %>%
                            arrange(pathway_name, LR)
LR_selected_genes <- unique(LR_selected_molecules$Genes) # all ligands and receptors

write.table(LR_selected_genes,'test.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)


netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object





############################################### Identify dysfunctional signaling (up- and down-regulated) by using differential expression analysis, from CellChat Tutorial
object.list <- list(mild = cellchat_mild, severe = cellchat_severe)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

c('DNT','CD4+T','CD8+T','B','NK','Basophil','Monocyte','pDC','Megaka','Erythroid')
c('Neutrophil_Imm_G1','Neutrophil_Imm_G2','Neutrophil_Mat_G3','Neutrophil_Mat_G4')

netVisual_bubble(cellchat,
                 sources.use = c('Neutrophil_Imm_G1','Neutrophil_Imm_G2','Neutrophil_Mat_G3','Neutrophil_Mat_G4'),
                 targets.use = c('DNT','CD4+T','CD8+T','B','NK','Basophil','Monocyte','pDC','Megaka','Erythroid'),
                 comparison = c(1, 2),
                 max.dataset = 2,
                 title.name = "Increased signaling in CCA", angle.x = 45, remove.isolate = T)

netVisual_bubble(cellchat,
                 sources.use = c('DNT','CD4+T','CD8+T','B','NK','Basophil','Monocyte','pDC','Megaka','Erythroid'),
                 targets.use = c('Neutrophil_Imm_G1','Neutrophil_Imm_G2','Neutrophil_Mat_G3','Neutrophil_Mat_G4'),
                 comparison = c(1, 2),
                 max.dataset = 2,
                 title.name = "Increased signaling in CCA", angle.x = 45, remove.isolate = T)

## differential signaling
pos.dataset = "severe" # the dataset with positive fold change
features.name = paste0(pos.dataset, ".merged") # define a char name used for storing the results of differential expression analysis
cellchat2 <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05) 
net <- netMappingDEG(cellchat2, features.name = features.name)

## signaling with upregulated and downregulated ligands
ligand_net.up <- subsetCommunication(cellchat2, net = net, datasets = "severe",ligand.logFC = 0.1, receptor.logFC = NULL) # 667, 37 signaling (logFC=0.05) --> 613, 33 signaling (logFC=0.1)
ligand_net.down <- subsetCommunication(cellchat2, net = net, datasets = "mild",ligand.logFC = -0.1, receptor.logFC = NULL) # 397, 21 signaling (logFC=0.05) --> 302, 17 signaling (logFC=0.1)

## signaling with upregulated and downregulated receptors
receptor_net.up <- subsetCommunication(cellchat2, net = net, datasets = "severe",ligand.logFC = NULL, receptor.logFC = 0.1) # 939, 45 signaling (logFC=0.05) --> 832, 43 signaling (logFC=0.1)
receptor_net.down <- subsetCommunication(cellchat2, net = net, datasets = "mild",ligand.logFC = NULL, receptor.logFC = -0.1) # 127, 20 signaling (logFC=0.05) --> 117, 16 signaling (logFC=0.1)

write_tsv(ligand_net.up, file = 'Genelist/Genelist3_CellChat_AllDiffSignaling_Ligand_Upregulate.tsv')
write_tsv(ligand_net.down, file = 'Genelist/Genelist3_CellChat_AllDiffSignaling_Ligand_Downregulate.tsv')
write_tsv(receptor_net.up, file = 'Genelist/Genelist3_CellChat_AllDiffSignaling_Receptor_Upregulate.tsv')
write_tsv(receptor_net.down, file = 'Genelist/Genelist3_CellChat_AllDiffSignaling_Receptor_Downregulate.tsv')


########################################################## Neutrophils associated differential signaling
## differential signaling with Neutrophils
Neu_clusters <- c('Neutrophil_Imm_G1','Neutrophil_Imm_G2','Neutrophil_Mat_G3','Neutrophil_Mat_G4')

## signaling with upregulated and downregulated ligands
Neu_ligand_net.up <- ligand_net.up %>% filter(source %in% Neu_clusters | target %in% Neu_clusters) %>% group_by(source) %>% arrange(target) %>% ungroup() # 297
Neu_ligand_net.down <- ligand_net.down %>% filter(source %in% Neu_clusters | target %in% Neu_clusters) %>% group_by(source) %>% arrange(target) %>% ungroup() # 133

## signaling with upregulated and downregulated receptors
Neu_receptor_net.up <- receptor_net.up %>% filter(source %in% Neu_clusters | target %in% Neu_clusters) %>% group_by(source) %>% arrange(target) %>% ungroup() # 412
Neu_receptor_net.down  <- receptor_net.down %>% filter(source %in% Neu_clusters | target %in% Neu_clusters) %>% group_by(source) %>% arrange(target) %>% ungroup() # 51

## combine
Neu_diffsignal_all <- rbind(Neu_ligand_net.up %>% select(source, target, ligand, receptor, pathway_name) %>% mutate(Regulate = 'Up', LRDiffSource = 'Ligand'),
                            Neu_ligand_net.down %>% select(source, target, ligand, receptor, pathway_name) %>% mutate(Regulate = 'Down', LRDiffSource = 'Ligand'),
                            Neu_receptor_net.up %>% select(source, target, ligand, receptor, pathway_name) %>% mutate(Regulate = 'Up', LRDiffSource = 'Receptor'),
                            Neu_receptor_net.down %>% select(source, target, ligand, receptor, pathway_name) %>% mutate(Regulate = 'Down', LRDiffSource = 'Receptor'))


Neu_net_selected <- rbind(Neu_diffsignal_all %>% filter(Regulate == 'Up', pathway_name %in% c('CCL','CXCL','COMPLEMENT','SPP1','VISFATIN','PD-L1','CD23')),
                          Neu_diffsignal_all %>% filter(Regulate == 'Down', pathway_name %in% c('CCL','CXCL', 'IL2', 'LCK')) %>% mutate(Regulate = 'Down'))
Neu_net_selected <- Neu_net_selected %>% 
                        select(pathway_name, ligand, receptor, Regulate) %>% 
                        unique() %>%
                        pivot_longer(cols = c(ligand, receptor), names_to = 'LR', values_to = 'molecules') %>%
                        group_by(pathway_name) %>%
                        arrange(pathway_name, LR)


Neu_net_selected %>% filter(Regulate=='Down', pathway_name=='CCL')
Neu_net_selected %>% filter(Regulate=='Down', pathway_name=='CXCL')
Neu_net_selected %>% filter(Regulate=='Down', pathway_name=='LCK')
