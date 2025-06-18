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


load('RData/Immune.combined3_Neutro.RData')

CellChatDB <- CellChatDB.mouse
future::plan(multisession, workers = 30)

my_level <- c('DNT','CD4+T','CD8+T','B','NK','Basophil','Neutrophil_Imm_G1','Neutrophil_Imm_G2','Neutrophil_Mat_G3','Neutrophil_Mat_G4','Monocyte','pDC','Megaka','Erythroid')
immune.combined3$Neutro_MatureGroups <- factor(immune.combined3$Neutro_MatureGroups, levels = my_level)
Idents(immune.combined3) <- 'Neutro_MatureGroups'


## mild and severe stages
cells_mild <- immune.combined3@meta.data[immune.combined3@meta.data$type=='mild',]
cells_severe <- immune.combined3@meta.data[immune.combined3@meta.data$type=='severe',]
immune_mild <- subset(immune.combined3, cells = rownames(cells_mild)) # 24877
immune_severe <- subset(immune.combined3, cells = rownames(cells_severe)) # 21300



##--------------------------------------------------------------------- 1. cellchat in mild and severe stages
## 1. immune, L-R count in mild stages
expr <- immune_mild@assays$RNA@data
meta <- immune_mild@meta.data
cellchat_mild <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_mild@DB <- CellChatDB
cellchat_mild <- subsetData(cellchat_mild)
cellchat_mild <- identifyOverExpressedGenes(cellchat_mild)
cellchat_mild <- identifyOverExpressedInteractions(cellchat_mild)
cellchat_mild <- computeCommunProb(cellchat_mild)
cellchat_mild <- filterCommunication(cellchat_mild, min.cells = 10)
cellchat_mild <- computeCommunProbPathway(cellchat_mild) # signaling pathway
cellchat_mild <- aggregateNet(cellchat_mild)


## 1.2 immune, L-R count in severe stages
expr <- immune_severe@assays$RNA@data
meta <- immune_severe@meta.data
cellchat_severe <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_severe@DB <- CellChatDB
cellchat_severe <- subsetData(cellchat_severe)
cellchat_severe <- identifyOverExpressedGenes(cellchat_severe)
cellchat_severe <- identifyOverExpressedInteractions(cellchat_severe)
cellchat_severe <- computeCommunProb(cellchat_severe)
cellchat_severe <- filterCommunication(cellchat_severe, min.cells = 10)
cellchat_severe <- computeCommunProbPathway(cellchat_severe) # signaling pathway
cellchat_severe <- aggregateNet(cellchat_severe)

save(cellchat_mild, cellchat_severe, file = 'RData/Immune.combined3_CellChat.RData')


##--------------------------------------------------------------------- 2. secreted signaling, mild and severe stages
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

## mild stage
expr <- immune_mild@assays$RNA@data
meta <- immune_mild@meta.data
cellchat_mild_secreted <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_mild_secreted@DB <- CellChatDB.use
cellchat_mild_secreted <- subsetData(cellchat_mild_secreted)
cellchat_mild_secreted <- identifyOverExpressedGenes(cellchat_mild_secreted)
cellchat_mild_secreted <- identifyOverExpressedInteractions(cellchat_mild_secreted)
cellchat_mild_secreted <- computeCommunProb(cellchat_mild_secreted)
cellchat_mild_secreted <- filterCommunication(cellchat_mild_secreted, min.cells = 10)
cellchat_mild_secreted <- computeCommunProbPathway(cellchat_mild_secreted)
cellchat_mild_secreted <- aggregateNet(cellchat_mild_secreted)


## severe stage
expr <- immune_severe@assays$RNA@data
meta <- immune_severe@meta.data
cellchat_severe_secreted <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_severe_secreted@DB <- CellChatDB.use
cellchat_severe_secreted <- subsetData(cellchat_severe_secreted)
cellchat_severe_secreted <- identifyOverExpressedGenes(cellchat_severe_secreted)
cellchat_severe_secreted <- identifyOverExpressedInteractions(cellchat_severe_secreted)
cellchat_severe_secreted <- computeCommunProb(cellchat_severe_secreted)
cellchat_severe_secreted <- filterCommunication(cellchat_severe_secreted, min.cells = 10)
cellchat_severe_secreted <- computeCommunProbPathway(cellchat_severe_secreted)
cellchat_severe_secreted <- aggregateNet(cellchat_severe_secreted)



##--------------------------------------------------------------------- 3. cell-cell contact signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")

## mild stage count
expr <- immune_mild@assays$RNA@data
meta <- immune_mild@meta.data
cellchat_mild_contact <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_mild_contact@DB <- CellChatDB.use
cellchat_mild_contact <- subsetData(cellchat_mild_contact)
cellchat_mild_contact <- identifyOverExpressedGenes(cellchat_mild_contact)
cellchat_mild_contact <- identifyOverExpressedInteractions(cellchat_mild_contact)
cellchat_mild_contact <- computeCommunProb(cellchat_mild_contact)
cellchat_mild_contact <- filterCommunication(cellchat_mild_contact, min.cells = 10)
cellchat_mild_contact <- computeCommunProbPathway(cellchat_mild_contact)
cellchat_mild_contact <- aggregateNet(cellchat_mild_contact)


## severe stage count
expr <- immune_severe@assays$RNA@data
meta <- immune_severe@meta.data
cellchat_severe_contact <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_severe_contact@DB <- CellChatDB.use
cellchat_severe_contact <- subsetData(cellchat_severe_contact)
cellchat_severe_contact <- identifyOverExpressedGenes(cellchat_severe_contact)
cellchat_severe_contact <- identifyOverExpressedInteractions(cellchat_severe_contact)
cellchat_severe_contact <- computeCommunProb(cellchat_severe_contact)
cellchat_severe_contact <- filterCommunication(cellchat_severe_contact, min.cells = 10)
cellchat_severe_contact <- computeCommunProbPathway(cellchat_severe_contact)
cellchat_severe_contact <- aggregateNet(cellchat_severe_contact)



##--------------------------------------------------------------------- 4. signaling in complete immune object
expr <- immune.combined3@assays$RNA@data
meta <- immune.combined3@meta.data
cellchat_all <- createCellChat(object = expr, meta = meta, group.by = "Neutro_MatureGroups")
cellchat_all@DB <- CellChatDB
cellchat_all <- subsetData(cellchat_all)
cellchat_all <- identifyOverExpressedGenes(cellchat_all)
cellchat_all <- identifyOverExpressedInteractions(cellchat_all)
cellchat_all <- computeCommunProb(cellchat_all)
cellchat_all <- filterCommunication(cellchat_all, min.cells = 10)
cellchat_all <- computeCommunProbPathway(cellchat_all) # signaling pathway
cellchat_all <- aggregateNet(cellchat_all)



save(cellchat_all, cellchat_mild,cellchat_severe,cellchat_mild_secreted,cellchat_severe_secreted,cellchat_mild_contact,cellchat_severe_contact,
     file = 'RData/Immune.combined3_CellChat.RData')