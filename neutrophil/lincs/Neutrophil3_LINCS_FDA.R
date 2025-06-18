## find LINCS drugs is FDA approved or not

library(tidyverse)

FDA_drugs <- read_tsv('LINCS/FDADrug/Products.txt')


############################ 1. Neu scDEG
## LINCS
## upregulated genes --> NES < 0
## downregulated genes --> NES > 0
LINCS_Neu_up_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_UP_sigAnno.txt', header = 1) %>% arrange(NES) %>% filter(NES < 0) # 360
LINCS_Neu_down_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_DOWN_sigAnno.txt', header = 1) %>% arrange(desc(NES)) %>% filter(NES > 0) # 123
## to upper for FDA matching
LINCS_Neu_up_sig$pert_iname <- toupper(LINCS_Neu_up_sig$pert_iname)
LINCS_Neu_down_sig$pert_iname <- toupper(LINCS_Neu_down_sig$pert_iname)
## filter BRD compounds perturbations
LINCS_Neu_up_sig <- LINCS_Neu_up_sig %>% filter(!startsWith(pert_iname, 'BRD')) # 263
LINCS_Neu_down_sig <- LINCS_Neu_down_sig %>% filter(!startsWith(pert_iname, 'BRD')) # 81
## only keep compounds and Peptides and other biological agents (e.g. cytokine)
LINCS_Neu_up_sig <- LINCS_Neu_up_sig %>% filter(pert_type %in% c('trt_cp','trt_lig')) # 36
LINCS_Neu_down_sig <- LINCS_Neu_down_sig %>% filter(pert_type %in% c('trt_cp','trt_lig')) # 11

## FDA
FDA_LINCS_Neu_up <- FDA_drugs %>% filter(DrugName %in% LINCS_Neu_up_sig$pert_iname | ActiveIngredient %in% LINCS_Neu_up_sig$pert_iname) # 13 drugNames with 6 active ingredients
FDA_LINCS_Neu_down <- FDA_drugs %>% filter(DrugName %in% LINCS_Neu_down_sig$pert_iname | ActiveIngredient %in% LINCS_Neu_down_sig$pert_iname) # 0

FDA_LINCS_Neu_up %>% select(DrugName, ActiveIngredient) %>% unique()
FDA_LINCS_Neu_down %>% select(DrugName, ActiveIngredient) %>% unique()

## FDA labeling
LINCS_Neu_up_sig2 <- LINCS_Neu_up_sig %>% mutate(FDA = ifelse(pert_iname %in% FDA_drugs$DrugName | pert_iname %in% FDA_drugs$ActiveIngredient, 'Approved', 'not Approved'))
LINCS_Neu_down_sig2 <- LINCS_Neu_down_sig %>% mutate(FDA = ifelse(pert_iname %in% FDA_drugs$DrugName | pert_iname %in% FDA_drugs$ActiveIngredient, 'Approved', 'not Approved'))

write.table(LINCS_Neu_up_sig2, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_UP_sigAnno_Selected.txt')
write.table(LINCS_Neu_down_sig2, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Neutrophil_DOWN_sigAnno_Selected.txt')

############################ 2. Immune scDEG
## LINCS immune scDEGs
LINCS_immune_up_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_UP_sigAnno.txt', header=1) %>% arrange(NES) %>% filter(NES < 0) # 1794
LINCS_immune_down_sig <- read.table('Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_DOWN_sigAnno.txt', header=1) %>% arrange(desc(NES)) %>% filter(NES > 0) # 211
## to upper for FDA matching
LINCS_immune_up_sig$pert_iname <- toupper(LINCS_immune_up_sig$pert_iname)
LINCS_immune_down_sig$pert_iname <- toupper(LINCS_immune_down_sig$pert_iname)
## filter BRD compounds perturbations
LINCS_immune_up_sig <- LINCS_immune_up_sig %>% filter(!startsWith(pert_iname, 'BRD')) # 937
LINCS_immune_down_sig <- LINCS_immune_down_sig %>% filter(!startsWith(pert_iname, 'BRD')) # 155
## only keep compounds and Peptides and other biological agents (e.g. cytokine)
LINCS_immune_up_sig <- LINCS_immune_up_sig %>% filter(pert_type %in% c('trt_cp','trt_lig')) # 151
LINCS_immune_down_sig <- LINCS_immune_down_sig %>% filter(pert_type %in% c('trt_cp','trt_lig')) # 17

## FDA
FDA_LINCS_immune_up <- FDA_drugs %>% filter(DrugName %in% LINCS_immune_up_sig$pert_iname | ActiveIngredient %in% LINCS_immune_up_sig$pert_iname) # 55 drugNames with 19 Active ingredients
FDA_LINCS_immune_down <- FDA_drugs %>% filter(DrugName %in% LINCS_immune_down_sig$pert_iname | ActiveIngredient %in% LINCS_immune_down_sig$pert_iname) # 2 drugNames with 2 Active ingredients

FDA_LINCS_immune_up %>% select(DrugName, ActiveIngredient) %>% unique()
FDA_LINCS_immune_down %>% select(DrugName, ActiveIngredient) %>% unique()

## FDA labeling
LINCS_immune_up_sig2 <- LINCS_immune_up_sig %>% mutate(FDA = ifelse(pert_iname %in% FDA_drugs$DrugName | pert_iname %in% FDA_drugs$ActiveIngredient, 'Approved', 'not Approved'))
LINCS_immune_down_sig2 <- LINCS_immune_down_sig %>% mutate(FDA = ifelse(pert_iname %in% FDA_drugs$DrugName | pert_iname %in% FDA_drugs$ActiveIngredient, 'Approved', 'not Approved'))

write.table(LINCS_immune_up_sig2, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_UP_sigAnno_Selected.txt')
write.table(LINCS_immune_down_sig2, 'Genelist/Genelist3_LINCS_GSEA_scDEG_Immune_DOWN_sigAnno_Selected.txt')
