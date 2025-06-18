## summary of counts of cell lines, perts

## files of LINCS phase 1 and 2
LINCS_phase1_sigdf <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_sig_info.txt.gz'), sep = '\t') # Metadata for each signature in the Level 5 matrix
LINCS_phase2_sigdf <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz'), sep = '\t')
LINCS_phase1_pert <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
LINCS_phase2_pert <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_pert_info.txt.gz'), sep = '\t')
LINCS_phase1_cellid <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase1/GSE92742_Broad_LINCS_cell_info.txt.gz'), sep = '\t')
LINCS_phase2_cellid <- read.csv(gzfile('/public/home/zhudf/LINCS/Dataset/LINCS_Phase2/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz'), sep = '\t')

## pert
LINCS_phase1_pert_trt <- LINCS_phase1_pert %>%  filter(startsWith(pert_type, 'trt')) # 51383-->51307
LINCS_phase2_pert_trt <- LINCS_phase2_pert %>%  filter(startsWith(pert_type, 'trt')) # 2170-->2149

## only trt kept
LINCS_phase1_sigdf_trt <- LINCS_phase1_sigdf %>% filter(startsWith(pert_type, 'trt')) # 473647-->451582
LINCS_phase2_sigdf_trt <- LINCS_phase2_sigdf %>% filter(startsWith(pert_type, 'trt')) # 118050-->111286

##-------------------------------------------------- cell lines
## cell lines in trt kept
celllines_phase1 <- unique(LINCS_phase1_sigdf_trt$cell_id) # 76 (98 in total from LINCS_phase1_cellid)
celllines_phase2 <- unique(LINCS_phase2_sigdf_trt$cell_id) # 41 (98 in total from LINCS_phase2_cellid)
unique(LINCS_phase1_cellid$cell_id)
unique(LINCS_phase2_cellid$cell_id)
unique(LINCS_phase1_cellid$base_cell_id)
unique(LINCS_phase2_cellid$base_cell_id)

## check cell_id or cell_base_line
LINCS_phase1_sigdf2 <- LINCS_phase1_sigdf %>% mutate(cell_id2 = str_extract(sig_id, "^[^:]+")) %>% mutate(cell_id2 = str_split_fixed(cell_id2, "_", 3)[, 2]) # add cell_id check if A375.311 exsited
LINCS_phase2_sigdf2 <- LINCS_phase2_sigdf %>% mutate(cell_id2 = str_extract(sig_id, "^[^:]+")) %>% mutate(cell_id2 = str_split_fixed(cell_id2, "_", 3)[, 2]) 
LINCS_phase1_sigdf_trt2 <- LINCS_phase1_sigdf_trt %>% mutate(cell_id2 = str_extract(sig_id, "^[^:]+")) %>% mutate(cell_id2 = str_split_fixed(cell_id2, "_", 3)[, 2]) # add cell_id from sig_id
LINCS_phase2_sigdf_trt2 <- LINCS_phase2_sigdf_trt %>% mutate(cell_id2 = str_extract(sig_id, "^[^:]+")) %>% mutate(cell_id2 = str_split_fixed(cell_id2, "_", 3)[, 2])
LINCS_phase1_sigdf_trt2 %>% select(cell_id, cell_id2) %>% unique()
LINCS_phase2_sigdf_trt2 %>% select(cell_id, cell_id2) %>% unique()
unique(LINCS_phase1_sigdf_trt2$cell_id)
unique(LINCS_phase1_sigdf_trt2$cell_id2)
unique(LINCS_phase2_sigdf_trt2$cell_id)
unique(LINCS_phase2_sigdf_trt2$cell_id2)

## check cell_id in phase1 identical to cell line from sig_id
identical(LINCS_phase1_sigdf_trt2$cell_id, LINCS_phase1_sigdf_trt2$cell_id2) # identical, and no A375.311 (cas9 modified) in Phase 1

## cell lines diff
cellline_phase1 <- unique(LINCS_phase1_sigdf_trt2$cell_id2)
cellline_phase2 <- unique(LINCS_phase2_sigdf_trt2$cell_id2)
intersect(cellline_phase1, cellline_phase2)
setdiff(cellline_phase1, cellline_phase2)
setdiff(cellline_phase2, cellline_phase1)

## cell lines combine
cellline_all <- union(LINCS_phase1_sigdf_trt2$cell_id2, LINCS_phase2_sigdf_trt2$cell_id2) # 98

##-------------------------------------------------- pert
pert_phase1 <- unique(LINCS_phase1_sigdf_trt$pert_id)
pert_phase2 <- unique(LINCS_phase2_sigdf_trt$pert_id)
intersect(pert_phase1, pert_phase2) # 910
setdiff(pert_phase1, pert_phase2) # 50232
setdiff(pert_phase2, pert_phase1) # 1239

pert_all <- union(pert_phase1, pert_phase2) # 52381 