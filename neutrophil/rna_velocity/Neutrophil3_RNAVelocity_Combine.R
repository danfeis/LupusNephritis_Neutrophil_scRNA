#### RNA Velocity
#### #### combined loom file to Seurat Object


library(velocyto.R)
library(dplyr)
library(Seurat)

load('RData/Neutrophil3.RData')
source('/public/home/zhudf/scRNA_kidney_mouse/RNAVelocity/RNAVelocity_Script/RNAVelocity_Func_CombineLoom.R')

velocity_dirs <- '/public/home/zhudf/scRNA_blood_mouse/CellRanger/'
velocity_samples <- c('MNLNBlood3181_cellranger','MNLNBlood4186_cellranger','MMLNBlood2182_cellranger','MMLNBlood41810_cellranger','MSLNBlood7182_cellranger','MSLNBlood4187_cellranger')
neutrophil_origident <- levels(neutrophil@meta.data$orig.ident) # ('MNL3181BL', 'MNL4186BL', 'MML2182BL', 'MML41810BL', 'MSL7182BL', 'MSL4187BL')
n_sample <- length(velocity_samples)

## extract velocity from cellranger
for(i in 1:n_sample){
    # velocity of sample
    f <- velocity_samples[i]
    f_dir <- paste0(velocity_dirs, f, '/velocyto/')
    f_velo <- paste0(f_dir, list.files(f_dir))
    
    # loom
    f_loom <- read.loom.matrices(f_velo, engine = "hdf5r")
    
    # orig.ident
    orig_ident <- neutrophil_origident[i]
    orig_ident_spliced <- paste(orig_ident, 'spliced', sep = '_')
    orig_ident_unspliced <- paste(orig_ident, 'unspliced', sep = '_')
    orig_ident_ambiguous <- paste(orig_ident, 'ambiguous', sep = '_')
    
    # subset veloity matrix
    f_intersect_spliced <- VeloMatSubset(f_loom$spliced, neutrophil, orig_ident)
    f_intersect_unspliced <- VeloMatSubset(f_loom$unspliced, neutrophil, orig_ident)
    f_intersect_ambiguous <- VeloMatSubset(f_loom$ambiguous, neutrophil, orig_ident)
    
    assign(orig_ident_spliced, f_intersect_spliced)
    assign(orig_ident_unspliced, f_intersect_unspliced)
    assign(orig_ident_ambiguous, f_intersect_ambiguous)
    
    print(paste('The velocity of sample:',orig_ident, 'cell count is', ncol(f_intersect_spliced)))
}

## concatenent of all samples
paste0(neutrophil_origident, '_spliced')
velocity_spliced <- cbind(MNL3181BL_spliced, MNL4186BL_spliced, MML2182BL_spliced, MML41810BL_spliced, MSL7182BL_spliced, MSL4187BL_spliced)
velocity_unspliced <- cbind(MNL3181BL_unspliced, MNL4186BL_unspliced, MML2182BL_unspliced, MML41810BL_unspliced, MSL7182BL_unspliced, MSL4187BL_unspliced)
velocity_ambiguous <- cbind(MNL3181BL_ambiguous, MNL4186BL_ambiguous, MML2182BL_ambiguous, MML41810BL_ambiguous, MSL7182BL_ambiguous, MSL4187BL_ambiguous)

## add to Seurat obj
neutrophil[['spliced']] <- CreateAssayObject(counts = as.matrix(velocity_spliced))
neutrophil[['unspliced']] <- CreateAssayObject(counts = as.matrix(velocity_unspliced))
neutrophil[['ambiguous']] <- CreateAssayObject(counts = as.matrix(velocity_ambiguous))


save(neutrophil, file = 'RData/Neutrophil3_RNAVelocity.RData')