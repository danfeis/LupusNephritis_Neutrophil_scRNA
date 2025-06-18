#### Run R velocity of neutrophil

export OMP_NUM_THREADS=1

library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(gridGraphics)
library(ggplotify)
library(ggplot2)
library(RColorBrewer)

load('RData/Neutrophil3_RNAVelocity.RData')
source('/public/home/zhudf/scRNA_kidney_mouse/RNAVelocity/RNAVelocity_Script/RNAVelocity_Func_CombineLoom.R')

neutrophil <- RunVelocity(object = neutrophil, deltaT = 1, kCells = 20, fit.quantile = 0.02)
neutrophil <- RunVelocity(object = neutrophil, ncores = 30)

## plot with color of mature groups
Idents(neutrophil) <- 'mature_clusters'
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = neutrophil)))
names(x = ident.colors) <- levels(x = neutrophil)
my_colorpalette <- my_color_neutrophil <- c('#f8eea4','#ecd416','#a8d1d1','#4d9393','#d66b78','#982a37','#c1c3e4', '#7b7fc4', '#363a7a','#181a36', '#627a9d')
cell.colors <- my_colorpalette[Idents(neutrophil)]
names(x = cell.colors) <- colnames(x = neutrophil)

p <- show.velocity.on.embedding.cor(emb = Embeddings(object = neutrophil, reduction = "umap"), 
                               vel = Tool(object = neutrophil, slot = "RunVelocity"), 
                               n = 200, 
                               scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 1), 
                               cex = 1, 
                               arrow.scale = 4, 
                               show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, 
                               grid.n = 20, 
                               arrow.lwd = 1, 
                               do.par = FALSE, 
                               cell.border.alpha = 0.1,
                               n.cores = 30)


## velo to ggplot
GrabGrob <- function(){
  grid.echo()
  grid.grab()
}
p_velo <- GrabGrob()
p_velo <- as.ggplot(p_velo)

ggsave('ImagesSmallSize3/Neutrophil_Fig3_RNAVelocity.pdf', p_velo, width = 10, height = 10)
ggsave('ImagesSmallSize3/Neutrophil_Fig3_RNAVelocity.png', p_velo, width = 10, height = 10)
save(neutrophil, p, file = 'RData/Neutrophil3_RNAVelocity_RunVelo.RData')
