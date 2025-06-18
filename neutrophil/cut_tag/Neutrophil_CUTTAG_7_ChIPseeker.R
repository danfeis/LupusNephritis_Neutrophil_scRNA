## ChIPseeker annotation
# conda activate R4.2

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(ggplot2)


###################################################################### 1. Annotation for peaks in HC, mild and severe
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
Peak_files <- list(Spi1_HC = 'BlacklistRemove/BlacklistRemove_Spi1_HC_peaks.bed',
                   Spi1_Mild = 'BlacklistRemove/BlacklistRemove_Spi1_Mild_peaks.bed',
                   Spi1_Severe = 'BlacklistRemove/BlacklistRemove_Spi1_Severe_peaks.bed',
                   Stat1_HC = 'BlacklistRemove/BlacklistRemove_Stat1_HC_peaks.bed',
                   Stat1_Mild = 'BlacklistRemove/BlacklistRemove_Stat1_Mild_peaks.bed',
                   Stat1_Severe = 'BlacklistRemove/BlacklistRemove_Stat1_Severe_peaks.bed')

## Annotation
peakAnno1 <- annotatePeak(Peak_files[[1]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Spi1_HC, 3836
peakAnno2 <- annotatePeak(Peak_files[[2]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Spi1_Mild, 9565
peakAnno3 <- annotatePeak(Peak_files[[3]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Spi1_Severe, 2770
peakAnno4 <- annotatePeak(Peak_files[[4]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Stat1_HC, 3731
peakAnno5 <- annotatePeak(Peak_files[[5]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Stat1_Mild, 19862
peakAnno6 <- annotatePeak(Peak_files[[6]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db") # Stat1_Severe, 6890

peakAnno_all <- rbind(data.frame(peakAnno1) %>% mutate(TF = 'Spi1', Condition = 'HC'),
                      data.frame(peakAnno2) %>% mutate(TF = 'Spi1', Condition = 'Mild'),
                      data.frame(peakAnno3) %>% mutate(TF = 'Spi1', Condition = 'Severe'),
                      data.frame(peakAnno4) %>% mutate(TF = 'Stat1', Condition = 'HC'),
                      data.frame(peakAnno5) %>% mutate(TF = 'Stat1', Condition = 'Mild'),
                      data.frame(peakAnno6) %>% mutate(TF = 'Stat1', Condition = 'Severe')) # 46654


write.table(data.frame(peakAnno1), 'ChIPseeker/ChIPseekerAnno_Spi1_HC.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(data.frame(peakAnno2), 'ChIPseeker/ChIPseekerAnno_Spi1_Mild.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(data.frame(peakAnno3), 'ChIPseeker/ChIPseekerAnno_Spi1_Severe.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(data.frame(peakAnno4), 'ChIPseeker/ChIPseekerAnno_Stat1_HC.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(data.frame(peakAnno5), 'ChIPseeker/ChIPseekerAnno_Stat1_Mild.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(data.frame(peakAnno6), 'ChIPseeker/ChIPseekerAnno_Stat1_Severe.txt', quote = FALSE, sep = '\t', row.names = FALSE)

write.table(peakAnno_all, 'ChIPseeker/ChIPseekerAnno_Peaks_All.txt', quote = FALSE, sep = '\t', row.names = FALSE)


###################################################################### 2. Annotation for Differential Peaks
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
DiffBindingRegions <- list(DESeq_Count_All = 'DESeq2/DESeq2_results_All.bed')

## Annotation
peakAnno7 <- annotatePeak(DiffBindingRegions[[1]], tssRegion=c(-500, 100), TxDb=txdb, annoDb="org.Mm.eg.db")

## Combind
DESeq2 <- read.table('DESeq2/DESeq2_results_All.tsv', header=1)[,-c(1:3)]
colnames(DESeq2) <- paste0('DESeq2_',colnames(DESeq2))
peakAnno_all2 <- cbind(data.frame(peakAnno7),DESeq2)

write.table(peakAnno_all2, 'ChIPseeker/ChIPseekerAnno_DiffBindingRegions_All.txt', quote = FALSE, sep = '\t', row.names = FALSE)

## Symbol for GO anlaysis
GO_Spi1_LNEst_up <- peakAnno_all2 %>% filter(DESeq2_TF=='Spi1',DESeq2_Condition=='LNEstablish',DESeq2_Regulate=='Upregulate') # 0
GO_Spi1_LNEst_down <- peakAnno_all2 %>% filter(DESeq2_TF=='Spi1',DESeq2_Condition=='LNEstablish',DESeq2_Regulate=='Downregulate') # 0
GO_Spi1_LNDev_up <- peakAnno_all2 %>% filter(DESeq2_TF=='Spi1',DESeq2_Condition=='LNProgression',DESeq2_Regulate=='Upregulate') # 1047
GO_Spi1_LNDev_down <- peakAnno_all2 %>% filter(DESeq2_TF=='Spi1',DESeq2_Condition=='LNProgression',DESeq2_Regulate=='Downregulate') # 20

GO_Stat1_LNEst_up <- peakAnno_all2 %>% filter(DESeq2_TF=='Stat1',DESeq2_Condition=='LNEstablish',DESeq2_Regulate=='Upregulate') # 0
GO_Stat1_LNEst_down <- peakAnno_all2 %>% filter(DESeq2_TF=='Stat1',DESeq2_Condition=='LNEstablish',DESeq2_Regulate=='Downregulate') # 73
GO_Stat1_LNDev_up <- peakAnno_all2 %>% filter(DESeq2_TF=='Stat1',DESeq2_Condition=='LNProgression',DESeq2_Regulate=='Upregulate') # 1720
GO_Stat1_LNDev_down <- peakAnno_all2 %>% filter(DESeq2_TF=='Stat1',DESeq2_Condition=='LNProgression',DESeq2_Regulate=='Downregulate') # 39

write.table(GO_Spi1_LNEst_up$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Spi1_LNEst_up.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Spi1_LNEst_down$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Spi1_LNEst_down.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Spi1_LNDev_up$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Spi1_LNDev_up.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Spi1_LNDev_down$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Spi1_LNDev_down.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(GO_Stat1_LNEst_up$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Stat1_LNEst_up.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Stat1_LNEst_down$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Stat1_LNEst_down.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Stat1_LNDev_up$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Stat1_LNDev_up.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(GO_Stat1_LNDev_down$SYMBOL, 'ChIPseeker/DiffBindingRegions_Symbol_Stat1_LNDev_down.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)