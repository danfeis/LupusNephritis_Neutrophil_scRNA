## DESeq2 differential analysis

library(DESeq2)
library(tidyverse)
library(cowplot)

DESeqColData <- function(samples, TF, Condition1, Condition2){

    ### samples: sample names, usually from colnames of count dataframe
    ### TF: TFs in comparison
    ### Condition1: Condition1 in sample comparison
    ### Condition2: Condition1 in sample comparison
    ### return: Coldata for DESeq2 analysis

    coldata <- data.frame(SampleID = samples,
                        Factor = rep(TF, 4),
                        condition = c(Condition1, Condition1, Condition2, Condition2),
                        Replicate = c(1,2,1,2))
    return(coldata)
}

DESeqDiffAnalysis <- function(Featurecount_df, Allsamplelist, contrast){

    ### Featurecount_df: count dataframe with each sample a column and each peak region a row
    ### Allsamplelist: sample information dataframe with samples, TF and conditions
    ### contrast: Treatment VS Reference
    ### return: DEGList of DESeq2

    # coldata
    samples <- colnames(Featurecount_df)
    samples_info <- Allsamplelist[samples,] # extract corresponding sample information
    TF <- unique(samples_info$TF)
    condition1 <- unique(samples_info$Condition)[1] # condition1 for first 2 samples
    condition2 <- unique(samples_info$Condition)[2] # condition2 for last 2 samples
    coldata <- DESeqColData(samples, TF, condition1, condition2)

    # # external spike-in normalization factors
    # spike_in_factors <- NormFactor[samples,]$NormFactors

    # Diff Analysis
    dds <- DESeqDataSetFromMatrix(Featurecount_df, coldata, design=~condition) # create DESeq2 dataset
    # dds <- estimateSizeFactors(dds) # Estimate DESeq2 size factors
    # original_size_factors <- sizeFactors(dds)
    # adjusted_size_factors <- spike_in_factors * spike_in_factors # Adjust DESeq2 size factors with spike-in factors
    # sizeFactors(dds) <- adjusted_size_factors
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", contrast))
    out <- data.frame(res) %>% arrange(padj) %>% filter(!is.na(padj)) # remove NA
    print(dim(out))
    print(dim(out %>% filter(padj<0.05, abs(log2FoldChange)>1)))
    print(dim(out %>% filter(padj<0.05, log2FoldChange>1)))
    print(dim(out %>% filter(padj<0.05, log2FoldChange<(-1))))

    return(out)
    # return(list(DESeq2_result = dds,
    #             DESeq2_stat = out))
}


DESeqVolcanoPlot <- function(out){
    p <- ggplot() +
            geom_point(data = out, aes(x = log2FoldChange, y = -log10(pvalue))) + 
            geom_point(data = out %>% filter(padj<0.05, log2FoldChange > 1), aes(x = log2FoldChange, y = -log10(pvalue)), color = 'red') +
            geom_point(data = out %>% filter(padj<0.05, log2FoldChange < (-1)), aes(x = log2FoldChange, y = -log10(pvalue)), color = 'blue')
    return(p)
}



## Sample information, Spi1 count and Stat1 count
Allsampleinfo <- read.table('AllSampleInfo.txt')
spi1_count_df <- read.table('FeatureCounts/FeatureCounts_Spi1.txt', header=2)
stat1_count_df <- read.table('FeatureCounts/FeatureCounts_Stat1.txt', header=2)
rownames(spi1_count_df) <- spi1_count_df$Geneid
rownames(stat1_count_df) <- stat1_count_df$Geneid
colnames(spi1_count_df) <- gsub("removeDuplicate\\.", "", colnames(spi1_count_df)) %>%
  gsub("_bowtie2_sort_dupRemove\\.bam", "", .) %>%
  gsub("\\.", "-", .)
colnames(stat1_count_df) <- gsub("removeDuplicate\\.", "", colnames(stat1_count_df)) %>%
  gsub("_bowtie2_sort_dupRemove\\.bam", "", .) %>%
  gsub("\\.", "-", .)


## 1. DESeq2 for count
DESeq_Spi1_LNEstablish <- DESeqDiffAnalysis(spi1_count_df %>% select(c('Spi1-W1','Spi1-W2','NL4-1S','NL7-1S')), Allsampleinfo, contrast = c('Mild','WT')) # 3271, 31 sig with 0 up and 6 down
DESeq_Spi1_LNDevelop <- DESeqDiffAnalysis(spi1_count_df %>% select(c('NL4-1S','NL7-1S','Spi1-L1','Spi1-L2')), Allsampleinfo, contrast = c('Severe','Mild')) # 5680, 789 sig with 772 up and 17 down
DESeq_Stat1_LNEstablish <- DESeqDiffAnalysis(stat1_count_df %>% select(c('Stat-W1','Stat-W2','NL4-2T','NL7-2T')), Allsampleinfo, contrast = c('Mild','WT')) # 8186, 33 sig with 0 up and 33 down
DESeq_Stat1_LNDevelop <- DESeqDiffAnalysis(stat1_count_df %>% select(c('NL4-2T','NL7-2T','Stat-L1','Stat-L2')), Allsampleinfo, contrast = c('Severe','Mild')) # 9845, 1124 sig with 1071 up and 53 down

p1_1 <- DESeqVolcanoPlot(DESeq_Spi1_LNEstablish)
p1_2 <- DESeqVolcanoPlot(DESeq_Spi1_LNDevelop)
p1_3 <- DESeqVolcanoPlot(DESeq_Stat1_LNEstablish)
p1_4 <- DESeqVolcanoPlot(DESeq_Stat1_LNDevelop)
p1 <- plot_grid(p1_1,p1_2,p1_3,p1_4)

ggsave('DESeq2/DESeq_Volcano_Count.pdf',p1, width = 10, height = 10)

#$ write significant results
DESeq_Count_All <- rbind(DESeq_Spi1_LNEstablish %>% filter(padj<0.05, abs(log2FoldChange)>1) %>% arrange(desc(log2FoldChange)) %>% mutate(TF = 'Spi1', Condition = 'LNEstablish') %>% mutate(Regulate = ifelse(log2FoldChange>1, 'Upregulate', 'Downregulate')),
                         DESeq_Spi1_LNDevelop %>% filter(padj<0.05, abs(log2FoldChange)>1) %>% arrange(desc(log2FoldChange)) %>% mutate(TF = 'Spi1', Condition = 'LNProgression') %>% mutate(Regulate = ifelse(log2FoldChange>1, 'Upregulate', 'Downregulate')),
                         DESeq_Stat1_LNEstablish %>% filter(padj<0.05, abs(log2FoldChange)>1) %>% arrange(desc(log2FoldChange)) %>% mutate(TF = 'Stat1', Condition = 'LNEstablish') %>% mutate(Regulate = ifelse(log2FoldChange>1, 'Upregulate', 'Downregulate')),
                         DESeq_Stat1_LNDevelop %>% filter(padj<0.05, abs(log2FoldChange)>1) %>% arrange(desc(log2FoldChange)) %>% mutate(TF = 'Stat1', Condition = 'LNProgression') %>% mutate(Regulate = ifelse(log2FoldChange>1, 'Upregulate', 'Downregulate')))
DESeq_Count_All <- cbind(rownames(DESeq_Count_All) %>% data.frame(name = .) %>% separate(name, into = c('chr', 'start', 'end'), sep = "\\.") %>% mutate(start = as.numeric(start), end = as.numeric(end)), DESeq_Count_All)
write_tsv(DESeq_Count_All, 'DESeq2/DESeq2_results_All.tsv')
write.table(DESeq_Count_All[,1:3], 'DESeq2/DESeq2_results_All.bed', col.names = FALSE, sep = '\t', row.names = FALSE, quote = FALSE) # bed files
