## Annotate protein coding genes in CUT&Tag Chipseeker anno
## Annotation includes:
## 1. Peaks in mild, moderate and severe states
## 2. Peaks in DiffBinding Anno

library(rtracklayer)
library(tidyverse)

## extract pcg from mm10 gtf file
gtf_file <- '/public/home/Shenglab/DataShare/Reference/mouse/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf'
gtf <- import(gtf_file)
pcg <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"] # protein coding genes: 21700
pcg_symbol <- pcg %>% data.frame() %>% select(gene_id, gene_name)

## CUT&Tag ChIPseeker Anno
Peaks <- read_tsv('ChIPseekerAnno_Peaks_All.txt')
DiffRegions <- read_tsv('ChIPseekerAnno_DiffBindingRegions_All.txt')

## match with ENSEMBL, add column to tell ProteinCodingGenes or not
Peaks$ProteinCodingGenes <- ifelse(Peaks$ENSEMBL %in% pcg_symbol$gene_id, 'Yes', 'No')
DiffRegions$ProteinCodingGenes <- ifelse(DiffRegions$ENSEMBL %in% pcg_symbol$gene_id, 'Yes', 'No')

write.csv(Peaks, 'ChIPseekerAnno_Peaks_All_PCG.csv')
write.csv(DiffRegions, 'ChIPseekerAnno_DiffBindingRegions_All_PCG.csv')





########################################### summarize
## several peaks could be annotated to a same gene, so the PCG count could be less than peak count in summarise
Peaks %>% filter(ProteinCodingGenes=='Yes') %>% group_by(TF, Condition) %>% summarise(Count = n(), PCGCount = n_distinct(SYMBOL))
DiffRegions %>%  filter(ProteinCodingGenes=='Yes') %>% group_by(DESeq2_TF, DESeq2_Condition) %>% summarise(Count = n(), PCGCount = n_distinct(SYMBOL))

### 1.1 peaks annotated with PCGs
# TF    Condition Count PCGCount
#   <chr> <chr>     <int>    <int>
# 1 Spi1  HC         3158     2461
# 2 Spi1  Mild       8734     6797
# 3 Spi1  Severe     2479     2009
# 4 Stat1 HC         3038     2387
# 5 Stat1 Mild      17855    10645
# 6 Stat1 Severe     6304     4612

### 1.2 Differential Peaks annotated with PCGs
# DESeq2_TF DESeq2_Condition Count PCGCount
#   <chr>     <chr>            <int>    <int>
# 1 Spi1      LNProgression      948      799
# 2 Stat1     LNEstablish         63       63
# 3 Stat1     LNProgression     1559     1158


Peaks %>% filter(ProteinCodingGenes=='Yes')  %>% group_by(TF) %>% summarise(Count = n(), PCGCount = n_distinct(SYMBOL)) # PCG count of each TF and LN condition
DiffRegions %>% filter(ProteinCodingGenes=='Yes')  %>% group_by(DESeq2_TF) %>% summarise(Count = n(), PCGCount = n_distinct(SYMBOL))

#   TF    Count PCGCount
#   <chr> <int>    <int>
# 1 Spi1  14371     8524
# 2 Stat1 27197    11863
# TotalPCGs = 13052

# DESeq2_TF Count PCGCount
#   <chr>     <int>    <int>
# 1 Spi1        948      799
# 2 Stat1      1622     1216
# TotalPCGs = 1529


########################################### PCG overlapped with SCENIC
### SCENIC inferred target genes of SPI1 and STAT1
SCENIC <- read_tsv('../../Genelist/Genelist_SCENIC_Targetgene.txt')
SCENIC_Target <- SCENIC %>% filter(TF %in% c('Spi1','Stat1'), TargetGeneCoef>0.01) %>% select(TF, TargetGene) %>% unique() # 4191 in total, coef=0.01 consistent with previous table
SCENIC_Stat1_targets <- unique((SCENIC_target %>% filter(TF=='Stat1', TargetGeneCoef>0.01))$TargetGene) # 1614
SCENIC_Spi1_targets <- unique((SCENIC_target %>% filter(TF=='Spi1', TargetGeneCoef>0.01))$TargetGene) # 2577


### overlap
# Stat1, 1293 out of 11863 CUT&Tag confirmed in SCENIC (0.1089944%)
intersect((SCENIC_Target %>% filter(TF=='Stat1') %>% select(TargetGene) %>% unique())$TargetGene, # unique target gene in SCENIC
          (Peaks %>% filter(ProteinCodingGenes=='Yes', TF == 'Stat1') %>% select(SYMBOL) %>% unique())$SYMBOL) # unique SYMBOL in CUT&Tag
# Spi1, 1536 out of 8524 CUT&Tag confirmed in SCENIC (0.1801971%)
intersect((SCENIC_Target %>% filter(TF=='Spi1') %>% select(TargetGene) %>% unique())$TargetGene, # unique target gene in SCENIC
          (Peaks %>% filter(ProteinCodingGenes=='Yes', TF == 'Spi1') %>% select(SYMBOL) %>% unique())$SYMBOL) # unique SYMBOL in CUT&Tag



########################################### PCG overlapped with ISG
ISG <- read.table('../../Genelist/Genelist_IFN100.txt')$V1
## ISG identified by CUT&Tag
DiffRegions_ISG <- DiffRegions_ISG %>% select(SYMBOL,DESeq2_TF,DESeq2_Condition) %>% unique() %>% data.frame()