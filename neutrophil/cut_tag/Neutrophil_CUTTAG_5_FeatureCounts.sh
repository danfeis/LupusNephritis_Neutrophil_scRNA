#!/bin/sh
#PBS -e CUTTAG_FeatureCounts.e
#PBS -o CUTTAG_FeatureCounts.o
#PBS -N CUTTAG_FeatureCounts
#PBS -l nodes=tc6000:ppn=15,walltime=10000:00:00
#######################################################
## FeatureCounts for counting

cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
removeDuplicate_dir='removeDuplicate'
SAF_dir='SAF'
FeatureCounts_dir='FeatureCounts'

## 1. Spi1
samples1=("Spi1-W1" "Spi1-W2" "NL4-1S" "NL7-1S" "Spi1-L1" "Spi1-L2")
bam_files1=$(for sample in "${samples1[@]}"; do echo -n "${removeDuplicate_dir}/${sample}_bowtie2_sort_dupRemove.bam "; done) # Generate the list of BAM files
featureCounts -p --countReadPairs -F 'SAF' -T 4 \
    -a ${SAF_dir}/All_merged_peaks_Spi1.saf \
    -o ${FeatureCounts_dir}/FeatureCounts_Spi1.txt \
    $bam_files1


## 2. Stat1
samples2=("Stat-W1" "Stat-W2" "NL4-2T" "NL7-2T" "Stat-L1" "Stat-L2")
bam_files2=$(for sample in "${samples2[@]}"; do echo -n "${removeDuplicate_dir}/${sample}_bowtie2_sort_dupRemove.bam "; done) # Generate the list of BAM files
featureCounts -p --countReadPairs -F 'SAF' -T 4 \
    -a ${SAF_dir}/All_merged_peaks_Stat1.saf \
    -o ${FeatureCounts_dir}/FeatureCounts_Stat1.txt \
    $bam_files2

