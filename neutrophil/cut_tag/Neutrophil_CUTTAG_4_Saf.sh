#!/bin/sh
#PBS -e CUTTAG_SAF.e
#PBS -o CUTTAG_SAF.o
#PBS -N CUTTAG_SAF
#PBS -l nodes=tc6000:ppn=15,walltime=10000:00:00
#######################################################
## narrowpeak file to SAF

cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
MergedPeaks='MergedPeaks'
SAF_dir='SAF'

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${MergedPeaks}/All_merged_peaks_Spi1.bed > ${SAF_dir}/All_merged_peaks_Spi1.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${MergedPeaks}/All_merged_peaks_Stat1.bed > ${SAF_dir}/All_merged_peaks_Stat1.saf