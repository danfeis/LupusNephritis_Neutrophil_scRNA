#!/bin/sh
#PBS -e CUTTAG_Blacklist.e
#PBS -o CUTTAG_Blacklist.o
#PBS -N CUTTAG_Blacklist
#PBS -l nodes=tc6000:ppn=15,walltime=10000:00:00
#######################################################
## remove peak resions in blacklist

cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
Blacklist_dir='BlacklistRemove'
MACS2_dir='MACS2'

## remove blasklist
for narrowPeak in ${MACS2_dir}/*.narrowPeak; do
    base_name=$(basename "$narrowPeak" | awk -F'_' '{print $3"_"$4}')
    echo "Remove blacklist regions with $base_name"
    echo "$(wc -l $narrowPeak)"
    bedtools intersect -v \
        -a $narrowPeak \
        -b ${Blacklist_dir}/ENCFF547MET.bed | cut -f1-3 > ${Blacklist_dir}/BlacklistRemove_${base_name}_peaks.bed
done