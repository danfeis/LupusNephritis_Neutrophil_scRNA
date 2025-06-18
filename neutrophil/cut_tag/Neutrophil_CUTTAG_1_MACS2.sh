#!/bin/sh
#PBS -e CUTTAG_MACS2.e
#PBS -o CUTTAG_MACS2.o
#PBS -N CUTTAG_MACS2
#PBS -l nodes=tc6000:ppn=15,walltime=10000:00:00
#######################################################
## peak calling

source activate MACS2
cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
removeDuplicate_dir='removeDuplicate'
MACS2_dir='MACS2'

## Spi1
samples1=("Spi1-W1" "NL4-1S" "Spi1-L1")
samples2=("Spi1-W2" "NL7-1S" "Spi1-L2")
TFs=("Spi1" "Spi1" "Spi1")
Conditions=("HC" "Mild" "Severe")

for i in "${!samples1[@]}"; do
    sample1=${samples1[$i]}
    sample2=${samples2[$i]}
    TF=${TFs[$i]}
    condition=${Conditions[$i]}
    bam_file1=${removeDuplicate_dir}/${sample1}_bowtie2_sort_dupRemove.bam
    bam_file2=${removeDuplicate_dir}/${sample2}_bowtie2_sort_dupRemove.bam
    echo "Peak Calling with MACS2 for TF $TF in $condition with files $bam_file1 and $bam_file2"
    
    # Call MACS2 for paired samples
    macs2 callpeak -t $bam_file1 $bam_file2 \
        -B -g mm -f BAMPE \
        -n MACS2_Peak_${TF}_${condition} \
        --outdir ${MACS2_dir} \
        -q 0.05 --keep-dup all > ${MACS2_dir}/MACS2_Peak_${TF}_${condition}_summary.txt
done


## Stat1
samples1=("Stat-W1" "NL4-2T" "Stat-L1")
samples2=("Stat-W2" "NL7-2T" "Stat-L2")
TFs=("Stat1" "Stat1" "Stat1")
Conditions=("HC" "Mild" "Severe")

for i in "${!samples1[@]}"; do
    sample1=${samples1[$i]}
    sample2=${samples2[$i]}
    TF=${TFs[$i]}
    condition=${Conditions[$i]}
    bam_file1=${removeDuplicate_dir}/${sample1}_bowtie2_sort_dupRemove.bam
    bam_file2=${removeDuplicate_dir}/${sample2}_bowtie2_sort_dupRemove.bam
    echo "Peak Calling with MACS2 for TF $TF in $condition with files $bam_file1 and $bam_file2"
    
    # Call MACS2 for paired samples
    macs2 callpeak -t $bam_file1 $bam_file2 \
        -B -g mm -f BAMPE \
        -n MACS2_Peak_${TF}_${condition} \
        --outdir ${MACS2_dir} \
        -q 0.05 --keep-dup all > ${MACS2_dir}/MACS2_Peak_${TF}_${condition}_summary.txt
done