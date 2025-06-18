#PBS -e CUTTAG_Bigwig.e
#PBS -o CUTTAG_Bigwig.o
#PBS -N CUTTAG_Bigwig
#PBS -l nodes=tc6000:ppn=15,walltime=10000:00:00
#######################################################
## generate bigwig files from MACS2 produced bedgraph files

source activate MACS2
cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
MACS2_dir='MACS2'
Bigwig_dir='Bigwig'
genome_file='Genome_mm10/genome.fa.fai'

for bedgraph_file in ${MACS2_dir}/*_treat_pileup.bdg; do
    base_name=$(echo $(basename $bedgraph_file) | cut -d'_' -f1-4) # Get the name
    bigwig_file=${Bigwig_dir}/${base_name}.bw
    bedGraphToBigWig $bedgraph_file $genome_file $bigwig_file # Convert bedGraph to BigWig
    echo "Converted $bedgraph_file to $bigwig_file"
done