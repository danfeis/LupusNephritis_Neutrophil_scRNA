## merge all peaks from HC, mild and severe together

cd /public/home/zhudf/scRNA_kidney_mouse/Share_Zhangjian3/CUTTAG3
Blacklist_dir='BlacklistRemove'
MergedPeaks='MergedPeaks'

## Spi1
bed_files=$(ls $Blacklist_dir/BlacklistRemove_Spi1_*)
cat ${bed_files[@]} | sort -k1,1 -k2,2n | bedtools merge > ${MergedPeaks}/All_merged_peaks_Spi1.bed

## Stat1
bed_files=$(ls $Blacklist_dir/BlacklistRemove_Stat1_*)
cat ${bed_files[@]} | sort -k1,1 -k2,2n | bedtools merge > ${MergedPeaks}/All_merged_peaks_Stat1.bed