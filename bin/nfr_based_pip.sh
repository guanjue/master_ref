
cd /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/

python /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/get_NFR.py -f Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt -w 100 -s 25 -l 79 -u 1000 -d 0

paste Xu_2009_ORF-Ts_V64_sort.TSS.2001.bed Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txt > Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmp

sort -k10,10n Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmp > Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmpsort

cat Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmpsort | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $7!="na") print $1,$2+$7-500,$2+$7+501,$4,$5,$6; else if ($6=="-" && $7!="na") print $1,$3-$7-500,$3-$7+501,$4,$5,$6; else if ($7=="na") print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed

cat Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmpsort | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $8!="na") print $1,$2+$8-500,$2+$8+501,$4,$5,$6; else if ($6=="-" && $8!="na") print $1,$3-$8-500,$3-$8+501,$4,$5,$6; else if ($8=="na") print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined_nfr.bed

cat Xu_2009_ORF-Ts_V64_sort_h3_merged_MN_seq_sort_read1_combined.cdt.nucleosome.txttmpsort | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $9!="na") print $1,$2+$9-500,$2+$9+501,$4,$5,$6; else if ($6=="-" && $9!="na") print $1,$3-$9-500,$3-$9+501,$4,$5,$6; else if ($9=="na") print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed

cat Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed | awk '{print $4}' > nfr_based_genename_order_nfrbased.txt



rm -r compare_plus1_nfrbased
mkdir compare_plus1_nfrbased 
mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed compare_plus1_nfrbased
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed compare_plus1_nfrbased
mv nfr_based_genename_order_nfrbased.txt compare_plus1_nfrbased
cp gene_group_list.txt compare_plus1_nfrbased
cp ScriptManager-v0.10.jar compare_plus1_nfrbased

cd compare_plus1_nfrbased
mkdir nfr_based_nfrsort
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed -s 0 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed -s 0 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_56422_T_read1_combined.tabular nfr_based_nfrsort
#mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_56422_T_pip_seq.txt nfr_based_nfrsort
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_56422_T_read1_combined.tabular nfr_based_nfrsort
#mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_56422_T_pip_seq.txt nfr_based_nfrsort

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_SRR3031844_1_read1_*.tabular nfr_based_nfrsort
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_SRR3031844_1_read1_*.tabular nfr_based_nfrsort

cd nfr_based_nfrsort
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/heatmap.R Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_56422_T_read1_combined.tabular Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_56422_T_read1_combined.png
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/heatmap.R Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_56422_T_read1_combined.tabular Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_56422_T_read1_combined.png
cd ..

mkdir nfr_based_nfrsort_genegroup
python /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/gene_group_split.py -t Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed -a 3 -g gene_group_list.txt -b 0 -i F
python /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/gene_group_split.py -t Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed -a 3 -g gene_group_list.txt -b 0 -i F

mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort.bed.genegroup Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup.bed
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort.bed.genegroup Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup.bed

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup.bed -s 0 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/56422_T.uniq.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup.bed -s 0 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup_56422_T_read1_combined.tabular nfr_based_nfrsort_genegroup
#mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup_56422_T_pip_seq.txt nfr_based_nfrsort_genegroup
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_56422_T_read1_combined.tabular nfr_based_nfrsort_genegroup
#mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_56422_T_pip_seq.txt nfr_based_nfrsort_genegroup

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup.bed -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup.bed -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false -o /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/compare_plus1_nfrbased/
mv Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup_SRR3031844_1_read1_*.tabular nfr_based_nfrsort_genegroup
mv Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_SRR3031844_1_read1_*.tabular nfr_based_nfrsort_genegroup

cd nfr_based_nfrsort_genegroup
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/heatmap.R Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup_56422_T_read1_combined.tabular Xu_2009_ORF-Ts_V64_sort_nfr_minus1_nfrbased_sort_genegroup_56422_T_read1_combined.png
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/heatmap.R Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_56422_T_read1_combined.tabular Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_56422_T_read1_combined.png
cd ..


mkdir old_plus1
cd old_plus1
echo 'Yeast_plus_one_sacCer3_TSS_distsort.bed' > gene_order.txt
python /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/gene_match.py -t Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_h3_merged_MN_seq_sort_read1_combined.cdt -a 0 -g gene_order.txt -b 3 -i T
python /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/master_ref/bin/gene_sort.py -t Yeast_plus_one_sacCer3_TSS_distsort_h3_merged_MN_seq_sort_read1_combined.cdt -a 0 -r Xu_2009_ORF-Ts_V64_sort_nfr_plus1_nfrbased_sort_genegroup_h3_merged_MN_seq_sort_read1_combined.cdt -b 0 -o Yeast_plus_one_sacCer3_TSS_distsort_h3_merged_MN_seq_sort_read1_combined_nfrbased_sort.cdt
cd ..
