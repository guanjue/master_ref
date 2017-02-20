script_bin1='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'
script_bin2='/Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq/pipseq_procap_scripts/bin/'
script_bin3='/Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/pipseq_procap_scripts/bin/'

cat SGD_features_ORF_Verified.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' > SGD_features_ORF_Verified_start.bed
#python pipseq_procap_scripts/bin/dist_sort_bed.py
sort -k1,1 -k2,2n Xu_2009_ORF-Ts_V64_sort.TSS.1bp.bed > sort_ref.bed
sort -k1,1 -k2,2n /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/pipseq_pk/sua7_pipseq_56422_56428_merged_sacCer3_T_filtered_samestrand_s10e60F1.bed > sua7_pipseq_merged_sacCer3_T_filtered_samestrand_s5e20F2_TWO_intersect_uniq_pkcenter_sort.bed
bedtools closest -a sort_ref.bed -b sua7_pipseq_merged_sacCer3_T_filtered_samestrand_s5e20F2_TWO_intersect_uniq_pkcenter_sort.bed -D a -t last > ref_pip_pkcenter_up1.txt
cat ref_pip_pkcenter_up1.txt | sort -k10,10n | awk -F '\t' -v OFS='\t' '{if ($8!="-1") print $7,$8,$9,$4,".",$6}' > ref_pip_pkcenter_up1.bed
cat ref_pip_pkcenter_up1.txt | sort -k10,10n | awk -F '\t' -v OFS='\t' '{if ($8!="-1") print $7,$8-500,$9+500,$4,".",$6}' > ref_pip_pkcenter_up1001.bed
cat ref_pip_pkcenter_up1.txt | sort -k10,10n | awk -F '\t' -v OFS='\t' '{if ($8!="-1") print $7,$8-50,$9+50,$4,".",$6}' > ref_pip_pkcenter_up101.bed
bedtools getfasta -fi sacCer3.fa -bed ref_pip_pkcenter_up101.bed -fo ref_pip_pkcenter_up101.fa -s 
cp ref_pip_pkcenter_up10*1.* pipseq_center


### sort by PIP-seq PRO-Cap distance
sort -u Xu_2009_ORF_Ts_V64_TSS_start_sort_1kb_procap_TSS_1kb.bed.1.bed > Xu_2009_ORF_Ts_V64_TSS_start_sort_1kb_procap_TSS_1kb_uniq.bed.1.bed
python pipseq_procap_scripts/bin/dist_sort_bed.py -t Xu_2009_ORF_Ts_V64_TSS_start_sort_1kb_procap_TSS_1kb_uniq.bed.1.bed -m 1 -o 3 -s ref_pip_pkcenter_up1.bed -n 1 -p 3 -a Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.bed -b ref_pip_pkcenter_up1_sort0.bed
#python pipseq_procap_scripts/bin/dist_sort_bed.py -t SGD_features_ORF_Verified_start.bed -m 1 -o 3 -s Xu_2009_ORF-Ts_V64_sort.TSS.1bp.bed -n 1 -p 3 -a SGD_features_ORF_Verified_start_sort.bed -b Xu_2009_ORF-Ts_V64_TSS_sort.bed
### sort SGD ORF AUG
python pipseq_procap_scripts/bin/vlookup.py -t SGD_features_ORF_Verified_start.bed -m 4 -s Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.bed -n 4 -o SGD_features_ORF_Verified_start_sort.txt
python pipseq_procap_scripts/bin/vlookup.py -t Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.bed -m 4 -s SGD_features_ORF_Verified_start_sort.txt -n 4 -o Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.txt
python pipseq_procap_scripts/bin/vlookup.py -t sort_ref.bed -m 4 -s SGD_features_ORF_Verified_start_sort.txt -n 4 -o sort_ref_sort.txt
python pipseq_procap_scripts/bin/vlookup.py -t ref_pip_pkcenter_up1_sort0.bed -m 4 -s SGD_features_ORF_Verified_start_sort.txt -n 4 -o ref_pip_pkcenter_up1_sort0.txt

### convert vlookup result to bed format
cat Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.bed
cat ref_pip_pkcenter_up1_sort0.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > ref_pip_pkcenter_up1_sort0.bed
cat SGD_features_ORF_Verified_start_sort.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > SGD_features_ORF_Verified_start_sort.bed
cat sort_ref_sort.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > sort_ref_sort.bed

### expand to 1001 bp window
cat SGD_features_ORF_Verified_start_sort.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > SGD_features_ORF_Verified_start_sort_1kb.bed
cat sort_ref_sort.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > sort_ref_sort_1kb.bed
cat ref_pip_pkcenter_up1_sort0.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > ref_pip_pkcenter_up1_sort.bed
cat Xu_2009_ORF_Ts_V64_TSS_start_procap_TSS_1_sort.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > ref_procap_pkcenter_up1_sort.bed


#mkdir pk_based_pksort
python $script_bin1'gene_group_split.py' -t SGD_features_ORF_Verified_start_sort_1kb.bed -a 3 -g $script_bin3'gene_group_list.txt' -b 0 -i F
python $script_bin1'gene_group_split.py' -t sort_ref_sort_1kb.bed -a 3 -g $script_bin3'gene_group_list.txt' -b 0 -i F
python $script_bin1'gene_group_split.py' -t ref_pip_pkcenter_up1_sort.bed -a 3 -g $script_bin3'gene_group_list.txt' -b 0 -i F
python $script_bin1'gene_group_split.py' -t ref_procap_pkcenter_up1_sort.bed -a 3 -g $script_bin3'gene_group_list.txt' -b 0 -i F

python $script_bin2'gene_group_split_notation.py' -t ref_pip_pkcenter_up1_sort.bed -a 3 -g $script_bin1'gene_group_list.txt' -b 0 -i F

echo sua7_hs0_Tfiltered_56422_56428_merged

paste ref_procap_pkcenter_up1_sort.bed.genegroup SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print ($2-$8); else print -($2-$8) }' > SGD_procap_dist.txt 
paste ref_pip_pkcenter_up1_sort.bed.genegroup ref_procap_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print -($2-$9); else print ($2-$9) }' > TSS_procap_dist.txt ### 1 tab after the Xu_2009_ORF-Ts_V64_TSS_sort_1kb.bed.genegroup file
paste ref_pip_pkcenter_up1_sort.bed.genegroup SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print -($2-$9); else print ($2-$9) }' > TSS_SGD_dist.txt
Rscript $script_bin2'hist.R' SGD_procap_dist.txt TSS_procap_dist.txt TSS_SGD_dist.txt






rm -r pipseq
mkdir pipseq
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
mv *_T*_combined.tabular pipseq
cd pipseq
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined.tabular -r 1 -c 5 -o sort_ref_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined green4 10
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular sort_ref_sort_1kb_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined green4 10
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular ref_pip_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined green4 10
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined_binned.tabular ref_procap_pkcenter_up1_sort_sua7_hs0_Tfiltered_56422_56428_merged_readc_combined green4 10

cd ..



rm -r rpb3_pipseq
mkdir rpb3_pipseq
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/rpb3_hs0_Tfiltered_56419.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 1 -t 3 -w 0 -h true -m false
mv *_T*_combined.tabular rpb3_pipseq
cd rpb3_pipseq
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined.tabular -r 1 -c 5 -o sort_ref_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined green4 10
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular sort_ref_sort_1kb_rpb3_hs0_Tfiltered_56419_readc_combined green4 10
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular ref_pip_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined green4 10
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined_binned.tabular ref_procap_pkcenter_up1_sort_rpb3_hs0_Tfiltered_56419_readc_combined green4 10

cd ..



rm -r procap
mkdir procap
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
mv *_anti.tabular procap
mv *_sense.tabular procap
cd procap
### sense strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_sense.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR3031844_1_readc_sense.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR3031844_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_sense.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_sense.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_sense_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR3031844_1_readc_sense_binned.tabular sort_ref_sort_1kb_SRR3031844_1_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.tabular ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.tabular ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned blue4 1
### anti strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_anti.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR3031844_1_readc_anti.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR3031844_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_anti.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_anti.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_anti_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR3031844_1_readc_anti_binned.tabular sort_ref_sort_1kb_SRR3031844_1_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.tabular ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.tabular ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned red4 1

### merge BR
composite -dissolve 50 -transparent-color white SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_sense_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1_readc_anti_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR3031844_1.png
composite -dissolve 50 -transparent-color white sort_ref_sort_1kb_SRR3031844_1_readc_sense_binned.png sort_ref_sort_1kb_SRR3031844_1_readc_anti_binned.png sort_ref_sort_1kb_SRR3031844_1.png
composite -dissolve 50 -transparent-color white ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.png ref_pip_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.png ref_pip_pkcenter_up1_sort_SRR3031844_1.png
composite -dissolve 50 -transparent-color white ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_sense_binned.png ref_procap_pkcenter_up1_sort_SRR3031844_1_readc_anti_binned.png ref_procap_pkcenter_up1_sort_SRR3031844_1.png
cd ..




rm -r proseq # 3' end
mkdir proseq
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/proseq_bam/SRR3031836_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
mv *_anti.tabular proseq
mv *_sense.tabular proseq
cd proseq
### sense strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_sense.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR3031836_1_read1_sense.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR3031836_1_read1_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_sense.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_sense.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_sense_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR3031836_1_read1_sense_binned.tabular sort_ref_sort_1kb_SRR3031836_1_read1_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.tabular ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.tabular ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned blue4 1
### anti strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_anti.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR3031836_1_read1_anti.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR3031836_1_read1_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_anti.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_anti.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_anti_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_anti_binned red4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR3031836_1_read1_anti_binned.tabular sort_ref_sort_1kb_SRR3031836_1_read1_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.tabular ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.tabular ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned red4 1

### merge BR
composite -dissolve 50 -transparent-color white SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_sense_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1_read1_anti_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR3031836_1.png
composite -dissolve 50 -transparent-color white sort_ref_sort_1kb_SRR3031836_1_read1_sense_binned.png sort_ref_sort_1kb_SRR3031836_1_read1_anti_binned.png sort_ref_sort_1kb_SRR3031836_1.png
composite -dissolve 50 -transparent-color white ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.png ref_pip_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.png ref_pip_pkcenter_up1_sort_SRR3031836_1.png
composite -dissolve 50 -transparent-color white ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_sense_binned.png ref_procap_pkcenter_up1_sort_SRR3031836_1_read1_anti_binned.png ref_procap_pkcenter_up1_sort_SRR3031836_1.png
cd ..




rm -r total_rna
mkdir total_rna
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/total_RNA/total_RNA_mapped_all_reads.7.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p false -a 0 -t 3 -w 0 -h true -m false
mv *_anti.tabular total_rna
mv *_sense.tabular total_rna
cd total_rna
### sense strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_sense.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_sense.tabular -r 1 -c 5 -o sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.tabular sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.tabular ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned blue4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.tabular ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned blue4 1
### anti strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_anti.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_anti.tabular -r 1 -c 5 -o sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.tabular sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.tabular ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned red4 1
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.tabular ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned red4 1

### merge BR
composite -dissolve 50 -transparent-color white SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.png SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.png SGD_features_ORF_Verified_start_sort_1kb_total_RNA_mapped_all_reads.png
composite -dissolve 50 -transparent-color white sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_sense_binned.png sort_ref_sort_1kb_total_RNA_mapped_all_reads_readc_anti_binned.png sort_ref_sort_1kb_total_RNA_mapped_all_reads.png
composite -dissolve 50 -transparent-color white ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.png ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.png ref_pip_pkcenter_up1_sort_total_RNA_mapped_all_reads.png
composite -dissolve 50 -transparent-color white ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_sense_binned.png ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads_readc_anti_binned.png ref_procap_pkcenter_up1_sort_total_RNA_mapped_all_reads.png
cd ..





rm -r mtif # require paired reads
mkdir mtif
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p true -a 0 -t 3 -w 0 -h true -m false

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/P_2013_full_RNAseq/SRR518891_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
mv *_anti.tabular mtif
mv *_sense.tabular mtif
cd mtif
### sense strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_sense.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR518891_1_readc_sense.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR518891_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR518891_1_readc_sense.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR518891_1_readc_sense.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_sense_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_sense_binned blue4 3
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR518891_1_readc_sense_binned.tabular sort_ref_sort_1kb_SRR518891_1_readc_sense_binned blue4 3
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.tabular ref_pip_pkcenter_up1_sort_SRR518891_1_readc_sense_binned blue4 3
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.tabular ref_procap_pkcenter_up1_sort_SRR518891_1_readc_sense_binned blue4 3
### anti strand
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_anti.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR518891_1_readc_anti.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR518891_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR518891_1_readc_anti.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR518891_1_readc_anti.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_anti_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_anti_binned red4 3
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR518891_1_readc_anti_binned.tabular sort_ref_sort_1kb_SRR518891_1_readc_anti_binned red4 3
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.tabular ref_pip_pkcenter_up1_sort_SRR518891_1_readc_anti_binned red4 3
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.tabular ref_procap_pkcenter_up1_sort_SRR518891_1_readc_anti_binned red4 3

### merge BR
composite -dissolve 50 -transparent-color white SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_sense_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1_readc_anti_binned.png SGD_features_ORF_Verified_start_sort_1kb_SRR518891_1.png
composite -dissolve 50 -transparent-color white sort_ref_sort_1kb_SRR518891_1_readc_sense_binned.png sort_ref_sort_1kb_SRR518891_1_readc_anti_binned.png sort_ref_sort_1kb_SRR518891_1.png
composite -dissolve 50 -transparent-color white ref_pip_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.png ref_pip_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.png ref_pip_pkcenter_up1_sort_SRR518891_1.png
composite -dissolve 50 -transparent-color white ref_procap_pkcenter_up1_sort_SRR518891_1_readc_sense_binned.png ref_procap_pkcenter_up1_sort_SRR518891_1_readc_anti_binned.png ref_procap_pkcenter_up1_sort_SRR518891_1.png
cd ..





rm -r mRNA # 3' reads
mkdir mRNA
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/mRNA_Booth/SRR3031848_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 0 -n 1 -e false -r 0 -p false -a 0 -t 3 -w 0 -h true -m false
mv *_anti.tabular mRNA
mv *_sense.tabular mRNA






rm -r chipexo
mkdir chipexo
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam.bai -c SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam.bai -c sort_ref_sort_1kb.bed.genegroup -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam.bai -c ref_pip_pkcenter_up1_sort.bed.genegroup -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/chipexo_bam/SRR1951311.uniq.bam.bai -c ref_procap_pkcenter_up1_sort.bed.genegroup -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
mv *_combined.tabular chipexo
cd chipexo
python $script_bin1'bin_row_col.py' -i SGD_features_ORF_Verified_start_sort_1kb_SRR1951311_read1_combined.tabular -r 1 -c 5 -o SGD_features_ORF_Verified_start_sort_1kb_SRR1951311_read1_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i sort_ref_sort_1kb_SRR1951311_read1_combined.tabular -r 1 -c 5 -o sort_ref_sort_1kb_SRR1951311_read1_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_pip_pkcenter_up1_sort_SRR1951311_read1_combined.tabular -r 1 -c 5 -o ref_pip_pkcenter_up1_sort_SRR1951311_read1_combined_binned.tabular
python $script_bin1'bin_row_col.py' -i ref_procap_pkcenter_up1_sort_SRR1951311_read1_combined.tabular -r 1 -c 5 -o ref_procap_pkcenter_up1_sort_SRR1951311_read1_combined_binned.tabular

Rscript $script_bin3'heatmap.R' SGD_features_ORF_Verified_start_sort_1kb_SRR1951311_read1_combined_binned.tabular SGD_features_ORF_Verified_start_sort_1kb_SRR1951311_read1_combined_binned green4 10
Rscript $script_bin3'heatmap.R' sort_ref_sort_1kb_SRR1951311_read1_combined_binned.tabular sort_ref_sort_1kb_SRR1951311_read1_combined_binned green4 10
Rscript $script_bin3'heatmap.R' ref_pip_pkcenter_up1_sort_SRR1951311_read1_combined_binned.tabular ref_pip_pkcenter_up1_sort_SRR1951311_read1_combined_binned green4 10
Rscript $script_bin3'heatmap.R' ref_procap_pkcenter_up1_sort_SRR1951311_read1_combined_binned.tabular ref_procap_pkcenter_up1_sort_SRR1951311_read1_combined_binned green4 10

cd ..




### 
#python $script_bin1'gene_group_split.py' -t SGD_features_ORF_Verified_start_sort.bed -a 3 -g $script_bin1'gene_group_list.txt' -b 0 -i F
#python $script_bin1'gene_group_split.py' -t Xu_2009_ORF-Ts_V64_TSS_sort.bed -a 3 -g $script_bin1'gene_group_list.txt' -b 0 -i F
#python $script_bin1'gene_group_split.py' -t Xu_2009_ORF_Ts_V64_TSS_start_sort_1kb_procap_TSS_1kb.bed.1001.bed -a 3 -g $script_bin1'gene_group_list.txt' -b 0 -i F

cat SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6 }' | uniq > SGD_features_ORF_Verified_start_sort_forfasta_1001.bed
cat sort_ref_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6 }' | uniq > sort_ref_sort_forfasta_1001.bed
cat ref_pip_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6 }' | uniq > ref_pip_pkcenter_up1_sort_forfasta_1001.bed
cat ref_procap_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-500,($2+$3-1)/2+501,$4,$5,$6 }' | uniq > ref_procap_pkcenter_up1_sort_forfasta_1001.bed

cat SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-250,($2+$3-1)/2+251,$4,$5,$6 }' | uniq > SGD_features_ORF_Verified_start_sort_forfasta_501.bed
cat sort_ref_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-250,($2+$3-1)/2+251,$4,$5,$6 }' | uniq > sort_ref_sort_forfasta_501.bed
cat ref_pip_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-250,($2+$3-1)/2+251,$4,$5,$6 }' | uniq > ref_pip_pkcenter_up1_sort_forfasta_501.bed
cat ref_procap_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-250,($2+$3-1)/2+251,$4,$5,$6 }' | uniq > ref_procap_pkcenter_up1_sort_forfasta_501.bed

cat SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-50,($2+$3-1)/2+51,$4,$5,$6 }' | uniq > SGD_features_ORF_Verified_start_sort_forfasta_101.bed
cat sort_ref_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-50,($2+$3-1)/2+51,$4,$5,$6 }' | uniq > sort_ref_sort_forfasta_101.bed
cat ref_pip_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-50,($2+$3-1)/2+51,$4,$5,$6 }' | uniq > ref_pip_pkcenter_up1_sort_forfasta_101.bed
cat ref_procap_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{ print $1,($2+$3-1)/2-50,($2+$3-1)/2+51,$4,$5,$6 }' | uniq > ref_procap_pkcenter_up1_sort_forfasta_101.bed

cat SGD_features_ORF_Verified_start_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,($2+$3-1)/2,($2+$3-1)/2+2,$4,$5,$6; else print $1,($2+$3-1)/2-1,($2+$3-1)/2+1,$4,$5,$6 }' | uniq > SGD_features_ORF_Verified_start_sort_forfasta_2.bed
cat sort_ref_sort_1kb.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,($2+$3-1)/2,($2+$3-1)/2+2,$4,$5,$6; else print $1,($2+$3-1)/2-1,($2+$3-1)/2+1,$4,$5,$6 }' | uniq > sort_ref_sort_forfasta_2.bed
cat ref_pip_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,($2+$3-1)/2,($2+$3-1)/2+2,$4,$5,$6; else print $1,($2+$3-1)/2-1,($2+$3-1)/2+1,$4,$5,$6 }' | uniq > ref_pip_pkcenter_up1_sort_forfasta_2.bed
cat ref_procap_pkcenter_up1_sort.bed.genegroup | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,($2+$3+1)/2,($2+$3+1)/2+2,$4,$5,$6; else print $1,($2+$3+1)/2-1,($2+$3+1)/2+1,$4,$5,$6 }' | uniq > ref_procap_pkcenter_up1_sort_forfasta_2.bed

###
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed SGD_features_ORF_Verified_start_sort_forfasta_101.bed -fo SGD_features_ORF_Verified_start_sort_forfasta_101.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed sort_ref_sort_forfasta_101.bed -fo sort_ref_sort_forfasta_101.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_pip_pkcenter_up1_sort_forfasta_101.bed -fo ref_pip_pkcenter_up1_sort_forfasta_101.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_procap_pkcenter_up1_sort_forfasta_101.bed -fo ref_procap_pkcenter_up1_sort_forfasta_101.fa -s

bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed SGD_features_ORF_Verified_start_sort_forfasta_1001.bed -fo SGD_features_ORF_Verified_start_sort_forfasta_1001.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed sort_ref_sort_forfasta_1001.bed -fo sort_ref_sort_forfasta_1001.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_pip_pkcenter_up1_sort_forfasta_1001.bed -fo ref_pip_pkcenter_up1_sort_forfasta_1001.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_procap_pkcenter_up1_sort_forfasta_1001.bed -fo ref_procap_pkcenter_up1_sort_forfasta_1001.fa -s 

bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed SGD_features_ORF_Verified_start_sort_forfasta_2.bed -fo SGD_features_ORF_Verified_start_sort_forfasta_2.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed sort_ref_sort_forfasta_2.bed -fo sort_ref_sort_forfasta_2.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_pip_pkcenter_up1_sort_forfasta_2.bed -fo ref_pip_pkcenter_up1_sort_forfasta_2.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed ref_procap_pkcenter_up1_sort_forfasta_2.bed -fo ref_procap_pkcenter_up1_sort_forfasta_2.fa -s 

Rscript $script_bin2'dimer_piechart.R' ref_procap_pkcenter_up1_sort_forfasta_2.fa ref_procap_pkcenter_up1_sort_forfasta_2.png

#rm -r fasta
mkdir fasta
mv *2.fa fasta
mv *01.fa fasta
mv *2.bed fasta
cp *01.bed fasta
mv *_2.png fasta
mv Xu_2009_ORF_Ts_V64_TSS_start_sort_1kb_procap_TSS_101.fa.* fasta
cd fasta
python $script_bin2'fasta2matrix.py' -i ref_pip_pkcenter_up1_sort_forfasta_101.fa
python $script_bin2'fasta2matrix.py' -i ref_procap_pkcenter_up1_sort_forfasta_101.fa
mkdir dnashape
cp SGD_features_ORF_Verified_start_sort_forfasta_10*1.fa dnashape
cp sort_ref_sort_forfasta_10*1.fa dnashape
cp ref_pip_pkcenter_up1_sort_forfasta_10*1.fa dnashape
cp ref_procap_pkcenter_up1_sort_forfasta_10*1.fa dnashape
cd dnashape
Rscript $script_bin3'plotDNAshape.R' SGD_features_ORF_Verified_start_sort_forfasta_1001.fa SGD_features_ORF_Verified_start_sort_forfasta_1001.pdf SGD_features_ORF_Verified_start_sort_forfasta_1001.png 502
Rscript $script_bin3'plotDNAshape.R' sort_ref_sort_forfasta_1001.fa sort_ref_sort_forfasta_1001.pdf sort_ref_sort_forfasta_1001.png 502
Rscript $script_bin3'plotDNAshape.R' ref_pip_pkcenter_up1_sort_forfasta_1001.fa ref_pip_pkcenter_up1_sort_forfasta_1001.pdf ref_pip_pkcenter_up1_sort_forfasta_1001.png 502
Rscript $script_bin3'plotDNAshape.R' ref_procap_pkcenter_up1_sort_forfasta_1001.fa ref_procap_pkcenter_up1_sort_forfasta_1001.pdf ref_procap_pkcenter_up1_sort_forfasta_1001.png 502

Rscript $script_bin3'plotDNAshape.R' SGD_features_ORF_Verified_start_sort_forfasta_101.fa SGD_features_ORF_Verified_start_sort_forfasta_101.pdf SGD_features_ORF_Verified_start_sort_forfasta_101.png 502
Rscript $script_bin3'plotDNAshape.R' sort_ref_sort_forfasta_101.fa sort_ref_sort_forfasta_101.pdf sort_ref_sort_forfasta_101.png 502
Rscript $script_bin3'plotDNAshape.R' ref_pip_pkcenter_up1_sort_forfasta_101.fa ref_pip_pkcenter_up1_sort_forfasta_101.pdf ref_pip_pkcenter_up1_sort_forfasta_101.png 502
Rscript $script_bin3'plotDNAshape.R' ref_procap_pkcenter_up1_sort_forfasta_101.fa ref_procap_pkcenter_up1_sort_forfasta_101.pdf ref_procap_pkcenter_up1_sort_forfasta_101.png 502

cd ..
cd ..

