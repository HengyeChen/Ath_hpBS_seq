#!/bin/bash

#######call p-value at 1kb resolution
fd=/home/chy/Ath/Ath/fq/ddm1/ddm1/bam
mkdir -p $fd/diff
cd $fd/diff
mkdir -p DMR
res='1kb'
region='/home/chy/ref/TAIR10_'${res}'.n.bed'
for mut in ddm1
do
    for motif in cg cwg cwwg
    do
        bedtools intersect -a /home/chy/Ath/Ath/bam/wt/Ath-wt_hpBS_comb_hp.${motif}.all.count -b $fd/${mut}/Ath-${mut}_hpBS_comb_hp.${motif}.all.count -wo | \
            awk '{a=($4+$5+$6+$7);b=($11+$12+$13+$14);print $1,$2,$3,($5+$6)/a,$7/a,($12+$13)/b,$14/b}' OFS='\t' > DMR/wt_${mut}_all_${motif}.perc_for_p.bed

        cat $region | awk '{print $0,FNR}' OFS='\t' | bedtools intersect -a - -b DMR/wt_${mut}_all_${motif}.perc_for_p.bed -wo | \
            awk '{print $1,$2,$3,$8,$9,$10,$11,$4}' OFS='\t' > DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.bed
        
        python3 /home/chy/software/call_pvalue_v5.py DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.bed DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.txt 4 6 8
        python3 /home/chy/software/call_pvalue_v5.py DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.bed DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.txt 5 7 8
    done
done