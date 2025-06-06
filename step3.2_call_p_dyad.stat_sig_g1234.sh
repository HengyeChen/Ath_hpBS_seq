fd='/home/chy/Ath/Ath/fq/ddm1/ddm1/bam/diff'
cd $fd
mkdir -p DMR
res='1kb'
region='/home/chy/ref/TAIR10_'${res}'.n.bed'
echo 1 | awk '{print "factor\ncg_hemi\ncg_me\ncwg_hemi\ncwg_me\ncwwg_hemi\ncwwg_me"}' > tmp1
mkdir -p DMR/join
mkdir -p DMR/sig
for mut in ddm1
do
	echo $mut > tmp2
    for motif in cg cwg cwwg
    do
        cat DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.txt | awk '$2<0.05 && $3<0.05 && ($4+$5)>0.1{print }' > DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.sig
        cat DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.txt | awk '$2<0.05 && $3<0.05 && ($4+$5)>0.1{print }' > DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.sig

	bedtools merge -i DMR/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.bed -c 4,5,6,7,8 -o mean | awk '{print $8,$0}' OFS='\t' | sort -k1,1 > tmps1
	sort -k1,1 DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.sig | join -t $'\t' -j 1 tmps1 - | awk '{print $2,$3,$4,$10,$11,$5,$7,$6,$8}' OFS='\t' | sort -k1,1 -k2,2n >  DMR/join/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.sig.join
	sort -k1,1 DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.sig | join -t $'\t' -j 1 tmps1 - | awk '{print $2,$3,$4,$10,$11,$6,$8,$5,$7}' OFS='\t' | sort -k1,1 -k2,2n >  DMR/join/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.sig.join

	cat DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.hemi.sig | wc -l >> tmp2
	cat DMR/sig/wt_${mut}_all_${motif}.perc_for_p.1kb_dyad.me.sig | wc -l >> tmp2
    done
	paste tmp1 tmp2 > tmp
	mv tmp tmp1
done
mv tmp1 DMR/sig/sig_summary.txt