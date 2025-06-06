#!/bin/bash

######set parameters#########
fd='/home/chy/Ath/Ath/fq/ddm1/ddm1/bam'
cd $fd
threads=20
mem=20G
index_CG='/home/chy/ref/Ath/Ath_CpG.bed'
index_CWG='/home/chy/ref/Ath/Ath_CWG.bed'
index_CWWG='/home/chy/ref/Ath/Ath_CWWG.fa.n.bed'
index_CHH_W='/home/chy/ref/Ath/Ath_CHH_Wat.n.bed'
index_CHH_C='/home/chy/ref/Ath/Ath_CHH_Cri.n.bed'
#####list all bam files
list=$(ls *.bam)
#####
for i in $list
do
    ######Convert all bam files to bed files(.cv)#########
	gene=$(echo $i | sed -e s/-/\\t/g -e s/_/\\t/g | awk '{print $2}')
	if [ ! -d $gene ]; then
	mkdir $gene
	fi
	cv=$i".cv"
	samtools view $i -@ ${threads} | awk '{print $1,$3,$4,$6,$14,$2}' OFS='\t' | sed -e s/M//g -e s/X:Z://g -e s/_1//g -e s/_2//g > $cv

    ###### pair reads in sorted bam files #########
    prefix='Ath-'${gene}'_hpBS_comb_hp'
    x=$gene
    cat Ath-${x}_hpBS_R*_1_hp.bam.cv > $x/Ath-${x}_hpBS_1_hp.bam.tmp
    cat Ath-${x}_hpBS_R*_2_hp.bam.cv > $x/Ath-${x}_hpBS_2_hp.bam.tmp
    paste $x/Ath-${x}_hpBS_1_hp.bam.tmp $x/Ath-${x}_hpBS_2_hp.bam.tmp | awk '{if($1==$7)print $2,$3,$4,$5,$6,$8,$9,$10,$11,$12}' > $x/Ath-${x}_hpBS_comb_hp.match.tmp
    paste $x/Ath-${x}_hpBS_1_hp.bam.tmp $x/Ath-${x}_hpBS_2_hp.bam.tmp | awk '{if($1!=$7)print $0}' > $x/Ath-${x}_hpBS_comb_hp.unmatch.tmp
    rm $x/*bam.tmp
    wc -l $x/Ath-${x}_hpBS_comb_hp.match.tmp
    wc -l $x/Ath-${x}_hpBS_comb_hp.unmatch.tmp

    prefix='Ath-'${x}'_hpBS_comb_hp'
    file=$x/$prefix.unmatch.tmp
    file1=$file.1
    file2=$file.2
    awk '{print $1,$2,$3,$4,$5,$6}' OFS='\t' $file > $file1
    awk '{print $7,$8,$9,$10,$11,$12}' OFS='\t' $file > $file2
    awk '{print $1}' $file1 > $file1.list
    awk '{print $1}' $file2 > $file2.list
    grep -f $file1.list $file2 > $file2.tmp
    grep -f $file2.list $file1 > $file1.tmp
    sort -S ${mem} --parallel=${threads} -k1,1 $file1.tmp > $file1.sort
    sort -S ${mem} --parallel=${threads} -k1,1 $file2.tmp > $file2.sort
    paste $file1.sort $file2.sort | awk '{if($1==$7)print $2,$3,$4,$5,$6,$8,$9,$10,$11,$12}' >> $x/$prefix.match.tmp
    rm $x/*.sort
    rm $x/*.list
    sort -S ${mem} --parallel=${threads} -k1,1 $x/$prefix.match.tmp > $x/$prefix.match.sorted
    rm $x/*.tmp

    ######call CG, CWG, CWWG, and CHH#########
    awk '{if(length($3)==2 && length($8)==2)print $1,$2-1,$2+$3-1,$4,$5,$6,$7-1,$7+$8-1,$9,$10}' OFS='\t' $x/$prefix.match.sorted | awk '{
        if($2==$7 && $3==$8)
        {
            f=$4
            r=$9
            st=$2
            ed=$3
        }
        else if($2>=$7)
        {
            st=$2
            if($3==$8)
            {
                f=$4
                r=substr($9,$2-$7+1)
                ed=$3
            }
            else if($3>$8)
            {
                f=substr($4,1,length($4)-$3+$8)
                r=substr($9,$2-$7+1)
                ed=$8
            }
            else if($3<$8)
            {
                f=$4
                r=substr($9,$2-$7+1,length($4))
                ed=$3
            }
        }
        else if($2<$7)
        {
            st=$7
            if($3==$8)
            {
                f=substr($4,$7-$2+1)
                r=$9
                ed=$3
            }
            else if($3>$8)
            {
                f=substr($4,$7-$2+1,length($9))
                r=$9
                ed=$8
            }
            else if($3<$8)
            {
                f=substr($4,$7-$2+1)
                r=substr($9,1,length($9)+$2-$7)
                ed=$3
            }
        }
        if($5==0)
        {
            print $1,st,ed,f,r,FNR
        }
        else if($5==16)
        {
            print $1,st,ed,r,f,FNR
        }
    }' OFS='\t' | sed s/chr//g > $x/$prefix.bed

    #call CG##################
    motif='cg'
    bedtools intersect -a $x/$prefix.bed -b $index_CG -wo -F 1 | awk '{print $7,$8,$9,substr($4,$8-$2+1,1),substr($5,$9-$2,1),$6}' OFS='\t' | sed -e s/z/0/g -e s/Z/1/g > $x/$prefix.cg
    awk '{if($4==1 && $5==1)print $0}' $x/$prefix.cg > $x/$prefix.cg.me
    awk '{if($4==1 && $5==0)print $0}' $x/$prefix.cg > $x/$prefix.cg.hemiW
    awk '{if($4==0 && $5==1)print $0}' $x/$prefix.cg > $x/$prefix.cg.hemiC

    cat $fd/$x/$prefix.$motif | awk '{if($4==0 && $5==0){print $1,$2,$3,1,0,0,0,$6}else if($4==1 && $5==0){print $1,$2,$3,0,1,0,0,$6}else if($4==0 && $5==1){print $1,$2,$3,0,0,1,0,$6}else if($4==1 && $5==1){print $1,$2,$3,0,0,0,1,$6}}' OFS='\t' > $fd/$x/$prefix.$motif.stat
    cat $fd/$x/$prefix.$motif.stat | sort -S ${mem} --parallel=${threads} -k1,1 -k2,2n | bedtools merge -c 4,5,6,7 -o sum,sum,sum,sum -d -2 -i -> $x/$prefix.${motif}.all.count

    #call CWG##################
    motif='cwg'
    bedtools intersect -a $x/$prefix.bed -b $index_CWG -wo -F 1 | awk '{print $7,$8,$9,substr($4,$8-$2+1,1),substr($5,$9-$2,1),$6}' OFS='\t' | sed -e s/x/0/g -e s/X/1/g > $x/$prefix.cwg
    awk '{if($4==1 && $5==1)print $0}' $x/$prefix.cwg > $x/$prefix.cwg.me
    awk '{if($4==1 && $5==0)print $0}' $x/$prefix.cwg > $x/$prefix.cwg.hemiW
    awk '{if($4==0 && $5==1)print $0}' $x/$prefix.cwg > $x/$prefix.cwg.hemiC
    awk '{if($4==0 && $5==1)print $0}' $x/$prefix.cwg > $x/$prefix.cwg.unme

    cat $fd/$x/$prefix.$motif | awk '{if($4==0 && $5==0){print $1,$2,$3,1,0,0,0,$6}else if($4==1 && $5==0){print $1,$2,$3,0,1,0,0,$6}else if($4==0 && $5==1){print $1,$2,$3,0,0,1,0,$6}else if($4==1 && $5==1){print $1,$2,$3,0,0,0,1,$6}}' OFS='\t' > $fd/$x/$prefix.$motif.stat
    cat $fd/$x/$prefix.$motif.stat | sort -S ${mem} --parallel=${threads} -k1,1 -k2,2n | bedtools merge -c 4,5,6,7 -o sum,sum,sum,sum -d -2 -i -> $x/$prefix.${motif}.all.count

    #call CWWG##################
    motif='cwwg'
    bedtools intersect -a ${fd}/${x}/$prefix.bed -b $index_CWWG -wo -F 1 | awk '{print $7,$8,$9,substr($4,$8-$2+1,1),substr($5,$9-$2,1),$6,$10}' OFS='\t' | sed -e s/h/0/g -e s/H/1/g > ${fd}/${x}/$prefix.$motif
    awk '{if($4==1 && $5==1)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.me
    awk '{if($4==1 && $5==0)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.hemiW
    awk '{if($4==0 && $5==1)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.hemiC
    awk '{if($4==0 && $5==1)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.unme

    cat $fd/$x/$prefix.$motif | awk '{if($4==0 && $5==0){print $1,$2,$3,1,0,0,0,$6,$7}else if($4==1 && $5==0){print $1,$2,$3,0,1,0,0,$6,$7}else if($4==0 && $5==1){print $1,$2,$3,0,0,1,0,$6,$7}else if($4==1 && $5==1){print $1,$2,$3,0,0,0,1,$6,$7}}' OFS='\t' > $fd/$x/$prefix.$motif.stat
    cat $fd/$x/$prefix.$motif.stat | sort -S ${mem} --parallel=${threads} -k1,1 -k2,2n | bedtools merge -c 4,5,6,7 -o sum,sum,sum,sum -d -2 -i -> $x/$prefix.${motif}.all.count

    #call CHH
    motif='chh'
    bedtools intersect -a $fd/${x}/$prefix.bed -b $index_CHH_W -wo -F 1 | awk '{print $7,$8,$9,substr($4,$8-$2+1,1),substr($5,$9-$2,1),$6,$10}' OFS='\t' | sed -e s/h/0/g -e s/H/1/g > $fd/${x}/$prefix.chh
    bedtools intersect -a $fd/${x}/$prefix.bed -b $index_CHH_C -wo -F 1 | awk '{print $7,$8,$9,substr($4,$9-$2,1),substr($5,$8-$2+1,1),$6,$10}' OFS='\t' | sed -e s/h/0/g -e s/H/1/g >> $fd/${x}/$prefix.chh
    awk '{if($4==1 || $5==1)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.me
    awk '{if($4==0 || $5==0)print $0}' $fd/$x/$prefix.$motif > $fd/$x/$prefix.$motif.unme

    cat $fd/$x/$prefix.$motif | awk '{if($4==0 || $5==0){print $1,$2,$3,1,0,$6,$7}else if($4==1 || $5==1){print $1,$2,$3,0,1,$6,$7}}' OFS='\t' > $fd/$x/$prefix.$motif.stat
    cat $fd/$x/$prefix.$motif.stat | sort -S ${mem} --parallel=${threads} -k1,1 -k2,2n | bedtools merge -c 4,5 -o sum,sum -d -1 -i -> $x/$prefix.${motif}.all.count
    rm $x/*tmp*
	rm $x/*.me $x/*.hemiW $x/*.hemiC $x/*.unme
done
