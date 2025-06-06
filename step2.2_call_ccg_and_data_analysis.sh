#!/bin/bash

#call ccg
path='/home/chy/Ath/Ath/fq/ddm1/ddm1/bam'

cd $path
mkdir -p ccg
for x in ddm1
do
    prefix='Ath-'${x}'_hpBS_comb_hp'
    index1='/home/chy/ref/Ath/Ath_CCG_Wat.n.bed'
    index2='/home/chy/ref/Ath/Ath_CCG_Cri.n.bed'
    #the order of three numbers is C1-C2-C3, no matter the CCG is on Wat or Cri strand. 0 and 1 represent unmethylated and methylated C.
    bedtools intersect -a $x/$prefix.bed -b $index1 -wo -F 1 | awk '{print $7,$8,$9,substr($4,$8-$2+1,1),substr($4,$8-$2+2,1),substr($5,$8-$2+3,1),$6}' OFS='\t' | grep -v '\.' | sed -e s/z/0/g -e s/Z/1/g -e s/x/0/g -e s/X/1/g > $x/$prefix.ccg_w
    bedtools intersect -a $x/$prefix.bed -b $index2 -wo -F 1 | awk '{print $7,$8,$9,substr($5,$9-$2,1),substr($5,$9-$2-1,1),substr($4,$9-$2-2,1),$6}' OFS='\t' | grep -v '\.' | sed -e s/z/0/g -e s/Z/1/g -e s/x/0/g -e s/X/1/g > $x/$prefix.ccg_c

    cat $x/$prefix.ccg_w $x/$prefix.ccg_c | awk '$6==0 || $6==1{print $4,$5,$6}' OFS='\t' > ccg/Ath-${x}_hpCCG_both.matrix
    cat ccg/Ath-${x}_hpCCG_both.matrix | awk '{print $1$2$3}' OFS='\t' | awk '{
        if($1=="011"){
            bc+=1
        }
        else if($1=="010"){
            b+=1
        }
        else if($1=="001"){
            c+=1
        }
        else if($1=="100"){
            a+=1
        }
        else if($1=="101"){
            ac+=1
        }
        else if($1=="110"){
            ab+=1
        }
        else if($1=="111"){
            abc+=1
        }
        n++
        }END{print abc/n,ab/n,ac/n,a/n,c/n,b/n,bc/n}' OFS='\t' > ccg/Ath-${x}_hpCCG_both.stat.perc
done