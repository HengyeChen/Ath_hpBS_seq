#!/bin/bash

#read processing:
tp='ddm1'
mkdir -p fastq
for i in Ath-${tp}_hpBS_R2
do
    trimmomatic PE -threads 48 $i'_'1.fq.gz $i'_'2.fq.gz $i'_'1_trim.fq s1 $i'_'2_trim.fq s2 ILLUMINACLIP:/home/chy/optionData/adapters/TruSeq3-PE-hp.fa:0:0:2 TRAILING:20 MINLEN:15
    fastqc -t 2 --nogroup $i'_'*_trim.fq
    mv $i'_'*_trim_fastqc.html fastq
done

#mapping:
for i in Ath-${tp}_hpBS_R2
do
    mkdir -p raw
    mkdir -p report
    bismark --multicore 15 -D 20 -R 3 --score_min L,0,-0.4 /ssd/index/bismark/Ath $i'_'1_trim.fq
    bismark --pbat --multicore 15 -D 20 -R 3 --score_min L,0,-0.4 /ssd/index/bismark/Ath $i'_'2_trim.fq
    samtools sort -m 1G -@ 48 $i'_'1_trim_bismark_bt2.bam > raw/$i'_'1.bam
    samtools sort -m 1G -@ 48 $i'_'2_trim_bismark_bt2.bam > raw/$i'_'2.bam
    rm $i'_'*_trim_bismark_bt2.bam $i'_'*_trim.fq
    mv $i'_'1_trim_bismark_bt2_SE_report.txt report/$i'_'1_bismark_report.txt
    mv $i'_'2_trim_bismark_bt2_SE_report.txt report/$i'_'2_bismark_report.txt
done

#calculate conversion rate:
mkdir -p report/lambda
for i in Ath-${tp}_hpBS_R2
do
    for j in 1 2
    do
        samtools view -H raw/$i'_'$j'.'bam > HD
        samtools view raw/$i'_'$j'.'bam | awk '$3=="chrL"' | cat HD - | samtools sort -@ 48 -m 1G - > $i'_'$j'_'lamb.bam
        bismark_methylation_extractor -s --parallel 48 --mbias_only -o report/lambda $i'_'$j'_'lamb.bam
        rm HD $i'_'$j'_'lamb.bam report/lambda/*png report/lambda/*bias.txt
    done
done

# select paired reads from the same
mkdir -p bam
for i in Ath-${tp}_hpBS_R2
do
    # Convert BAM to BED format for read 1, filter out organellar chromosomes, extract read name and coordinates, remove read suffix, sort by read name
    bamToBed -i raw/$i'_'1.bam | awk '{if($1!~/chr[CLMT]/) print $4,NR,$1,$2,$3,$6}' OFS='\t' | sed 's/_1:\w:0:\w*//g' | sort -k1,1 -S 9G > R1
    # Same as above for read 2, then join with read 1 data and identify hairpin reads based on matching coordinates and opposite strands
    bamToBed -i raw/$i'_'2.bam | awk '{if($1!~/chr[CLMT]/) print $4,NR,$1,$2,$3,$6}' OFS='\t' | sed 's/_2:\w:0:\w*//g' | sort -k1,1 -S 9G | join -j 1 R1 - | awk '{if($3==$8 && (($6=="+" && $4==$9) || ($6=="-" && $5==$10)) && $6!=$11){print $2 > "LN1"; print $7 > "LN2"}}'
    # Extract header from read 1 BAM file
    samtools view -H raw/$i'_'1.bam > HD1
    # Extract header from read 2 BAM file
    samtools view -H raw/$i'_'2.bam > HD2
    # Filter read 1 BAM file to keep only hairpin reads, add header, and sort
    samtools view raw/$i'_'1.bam | awk 'FNR==NR {h[$1];next} (FNR in h)' LN1 - | cat HD1 - | samtools sort -@ 48 -m 1G - > $i'_'1_hap.bam
    # Filter read 2 BAM file to keep only hairpin reads, add header, and sort
    samtools view raw/$i'_'2.bam | awk 'FNR==NR {h[$1];next} (FNR in h)' LN2 - | cat HD2 - | samtools sort -@ 48 -m 1G - > $i'_'2_hap.bam
    # Remove PCR duplicates from read 1 hairpin BAM file
    deduplicate_bismark -s --bam $i'_'1_hap.bam
    # Remove PCR duplicates from read 2 hairpin BAM file
    deduplicate_bismark -s --bam $i'_'2_hap.bam
    # Clean up intermediate files
    rm R1 LN* HD* $i'_'*_hap.bam
    # Rename deduplicated BAM files
    rename 'hap.deduplicated.bam' 'hp.bam' *
	mv *hp.bam bam/
    # Move deduplication reports to report directory
    mv $i'_'1_hap.deduplication_report.txt report/$i'_'1_dedup_report.txt
    mv $i'_'2_hap.deduplication_report.txt report/$i'_'2_dedup_report.txt
done