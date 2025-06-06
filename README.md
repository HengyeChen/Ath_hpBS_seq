# Ath_hpBS_seq
pipeline for analyzing hpBS-seq data

## Step1. Mapping and pairing
### Input files
Input files are regular fastq.gz files
### Software requirement
This step requires bismark, bedtools, trimmomatic, fastqc, and samtools
### Output files
The final output files of this step are named as *_1_hp.bam and *_2_hp.bam

## Step2. Call methylation status at CHH, CG, CWG, CWWG, and CSG contexts
### Input files
_hp.bam files from step1
### Software requirement
bedtools
### Output files
Counts of unme, hemiW, hemiC, and fully-methylated dyads at each position are saved in *.all.count files. These files are in bed format. The four column after #chr #pos1 #pos2 are the count of unme, hemiW, hemiC, and fully-methylated dyads, respectively.

## Step3. Call DMRs
### Input files
all.count files
### Software and package requirement
bedtools, Python 3.12, Scipy, and Numpy
### Output files
p-values are saved in .me.txt and .hemi.txt for full- and hemi-methylation in all 1kb bins. 
The number of DMRs are saved in sig_summary.txt. The indexes of 1kb bin, p-values, and methylation frequencies of DMRs are save in .sig files.


## Other information 
sample reference files are uploaded in bed.gz format. Some files are not uploaded due to the large sizes.







