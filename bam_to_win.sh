#!/bin/sh
# usage: insert bam fils into BAM_PATH dir.  output is writen into RES_PATH.  
WIN_SIZE=10000
GNM_SIZE=3095693983
BAM_PATH=bam_files/
BAM_FILES=$BAM_PATH"*.bam" 
RES_PATH=bam_to_window/
		
#get number of bam files in reference path
J=$(ls $BAM_FILES | grep ".bam" | wc -l);
echo $J

for f in $BAM_FILES
do
# create new name for result file 
RES_NAME=$(echo $f | sed 's/.bam/.win/' | sed 's,'"$BAM_PATH"',,')
echo $RES_PATH$RES_NAME
# number of maped reads in file (ignore random reads and unmaped reads)
N_READS=$(samtools view $f | grep -v 'chrUn' | grep -v 'random' | wc -l)
echo nmber of reads: $N_READS
# length of read in file
R_LEN=$(samtools view $f | head -1 | awk '{print length($10)}')
echo read length: $R_LEN
# calculate RPGC (reads per genomic content) normalization factor. Normalizes reads to 1X depth coverage.

NORM_FAC=$(echo $N_READS*$R_LEN/$GNM_SIZE | bc -l) # was NORM_FAC=$(echo $N_READS*$R_LEN/$GNM_SIZE | bc -l) RPGC from github.com/fidelram/deepTools/wiki/Normalizations
echo norm factor: $NORM_FAC 
# read through bam file, sum up number of reads in every window of size WIN_SIZE, and write results into result directory
# TODO delete 'head -5000000'
samtools depth $f | grep -v 'chrUn' | grep -v 'random' | sed 's/chr//' | sed 's/X/23/' | sed 's/Y/24/' | sed 's/M/25/' | awk -v WIN_SIZE=$WIN_SIZE '{a[$1][int($2/WIN_SIZE)]+=$3}END{for (chr in a)\
	{for (location in a[chr]) print chr, location*WIN_SIZE, (location+1)*WIN_SIZE, a[chr][location]/($NORM_FAC*$R_LEN)}}' > $RES_PATH$RES_NAME
done


#sed 's/chrX/22/' | sed 's/chrY/23/sa' | sed 's/chrYM/24' | awk '{for (i=1;i++;i<23) {sed 's/chr$i/$i/'}}'

#samtools depth $f | awk -v WIN_SIZE=1000 -v chr=chr -v loc=1000 -v sum=0 -v i=1 '{if ($1 != (chr i)) {i+=1; loc = WIN_SIZE; sum = 0}}{while ($1==(chr i))if($2<loc)\
#{sum+=$3; break}else {print chr i, loc-WIN_SIZE, loc, sum; loc += WIN_SIZE; sum = 0}}' > $RES_PATH$RES_NAME



