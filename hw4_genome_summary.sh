#!/bin/bash

conda activate ee282
cd /pub/jenyuw/classrepos/EE282
# saving the sequences
faFilter -maxSize=100000 dmel-all-chromosome-r6.48.fasta le_100k.fasta
faFilter -minSize=100001 dmel-all-chromosome-r6.48.fasta gt_100k.fasta
# checking the results
#bioawk -c fastx '{print $seq}' le_100k.fasta |awk '{print length}'|sort -r -n|head -100|less

faSize le_100k.fasta
faSize gt_100k.fasta
# a rapid way
#faFilter -minSize=100001 <(zcat dmel.chr.gz) /dev/stdout | faSize /dev/stdin

faSize -detailed le_100k.fasta >le_100k_length.txt
faSize -detailed gt_100k.fasta >gt_100k_length.txt
#another way
#bioawk -c fastx ' {print $name "\t" length($seq) "\t" gc($seq)}' le_100k.fasta >le_100k_length.txt

bioawk -c fastx ' {print $name "\t" gc($seq)}' le_100k.fasta >le_100k_gc.txt
bioawk -c fastx ' {print $name "\t" gc($seq)}' gt_100k.fasta >gt_100k_gc.txt

./plotting_hw4.r


sort -k 2,2 -n -r  le_100k_length.txt |gawk 'BEGIN {print "Name \t Length \t Assembly"} {print $1 "\t" $2 "\t" "le_100k"}' \
|plotCDF2 /dev/stdin le.100k.cdf.png

sort -k 2,2 -n -r  gt_100k_length.txt |gawk 'BEGIN {print "Name \t Length \t Assembly"} {print $1 "\t" $2 "\t" "gt_100k"}' \
|plotCDF2 /dev/stdin gt.100k.cdf.png