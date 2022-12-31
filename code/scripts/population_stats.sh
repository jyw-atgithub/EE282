#!/bin/bash
conda activate everything
cd /home/jenyuw/DSPR_snp/snp_processing

basepath=/home/jenyuw/DSPR_snp/snp_processing
raw="/home/jenyuw/Reference_genome"
processed=$basepath/processed
results=/home/jenyuw/DSPR_snp/snp_processing/results
mergedSNP=/home/jenyuw/DSPR_snp/raw_fq/together/parr_out.annotated.vcf

#getting the Tajima's D in sliding windows
## with VCFtools. Many "nan" values were generated. non-overlaping windows
#vcftools --gzvcf dmel-all-r6.46.digSNPs.vcf.gz --out tajimasd --TajimaD 10000
## with vcf-kit. non-overlaping windows
#vk tajima 10000 10000 dmel-all-r6.46.digSNPs.vcf.gz > tajimasd.vk

for i in $processed/*.vcf.gz
do
#echo $i
tempname=$(basename ${i})
#echo $(basename ${i})
name=${tempname:0:-7}
#echo $name
vk tajima 10000 10000 $i > $results/$name.tajimasd.vk
done

vk tajima 10000 10000 A1.annotated.vcf.gz > A1.annotated.tajimasd.vk
vk tajima 100000 100000 A1.annotated.vcf.gz > A1.annotated.tajimasd.vk

vk tajima 10000 10000 DSPR.r6.SNPs.vcf.gz > DSPR.r6.tajimasd.vk
vk tajima 10000 10000 parr_out.annotated.vcf.gz > parr_out.tajimasd.vk


#vcftools --vcf - --freq  --out $processed/$dmel.${i}SNPs

#use --freq2 and bioawk to retain bi-allelic sites and remove allele freq

: <<'SKIP'
vcftools --gzvcf dmel-all-r6.46.synSNPs.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > dmel-all-r6.46.synSNPs.freq2
SKIP

for i in $processed/*.vcf.gz
do
#echo $i
intm=$(basename $i)
name=$(echo $intm | sed 's/dmel-all-r6.46.//; s/.vcf.gz//')
#echo "name is $name"
vcftools --gzvcf $i --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > $results/dmel-all-r6.46.$name.freq2 &
done

vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > DSPR.r6.SNPs.freq2

vcftools --gzvcf parr_out.annotated.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > parr_out.annotated.freq2

vcftools --gzvcf A1.annotated.vcf.gz --stdout --freq2 \
| bioawk -t -c header 'BEGIN {print "CHROM" "\t" "POS" "\t" "N_ALLELES" "\t" "N_CHR"} \
$3=="2" {print $1 "\t" $2 "\t" $3 "\t" $4}' > A1.annotated.freq2



echo "start to wait"
wait
echo "waiting ends"

for i in $processed/*.vcf.gz
do
intm=$(basename $i)
name=$(echo $intm| sed 's/dmel-all-r6.46.//; s/.vcf.gz//')
vcftools --gzvcf $i --stdout --window-pi 10000 > dmel-all-r6.46.$name.pi
done

vcftools --gzvcf parr_out.annotated.vcf.gz --stdout --window-pi 10000 > parr_out.annotated.pi
vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 10000 > A1.annotated.pi
#vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 100000 > A1.annotated.pi2
#vcftools --gzvcf A1.annotated.vcf.gz --stdout --window-pi 1000000 > A1.annotated.pi3
vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --window-pi 10000 > r6.annotated.pi
#vcftools --gzvcf DSPR.r6.SNPs.vcf.gz --stdout --window-pi 100000 > r6.annotated.pi2