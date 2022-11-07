# Homework 3

## First Part--Summarize Genome Assembly

Downloading the genome

    wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz

Checking the file integrity

    wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt
    md5sum -c md5sum.txt | grep dmel-all-chromosome-r6.48.fasta.gz |less

Using ucsc-faSize

    faSize dmel-all-chromosome-r6.48.fasta.gz |less

Here is the output:

    143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
    Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
    N count: mean 616.6 sd 6960.7
    U count: mean 76242.3 sd 1379508.4
    L count: mean 0.0 sd 0.0
    %0.00 masked total, %0.00 masked real


### **Answers to the questions:**

1. Total number of nucleotides is *143726002* bases.
2. Total number of Ns is *1152978*.
3. Total number of sequences is *1870*

-------

## Second Part--Summarize an Annotation File

Downloading the annotation file

    wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz


Checking the file integrity

    curl http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt -0 gtfmd5sum
    md5sum -c gtfmd5sum

Counting the number of each feature.

This is the first way which I figured out. Here, we use `cut` to isolate the feature column and sort them. After sorting, we can use `uniq` to count each entries. In order to present the result by numerical order `sort -h` or `-n` is required.

    gunzip -k dmel-all-r6.48.gtf.gz
    cut -f 3 dmel-all-r6.48.gtf| sort |uniq -c |sort -k 1 -h|less

The output:

    32 snRNA
    115 rRNA
    262 pre_miRNA
    300 snoRNA
    312 tRNA
    365 pseudogene
    485 miRNA
    3053 ncRNA
    17896 gene
    30799 mRNA
    30825 stop_codon
    30885 start_codon
    33738 3UTR
    46802 5UTR
    163242 CDS
    190050 exon

Counting the number of each feature with `bioawk`

To match the goal of this homework, `bioawk` is used to isolate the feature column. The following process are similar then previous. 

    bioawk -c gff '{ print $3 }' <dmel-all-r6.48.gtf.gz|sort|uniq -c|sort -n|less

The output is identical:

        32 snRNA
       115 rRNA
       262 pre_miRNA
       300 snoRNA
       312 tRNA
       365 pseudogene
       485 miRNA
      3053 ncRNA
     17896 gene
     30799 mRNA
     30825 stop_codon
     30885 start_codon
     33738 3UTR
     46802 5UTR
    163242 CDS
    190050 exon

Counting the number of genes on chromosomes.

First, we have to screen out the rows containing gene as the feature so the condition is `$3=="gene"`. Then, we print the column of sequence name as `{ print $1 }`. However, we only need the results from chromosome X, Y, 2L, 2R, 3L, 3R and 4. Additional conditions can be added in `bioawk` or can be executed with `grep`. Because there are many "4" in the strings, we have to anchor 4 at the boundary if `grep` is used.. 

    bioawk -c gff ' $3=="gene"&&($1=="X"||$1=="Y"||$1=="2L"||$1=="2R"||$1=="3L"||$1=="3R"||$1=="4")  { print $1 }' <dmel-all-r6.48.gtf | sort| uniq -c|sort -h|less

or

    bioawk -c gff ' $3=="gene" { print $1 }' <dmel-all-r6.48.gtf | sort | uniq -c | sort -k 2 | grep -E 'X|Y|2L|2R|3L|3R|\b4' |less

The prior commands give the following results:

    3515 2L
    3653 2R
    3489 3L
    4227 3R
    114 4
    2708 X
    113 Y



### **Answers to the questions:**

1. Total number of features of each type, sorted from the most common to the least common

|Number|feature type|
|---:|:---|
|32 |snRNA|
|115 |rRNA|
|262 |pre_miRNA|
|300 |snoRNA|
|312 |tRNA|
|365 |pseudogene|
|485 |miRNA|
|3053 |ncRNA|
|17896 |gene|
|30799 |mRNA|
|30825 |stop_codon|
|30885 |start_codon|
|33738 |3UTR|
|46802 |5UTR|
|163242 |CDS|
|190050 |exon|

2. Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)

|Number of genes|chromosome arm|
|---:|:---|
|3515 |2L|
|3653 |2R|
|3489 |3L|
|4227 |3R|
|114 |4|
|113 |Y|
|2708 |X|