
### There were five scripts involved. Here, we will explain them one by one.

## mapping.sh

One Anaconda environment includinf all required packages was used. First, we had to index the genome for `bwa` and `samtools`. The index files were placescollaterally with the reference file. Second, in a loop, we finished all the mapping process. Everything started from **fastq** files, named as `strain_lib_orientation.fq`. E.g., `A8_1_1.fq`. We extracted each elements of the **fastq** to name following files. `fastp` was used for cleaning up the reads, removing primers and generating an **html** report. The processed **fastq** files were passed to `bwa aln` so **sai** files were generated. Then, `bwa sampe` finished the basic mapping while adding  `@RG` information is essential for merging **sam** files later. **Sam** files from `bwa` were compressed, sorted, marked duplicates and filtered with `samtools` and then `qualimap bamqc`. Finally, three **bam** files of each strain were merged into one. Resulted 15 **bam** alignments were used in following applications.

## call.annotation.sh

This sript was basicallycomposed of two loops. The first for SNP calling and the second for annotation. We set "-g 200" in `freebayes` to avoid hightly repetitive resions which caused unusually high mapping read numbers. Before annotation, SNPs were filtered again with `rtg` and `vcffilter` to discard SNPs with poor quality scores. Default setting of `snpEff` were used. 

We tried another way to call variants. All the **bam** files were merged together at the beginning and then it was indexed (essential). Since the aligned reads were too huge, paralle version of `Freebayes` were used or the program may fail prematurly. `rtg vcffilter --all-samples -d 10 ` was a new part of the command to process all strains and to only allow reads with more than 10X depth. 

## coordinate.sh

We abtained the coordinates of each feature from a annotation file (**gff**) at Flybase when the position was minus one because of 1-based to 0-based transition. However, some features overlap with others so we had to process them to get regions with only one type of feature. Also, we only investigated euchromatic regions. Five prime untranslated regions (five_prime_UTR), three prime untranslated regions (three_prime_UTR), codons (CDS) and introns sometimes pverlap so the **bed** file of one of the features was substracted (`bedtools substract`) by all the other feature and then intersected with ranges of euchromatin (`bedtools intersect`.) For introns, we only included short (<65 bp) or long (>120 bp) introns when miltple times of `bedtools substract` function were used. The intergenic regions were divided into distal (>4000 bp away) or proximal ones. 

## calssifying.sh

With the coordinates (**bed** files) retirved from the previous script, we could intersect them with a **vcf** file. Sonymous variants and missense variants were differentiated in coordinates of CDS. 


## population_stats.sh

Here, we calculated theta, pi and Tajima's D with `vcftools` or `vcf kit`. To obtain the site frequency spectrum, `vcftools --freq2` was first used to extract the allele frequency and `bioawk` was used to organize the output.