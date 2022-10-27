# Project Analysis Proposal

-------


### Goals and data source.

In this research, we aim to estimate the effect from types of genetic variations on population statistics, especially on selection signals. We use *Drosophila melanogaster* as our research system. Because there are purly isogenic inbred lines, the influence from demographic history can be eliminated easily. The research materials, DNA sequences, come form public datasets, such as [DSPR](https://wfitch.bio.uci.edu/~dspr/index.html) and [Flybase](http://flybase.org/). When the pipeline is built well, sequences from other database, including [DGRP](http://dgrp2.gnets.ncsu.edu/) and NCBI (PRJNA418342, PRJNA618654, PRJNA559813), will be tested. The sequence of *D. simulans*, acting as an outgroup, is also acquired from NCBI. Thus, the SNPs and structural variations (SVs) are polarized. 

### Preparation of data for analysis

The process of analysis is described as following. First, short Illumina reads from DSPR are filtered according to quality and mapped against latest *D. melanogaster* reference genome with **BWA-aln**. Quality control is conducted after mapping again with **samtools**. Then, the SNPs are called with **Freebayes**, while **GATK** is another option. Crude SNPs are filtered with **vcffilter** and finally annotated by **SnpEff**. After annotaing the SNPs 
(and then SVs), the key step is to separate the SNPs from coding regions, non-coding genes, introns, 5'UTRs, 3'UTR, promoters, and intergenic regions. To process the vcf file, a customed bash scripts with **bedtools** is modified from existing work. All the process for calling structural variations are similar to calling SNPs. After mapping and quality control, **SVMU** is used for SV calling. The annotation can be done by **CNVnator** & **nanotatoR** or, **AnnotSV**. The SVs are classified into short indels, long indels \( \>200 bp \), copy number variants \( \>1 copy \), inversions, and translocations. After classification, we can perform downstream analysis. 

### Statistical analisis and visualization

By isolating SNPs \( and then SVs\) in different functional regions, we can calculate the pi, theta-w, Tajima's D and mutation rate. These values are calculated right at the sites and in both upstream plus down streams by sliding windows of 5 kbp with 1 kbp increment. The statistics across the genome is also calculated by the sliding wndow method. To visualize our results, **ggPlot2** is the main package to draw line graphs, bar charts and box plots. The statistics across a certain region are presented as line graphs. To compare the distributions, box plots representing the value of each group are created. The site frequency spectrum can be extracted from vcf files with **vcftools** and **bcftools**. With different input, we can generate the site frequency spectrum in different regions or conditioned by different Tajima's D range. A special package AlleleShift can also be used to draw site frequency spectrum.  

### Conclusion
There are many similar graphs in research articles (Daniel K Fabian et al., 2012, Mol Ecol; Stefano Mona et al., 2010, BMC Evolutionary Biology) so we can take their graphs and methods as references. With current scripts, we are able to process the sequences and will spend some time on visualization. After visualizing our data, we can compare the distribution and trends of population statistics under different conditions. Statistical tests, such as &Chi;<sup>2</sup> test and linear regression, may help. Besides, we can understand how the SNPs influence the selection and  population statistics in neighboring regions. When the pipelines is sophiscated, similar analysis will be conducted for structural variations.


---------------------

### List of sofeware and packages

| purpose | names of softwares and packages |
|--------|----------|
| genome mapping| bwa aln/sampe, minimap2 |
| SNP calling & annotation| GATK, Freebayes, SnpEff |
| SV calling & annotation| SVMU, CNVnator, nanotatoR, AnnotSV |
| population statistics| Egglib-3, Arlequin-3.5, poppr in R, vcf2sfs in R |
| Visualization | ggPlot2 in R, AlleleShift in R |
| miscelleneous | vcf2bed of BEDOPS, histlib, vcftools, bcftools, vcffilter, bedtools |
