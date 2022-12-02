#!/bin/bash

conda activate ee282
cd /pub/jenyuw/classrepos/EE282
minimap2 -t 60 -x ava-ont iso1_onp_a2_1kb.fastq iso1_onp_a2_1kb.fastq > overlaps.paf
miniasm -f iso1_onp_a2_1kb.fastq overlaps.paf > iso1_self_assemble.gfa
#transformation
awk '/^S/{print ">"$2"\n"$3}' iso1_self_assemble.gfa | fold > iso1_self_assemble.fa

faSize -detailed iso1_self_assemble.fa | gawk ' {n=n+$2; print $2} END { print n; }'|sort -rn|\
gawk ' NR==1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '

#installing the dependency
conda install -c conda-forge perl
conda install -c bioconda perl-bioperl
#making contigs
#the legands are given by aplhabetic orders istead of its original assignemnt
faSplitByN dmel-all-chromosome-r6.48.fasta dmel-all-contig.fa 20

faSize -detailed iso1_self_assemble.fa|sort -k 2,2 -n -r |gawk 'BEGIN {print "Name \t Length \t Assembly"} {print $1 "\t" $2 "\t" "A_iso1_self_assemble"}' >iso1_self_assemble_length.txt

faSize -detailed dmel-all-chromosome-r6.48.fasta|sort -k 2,2 -n -r |gawk 'BEGIN {print "Name \t Length \t Assembly"} {print $1 "\t" $2 "\t" "B_NCBI_chromosome"}' >r648_chr_length.txt

faSize -detailed dmel-all-contig.fa|sort -k 2,2 -n -r |gawk 'BEGIN {print "Name \t Length \t Assembly"} {print $1 "\t" $2 "\t" "C_NCBI_contig"}' >r648_contig_length.txt

plotCDF2 iso1_self_assemble_length.txt r648_chr_length.txt r648_contig_length.txt assembles.png
#plotCDF3 iso1_self_assemble_length.txt r648_chr_length.txt r648_contig_length.txt assembles3.png

conda install -c bioconda busco=5.4.3
busco --list-datasets
busco -i dmel-all-chromosome-r6.48.fasta -o r648_BUSCO -m genome --cpu 20 --auto-lineage-euk -f
#using specific database according to the original suggestions and results from automatic search. 
busco -i iso1_self_assemble.fa -o iso1_BUSCO2 -m genome --cpu 20 -l diptera_odb10

nucmer --mumreference --prefix=r648-chr_contig dmel-all-chromosome-r6.48.fasta dmel-all-contig.fa

#delta-filter -g r648-chr_contig.delta > r648-chr_contig-filtered.delta
#delta-filter -r -u 50 r648-chr_contig.delta > r648-chr_contig-filtered.delta
delta-filter -r -u 10 r648-chr_contig.delta > r648-chr_contig-filtered.delta

#mummerplot --png --large --prefix=r648-chr_contig r648-chr_contig-filtered.delta
mummerplot --png --large --prefix=r648-chr_contig_2 -R r648.Rfile r648-chr_contig-filtered.delta