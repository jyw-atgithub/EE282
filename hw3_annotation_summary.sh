#!/bin/bash
srun -A class-ee282 --pty --x11 bash -i
conda activate ee282
cd /data/homezvol2/jenyuw/classrepos/EE282

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz
curl http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt -0 gtfmd5sum
md5sum -c gtfmd5sum

#gunzip -k dmel-all-r6.48.gtf.gz
#cut -f 3 dmel-all-r6.48.gtf| sort |uniq -c |sort -k 1 -h|less
bioawk -c gff '{ print $3 }' <dmel-all-r6.48.gtf.gz|sort|uniq -c|sort -n|less

bioawk -c gff ' $3=="gene"&&($1=="X"||$1=="Y"||$1=="2L"||$1=="2R"||$1=="3L"||$1=="3R"||$1=="4")  { print $1 }' \
<dmel-all-r6.48.gtf | sort| uniq -c|sort -h|less

#bioawk -c gff ' $3=="gene" { print $1 }' <dmel-all-r6.48.gtf | sort | uniq -c | sort -k 2 | grep -E 'X|Y|2L|2R|3L|3R|\b4' |less

