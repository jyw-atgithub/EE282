#!/usr/bin/env bash

srun -A class-ee282 --pty --x11 bash -i
conda activate ee282
cd cd /data/homezvol2/jenyuw/classrepos/EE282

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt

md5sum -c md5sum.txt|grep "dmel-all-chromosome-r6.48" |less
# md5sum -c <(grep "dmel-all-chromosome-r6.48" md5sum.txt)
# md5sum -c md5sum.txt | grep dmel-all-chromosome-r6.48.fasta.gz

faSize dmel-all-chromosome-r6.48.fasta.gz |less
