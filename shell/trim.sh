#!/bin/bash

for FILE in `ls -1 *.fastq.gz | sed 's/.fastq.gz//'`
do
    # 1134 1740 based on NCBI calculation, 1300-1800 based on paper 
    gunzip -c ${FILE}.fastq.gz | chopper -q 10 --minlength 1300 --maxlength 1800 --threads 12 | gzip > ${FILE}_filtered.fastq.gz;
done
