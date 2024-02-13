#!/bin/bash

for FILE in `ls -1 *.fastq.gz | sed 's/.fastq.gz//'`
do
    # 1134 1740 based on NCBI calculation
    seqkit seq ${FILE}.fastq.gz -M 1900 -m 1350 -j 20 -o ${FILE}_trimmed.fastq.gz
done
