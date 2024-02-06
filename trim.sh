#!/bin/bash

for FILE in `ls -1 *.fastq.gz | sed 's/.fastq.gz//'`
do
    seqkit seq ${FILE}.fastq.gz -M 1900 -m 1350 -j 20 -o ${FILE}_trimmed.fastq.gz
done