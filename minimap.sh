#!/bin/bash

for FILE in `ls -1 *_trimmed.fastq.gz | sed 's/_trimmed.fastq.gz//'`
do
    echo "minimap2 for ${FILE}"
    minimap2 -ax map-ont -t 24 ../databases/ncbi16s18srRNA.bacteria_archaea_fungi.mmi ${FILE}_trimmed.fastq.gz > ${FILE}.sam
    echo "doing samtools"
    samtools view -o ${FILE}.bam ${FILE}.sam
    samtools sort -o ${FILE}.sorted.bam ${FILE}.bam
    samtools index ${FILE}.sorted.bam
    samtools view -F 0x900 -o ${FILE}_primary.sorted.bam ${FILE}.sorted.bam
    samtools index ${FILE}_primary.sorted.bam
    seqkit bam -C ${FILE}_primary.sorted.bam 2> ${FILE}_prialignment.tsv
done