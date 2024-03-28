#!/bin/bash

wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz  # Archaea 16S
wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz  # Bacteria 16S
cat *rRNA.fna.gz > ncbi16s18srRNA.bacteria_archaea.fna.gz 
gunzip ncbi16s18srRNA.bacteria_archaea.fna.gz 
minimap2 -x map-ont -d ncbi16s18srRNA.bacteria_archaea.mmi ncbi16s18srRNA.bacteria_archaea.fna
wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
grep -oP '^>[^\s]+' ncbi16s18srRNA.bacteria_archaea.fna | sed 's/^>//' > accessions.txt
zcat nucl_gb.accession2taxid.gz | grep -w -f accessions.txt > hit.acc2taxid.tsv
cut -f2,3 hit.acc2taxid.tsv > ref2taxid.targloci.tsv
gzip ncbi16s18srRNA.bacteria_archaea.fna
# rm archaea.16SrRNA.fna.gz bacteria.16SrRNA.fna.gz nucl_gb.accession2taxid.gz hit.acc2taxid.tsv accessions.txt
