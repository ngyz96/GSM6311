#!/bin/bash

for DIRECTORY in `ls -1`;
do
	cd ${DIRECTORY};
	cat *.fastq.gz > ${DIRECTORY}.fastq.gz;
	cp ${DIRECTORY}.fastq.gz ../${DIRECTORY}.fastq.gz
	cd -;
done