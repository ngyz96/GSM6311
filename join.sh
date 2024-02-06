#!/bin/bash

for DIRECTORY in `ls -1 -d`;
do
	cd ${DIRECTORY};
	cat *.fastq.gz > ../../combined/${DIRECTORY}.fastq.gz;
	cd -;
done