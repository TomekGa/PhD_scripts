#!/bin/bash

#requires:
#1. paired .gz files with all reads from a single run
#2. descriptive file for amplisas with primers, individual IDs and tags
#3. installed pear software for mearging overlapping reads

#example files Run1_R1.gz & Run1_R2.gz & Run1_descriptive.csv

#decompress
gzip -d Run1_R1.gz
gzip -d Run1_R2.gz

#merge
pear -f Run1_R1.fastq -r Run1_R2.fastq -o Run1.fastq -v 5 -j 20
# ^ v - overlap parameter (see pear documentation)
# ^ j - number of cores used

#convert to fasta
sed -n '1~4s/^@/>/p;2~4p' Run1.fastq > Run1.fasta
rm *.fastq #clean up

#running amplisas

perl ~/amplisat/ampliSAS.pl -i Run1.fasta -d Run1_descriptive.csv -o ./ -t Illumina -thr 17


