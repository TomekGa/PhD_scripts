#!/bin/bash

#Runs blast search with DNA query
#REQUIRES installed blast
#OUTPUT - blast output in pairwise and tabular format

#RUN
#bash BLAST_script_nucl.sh -i input.fasta -q query.fasta -k TRUE -o FALSE

while getopts "i:q:k:o:" option
do
case "${option}"
in
i) in_name=$OPTARG;; #input file
q) query_name=$OPTARG;; #query file
k) keep_db=$OPTARG;; #'FALSE' if delete created databases 
o) only_BLAST=$OPTARG;; #'FALSE' if not making BLAST on existing database - then name of database is the input file
esac
done

db_name=${in_name%.*}
out_name=$db_name'_BLAST'
tabular_name=$db_name'_tabular.txt'
pairwise_name=$db_name'_pairwise.txt'

if [ "$only_BLAST" == "FALSE" ]; then
makeblastdb -in $in_name -out $db_name -parse_seqids -dbtype 'nucl'
else db_name=$in_name
fi

blastn -db $db_name -query $query_name -evalue 10e-10 -outfmt 11 -word_size 9 -gapopen 3 -gapextend 2 -reward 1 -penalty -1 -out $out_name #ADJUST

blast_formatter -archive $out_name -outfmt "7 sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out $tabular_name #ADJUST

blast_formatter -archive $out_name -outfmt 0 -out $pairwise_name #ADJUST

if [ "$keep_db" == "FALSE" ]; then
find . -type f -name "*.n*" -exec rm -f {} \;
fi

