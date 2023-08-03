#!/bin/bash

#you must have seqkit installed via conda
#source /home/tomasz.gaczorek/miniconda3/etc/profile.d/conda.sh   # ADJUST and run if needed (bash might not see programs inside activated conda environment)

while getopts "i:p:m:F:R:" option
do
case "${option}"
in
m) mhc_type=$OPTARG;;
i) IDs=$OPTARG;; #file with IDs
p) path_var=$OPTARG;; # path to the folder with reads (do not end with '/')
# ^ /mnt/matrix/projects/backup_miseq/201105_M01530_0062_000000000-J9L75/Data/Intensities/BaseCalls
F) forward=$OPTARG;; #forward primer
R) reverse=$OPTARG;; #reverse primer
esac
done

mkdir "MHC_$mhc_type"
# copying proper files
while IFS= read -r line #read the file line by line
do
ID="$(echo $line | cut -f1 -d,)" # get ID
echo "$ID"
$(cp ${path_var}/*$ID*-${mhc_type}_*R*.fastq.gz "./MHC_${mhc_type}/") # copy
done < "$IDs"
cd "MHC_$mhc_type"

#merging reads
mkdir merged
gzip -d *.gz #unzip
for f in *_R1_001.fastq
    do
    left=$f
    right=$(echo $f | sed s/_R1_/_R2_/)
    output="./merged/$(echo $f | sed s/_R1_001.fastq/_merged/)"
    echo $left
    echo $right
    echo $output
    pear -f $left -r $right -o $output -v 5 -j 20 #mearging
    done
rm *.fastq
cd merged
ls | grep -v assembled.fastq | xargs rm

#trimming
conda activate SeqKit
mkdir ../trimmed
reverse=$(echo $reverse | tr ACGTMRWSYKVHDBN TGCAKYWSRMBDHVN | rev)
primers_length="$((${#forward}+1)):-$((${#reverse}+1))"
echo "Chosen range: $primers_length"
for f in $(ls *.fastq)
do
echo $f
base=$(echo ${f} | cut -f1 -d'_')
#cutadapt --cut 18 --cut -18 -m 200 -o "../trimmed/"${base}"_trimmed.fastq" ${f} #trimming
seqkit amplicon ${f} -F $forward -R $reverse -r $primers_length -o "../trimmed/${base}_trimmed.fastq"
done
conda deactivate
cd ../trimmed
rm -r ../merged
zip -m "../MHC_${mhc_type}.zip" *
#zip -m ../MHC_II.zip $(ls | grep II)
cd ..
rm -r ./trimmed

#running amplisas - ADJUST path to the ampliSAS scripts
if [ "$mhc_type" == "I" ]; then
perl ~/amplisat/ampliSAS_mhci_old.pl -i "MHC_${mhc_type}.zip" -o ./ -t Illumina -thr 17
else perl ~/amplisat/ampliSAS_mhcii.pl -i "MHC_${mhc_type}.zip" -o ./ -t Illumina -thr 17
fi
