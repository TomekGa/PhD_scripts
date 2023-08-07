#creates de novo RNA assembly

#REQUIRES:
# - path to raw reads
# - installed TRINITY software (either via mamba, via conda or from source)
# if installed from source it might be needed to specify the instalation location of other components e.g. 
#os.system('export PATH="$PATH:/home/tomasz.gaczorek/Transcriptomes/jellyfish-2.3.0/bin"')
#os.system('export PATH="$PATH:/home/tomasz.gaczorek/Transcriptomes/Salmon-latest_linux_x86_64/bin"')

#OUTPUT 
#de novo assembly

#RUNNING
#python3 trinity_tool.py ~/path_to_raw/

import os
import sys

path_in = sys.argv[1]

def list_files1(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith('.' + extension)]

fastq_list = list_files1(path_in,'fastq')
fastq_list.sort()
print(fastq_list)

L = []
R = []
for i in range(len(fastq_list)):
	if i % 2 == 0:
		whole_path = path_in+'/'+fastq_list[i]
		L.append(whole_path)
	else:
		whole_path = path_in+'/'+fastq_list[i]
		R.append(whole_path)

L = ','.join(L)
R = ','.join(R)
		
os.system('/home/tomasz.gaczorek/Transcriptomes/trinityrnaseq-v2.8.6/Trinity --seqType fq --max_memory 450G --left '+L+' --right '+R+' --CPU 60 --trimmomatic --quality_trimming_params "ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:75"')

#if new conda env created, remember to install also bowtie (not bowtie2)