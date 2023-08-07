#perform mapping of your original reads to reference sequences

#REQUIRES:
# - fasta file with reference sequences
# - path to directory with raw reads

#OUTPUT
# - sorted bam
# - index file
# (can be visualised with Tablet software)

#RUNNING
# python3 raw_mapping.py reference.fasta ~/path_to_raw/

import os
import sys

reference_file = sys.argv[1]
trans_path = sys.argv[2]
align_type = 'very-sensitive'
#in_dir = 'Mapping_files'
#os.system('cd '+in_dir)
#os.system('indexes_compressed=$(ls *.fastq* | uniq')
#os.system('printf "%s\n" "$indexes_compressed"')
#os.system('gzip -d $indexes_compressed')

def list_files1(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith('.' + extension)]
 

#cwd = os.getcwd()  
#os.chdir(cwd+'/'+in_dir)
#cwd = os.getcwd() 
fastq_list = list_files1(trans_path,'fastq')
fastq_list.sort()
print(fastq_list)

L=[]
L_cut_names = []
R=[]
for i in range(len(fastq_list)):
	name_cut = fastq_list[i].split('.fas')[0]
	if (i+1) % 2 == 0:
		#new_name = name_cut+'_R_.fastq'
		new_name = fastq_list[i] #######################
		#os.rename(fastq_list[i],new_name)
		R.append(trans_path+'/'+new_name)
	else:
		#new_name = name_cut+'_L_.fastq'
		new_name = fastq_list[i] #######################
		#os.rename(fastq_list[i],new_name)
		L.append(trans_path+'/'+new_name)
		L_cut_names.append(name_cut)
		
#L=','.join(L)
#R=','.join(R)

### creating reference file ###
os.system('bowtie2-build ./'+reference_file+' ref_idx') #creates 6 files - compress sequence information for faster alignment
os.system('samtools faidx ./'+reference_file)

out_names = []
for i in range(len(L)):
	out_name_sam = L_cut_names[i]+'.sam'
	out_name_bam = L_cut_names[i]+'.bam'
	out_name_bam1 = L_cut_names[i]+'_sorted.bam'
	os.system('bowtie2 --'+align_type+' --no-unal --no-mixed -q --phred33 -p 60 -x ref_idx -1 '+L[i]+' -2 '+R[i]+' -S '+out_name_sam+' --rg-id '+L_cut_names[i]+' --rg SM:'+L_cut_names[i])
	os.system('samtools view -bS '+out_name_sam+' > '+out_name_bam)
	os.system('samtools sort '+out_name_bam+' -T aaa -o '+out_name_bam1)
	os.system('samtools index '+out_name_bam1)
	#out_names.append(out_name_bam1)
	os.system('rm '+out_name_bam)

os.system('samtools merge all.bam *.bam')
os.system('samtools sort all.bam -T aaa -o all_sorted.bam')
os.system('samtools index all_sorted.bam')