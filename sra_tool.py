#download raw reads from SRA based on the text file with accessions

import os

#it requires:
# #path to your extract SRAtoolkit folder
path_SRAtoolkit = '/home/tomasz.gaczorek/Transcriptomes/sratoolkit.2.10.2-centos_linux64/'
# #path to list of accession numbers dowloaded from NCBI
path_accessions = './SraAccList.txt'

f1 = open(path_accessions,"r")
acc_list = f1.readlines()
f1.close()

for i in range(len(acc_list)):
	given_number = acc_list[i].strip()
	os.system(path_SRAtoolkit+"bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files "+given_number)
	#os.system('mv ./*.fastq /mnt/matrix/temp/Tomek/podacris/cretensis_single')

