#Retrieve sequences from assembly base on pairwise BLAST output (based on DNA)
#The code controls for the strand direction
#Note it relies on specific file naming convention - core_pairwise.txt (BLAST) & core.fasta (assembly) - both files in the same directory
#It requires additional script "fasta_id_recov_LAST.py"
#You can choose desired sequences after the code is executed (prompt) by using similar convention as in printing e.g., 1-10,17

#RUN
#python3 core_pairwise.txt

import sys
import os

filename_pairwise = sys.argv[1]
#filename_pairwise = 'total_assembly_pairwise.txt' #to be deleted - previous will be kept
assembly_filename = filename_pairwise.split('_pairwise.txt')[0]+'.fasta'
d_pairwise = open(filename_pairwise,'r').readlines()

def get_seq_numbered(data_pairwise):
	given_line = len(data_pairwise)
	number = 0
	fasta_number = 0
	number_of_queries = 0
	beginnings = []
	for i in range(len(data_pairwise)):
		if data_pairwise[i].startswith('Sequences producing significant alignments'):
			given_line = i
			number_of_queries = number_of_queries + 1
			beginnings.append(i)
			fasta_number = 0
		if i > (given_line+1) and data_pairwise[i] != '\n':
			number = number+1
			num_write = "{:03d}".format(number)
			data_pairwise[i] = str(num_write)+' '+data_pairwise[i]
		if data_pairwise[i] == '\n' and given_line < len(data_pairwise) and i != (given_line+1):
			given_line = len(data_pairwise)
			number = 0
		if data_pairwise[i].startswith('>'):
			fasta_number = fasta_number+1
			data_pairwise[i-1] = str(fasta_number)+'\n'
	beginnings.append(len(data_pairwise))
	return [data_pairwise,number_of_queries,beginnings]

output = get_seq_numbered(d_pairwise)
number_of_queries = output[1]
beginnings = output[2]

out_file = open(filename_pairwise.split('_pairwise.txt')[0]+'_numbered.txt','w')
out_file.write("".join(output[0]))
out_file.close()

filename_pairwise_changed = filename_pairwise.split('_pairwise.txt')[0]+'_numbered.txt'
d_pairwise1 = open(filename_pairwise_changed,'r').readlines()

for k in range(number_of_queries):
	template_text = 'Look at the file called: '+filename_pairwise_changed+'\nLines: '+str(beginnings[k])+'-'+str(beginnings[k+1])+'\n'
	prompt_text = template_text+'Which sequences do you want to maintain for query #'+str(k+1)+'? (use commas and hyphens) '
	wektor = input(prompt_text)
	wektor = wektor.split(",")
	
	indexes_of_seq = []
	for i in wektor:
		if i.count('-') == 0:
			indexes_of_seq.append(int(i))
		else:
			ranging = i.split('-')
			ranging = [j for j in range(int(ranging[0]),int(ranging[1])+1)]
			indexes_of_seq = indexes_of_seq + ranging


	d_pairwise = d_pairwise1[beginnings[k]:beginnings[k+1]]
	fasta_output = []
	reversing = ''			
	for i in range(len(d_pairwise)):
		if d_pairwise[i].startswith('>'):
			number = int(d_pairwise[i-1].split('\n')[0])
			if number in indexes_of_seq:
				seq_string = ''
				j = i
				while not d_pairwise[j].startswith('Length'):
					seq_string = seq_string + d_pairwise[j].strip('\n')
					j = j+1
				seq_string = seq_string+'\n'
				fasta_output.append(seq_string)
				reversing_label = d_pairwise[j+4].split('=')[1]
				if reversing_label == 'Plus/Minus\n':
					reversing = reversing+'1'
				else:
					reversing = reversing+'0'

	fasta_output.insert(0,reversing+'\n')
	out_file = open('temporary_IDs.txt','w')
	out_file.write("".join(fasta_output))
	out_file.close()
	
	os.system('python3 fasta_id_recov_LAST.py '+assembly_filename+' temporary_IDs.txt '+str(k+1))
	new_file_name = assembly_filename.split('.fasta')[0]+'_given_IDs.fasta'
	newer_file_name = assembly_filename.split('.fasta')[0]+'_given_IDs_query'+str(k+1)+'.fasta'
	os.system('mv '+new_file_name+' '+newer_file_name)

	
