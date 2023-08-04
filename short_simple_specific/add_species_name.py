#Adds species (or any other) name to each sequence header in fasta file

#RUN
#python3 add_species_name.py file.fasta species_name

import os
import sys

file_to_use = sys.argv[1]
given_name = sys.argv[2]

data = open(file_to_use,'r').readlines()

new_data = []
for i in range(len(data)):
	if data[i].startswith('>'):
		new_line = data[i].split('>')
		#print(new_line)
		new_data.append('>'+given_name+' '+new_line[1])
	else:
		new_data.append(data[i])
		
changed_file = open(file_to_use,'w')
changed_file.write(''.join(new_data))
changed_file.close()
