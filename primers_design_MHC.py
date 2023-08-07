###remember to have fasta file in format not splited to 60 characters rows !!!

#Tomek Gaczorek
#tomek.gaczorek@gmail.com

#script to design multiple primers for highly variable regions (e.g. MHC) 
#it tries to provide the best combination of primers based on desired primer length, maximum degeneration and shift of achoring region  
#it works well for big number of sequences
#note that it might not provide the optimal primer combination in all cases as the algorithm cannot check all combinations
#the script ignores all sequences with gaps in the analysed window

#REQUIRES:
# - aligned fasta file - each second line must be a sequence (make sure the sequence is not splitted into multiple rows)
# - beginning and end of the marker

#OUTPUT:
# text files with the optimal primer combinations for each primer position 

#RUNNING
# python3 primers_design_new_approach.py alignment.fasta <start> <end>

import sys
from functools import reduce
import math
from more_itertools import distinct_combinations

#filename = "MHC_II_references_Proteus_vs_Ambystoma_exon2.fas"
filename = sys.argv[1]
#start = 148
#end = 418 
#MHC I
start = int(sys.argv[2])
end = int(sys.argv[3])
min_primer_length = 18
max_primer_length = 18
shift = 30 #ADJUST - maximum shift of primer from the start and from the end
combinations_length = [i for i in range(min_primer_length,max_primer_length+1)]
combinations_startpoint_start = [i for i in range(start-1,(start+shift+1))]
combinations_startpoint_end = [i for i in range((end-shift-max_primer_length),(end-min_primer_length+1))]
#combinations_startpoint_end = [end - max_primer_length]
#print(len(combinations_startpoint_end))
#print(len(combinations_startpoint_start))
#print(len(combinations_length))
max_degeneration = 12 #ADJUST
min_templates_number = 1
primers_types = ['start','end'] #ADJUST - it's possible to design only F or R primers
#primers_types = ['end']

###############################
		
class Sequence(object):
	def __init__(self,number,start,stop,sequence):
		self.number = number
		self.start = start
		self.stop = stop
		self.bases = []
		self.get_bases_from_string(sequence)
		
	def get_bases_from_string(self,string):
		for i in string:
			self.bases.append(i)
			
class Combination_place(object):
	def __init__(self,number,location):
		self.number = number
		self.location = location
		self.sequences = []
		#self.primers_bases = []
		
class Newton_combination(object):
	def __init__(self,degree):
		self.degree = degree
		self.newton_sequences = []

####################################



def unique(list1): 
	unique_list = []
	for x in list1: 
		if x not in unique_list: 
			unique_list.append(x) 
	return unique_list
	
def get_rid_of_hyphae(list1):
	hyphae_list = []
	for x in list1:
		if '-' not in x:
			hyphae_list.append(x)
	return hyphae_list

def get_seqs_from_file(file_name):
	input_sequences = []
	fasta_file = open(file_name,'r').readlines()
	count = 0
	for i in range(len(fasta_file)):
		if (i+1) % 2 == 0:
			count = count + 1
			fasta_file[i] = fasta_file[i].strip()
			input_sequences.append(fasta_file[i])
	return input_sequences

def get_combinations_place(location):
	list_of_combinations = []
	if location == 'start':
		combinations_startpoint = combinations_startpoint_start
	else:
		combinations_startpoint = combinations_startpoint_end
	count = 0
	for i in range(len(combinations_length)):
		for j in range(len(combinations_startpoint)):
			count = count + 1
			comb = Combination_place(count,location)
			for k in range(len(input_sequences)):
				primer = Sequence(k+1,combinations_startpoint[j],combinations_startpoint[j]+combinations_length[i],input_sequences[k][combinations_startpoint[j]:(combinations_startpoint[j]+combinations_length[i])])
				comb.sequences.append(primer)
			comb.seq_length = len(comb.sequences[0].bases)
			comb.seq_start = comb.sequences[0].start
			list_of_combinations.append(comb)
	return list_of_combinations

def check_for_duplicates(combinations):
	to_maintain = []
	for i in combinations:
		sequences = []
		for j in i.sequences:	
			sequences.append(''.join(j.bases))
		unique_seq = unique(sequences)
		unique_seq = get_rid_of_hyphae(unique_seq)
		preserved = []
		for l in unique_seq:
			index = sequences.index(l)
			preserved.append(index)
		new_sequences = []
		for j in range(len(i.sequences)):
			if j in preserved:
				new_sequences.append(i.sequences[j])
		i.sequences = new_sequences
		i.number_of_templates = len(i.sequences)
		if i.number_of_templates >= min_templates_number:
			to_maintain.append(i.number)
	new_combinations = []
	for i in to_maintain:
		for j in combinations:
			if i == j.number:
				new_combinations.append(j)
	return new_combinations
	
def get_right_letter(list1):
	good_list = []
	for i in list1:
		if i != '-':
			good_list.append(i)
	if len(good_list) == 0:
		letter = 'N'
	elif len(good_list) == 1:
		letter = good_list[0]
	elif len(good_list)==2:
		if 'A' in list1 and 'G' in list1:
			letter = 'R'
		elif 'A' in list1 and 'C' in list1:
			letter = 'M'
		elif 'C' in list1 and 'G' in list1:
			letter = 'S'
		elif 'A' in list1 and 'T' in list1:
			letter = 'W'
		elif 'C' in list1 and 'T' in list1:
			letter = 'Y'
		elif 'G' in list1 and 'T' in list1:
			letter = 'K'
	elif len(good_list)==3:
		if 'A' not in list1:
			letter = 'B'
		elif 'C' not in list1:
			letter = 'D'
		elif 'G' not in list1:
			letter = 'H'
		elif 'T' not in list1:
			letter = 'V'
	elif len(good_list)==4:
		letter = 'N'
	return letter
        
def check_degeneration(list1):
	#print(list1)
	numbers = []
	for i in list1:
		if i in ['A','C','G','T']:
			numbers.append(1)
		elif i in ['R','M','S','W','Y','K']:
			numbers.append(2)
		elif i in ['B','D','H','V']:
			numbers.append(3)
		else:
			numbers.append(4)
	outcome = reduce((lambda x,y: x*y),numbers)
	return outcome
			
def get_primers(newton):
	upper_names_list = []
	for j in range(len(newton.newton_sequences[0].bases)):
		names_list = []
		for k in newton.newton_sequences:
			names_list.append(k.bases[j])
		names_list = unique(names_list)
		if len(names_list) > 1:
			upper_names_list.append(get_right_letter(names_list))
		else:
			upper_names_list.append(names_list[0])
	newton.primers_bases = upper_names_list
	newton.degeneration = check_degeneration(upper_names_list)

def newtons(i):
	global do_it
	for j in range(i.number_of_templates):
			newtons_list = []
			degeneration_step = i.number_of_templates-j
			#print(degeneration_step)
			print('check 1')
			combinatorix = list(distinct_combinations(i.sequences,degeneration_step))
			new_combinatorix = []
			for k in combinatorix:
				n_x = []
				for l in k:
					n_x.append(l)
				new_combinatorix.append(n_x)
			combinatorix = new_combinatorix
			###remember - if all or 1, no lists inside other lists!!!!!!
			print('check 2')
			for k in combinatorix:
				newton = Newton_combination(degeneration_step)
				newton.newton_sequences = k
				newtons_list.append(newton)
			dagen_list = []
			dagen_list_values = []
			print('check 3')
			for k in newtons_list:
				get_primers(k)
				if k.degeneration <= max_degeneration:
					dagen_list.append(k)
					dagen_list_values.append(k.degeneration)
			if len(dagen_list_values)>0:
				minimum = min(dagen_list_values)
				index = dagen_list_values.index(minimum)
				given_newton = dagen_list[index]
				new_seq = []
				for l in i.sequences:
					if l not in given_newton.newton_sequences:
						new_seq.append(l) 
				out_out.append(given_newton)
				i.number_of_templates = len(new_seq)
				i.sequences = new_seq
				if len(i.sequences)==0:
					do_it = False
				break
				

##########################################

input_sequences = get_seqs_from_file(filename)
for ptype in primers_types:
	combinations = get_combinations_place(ptype)
	combinations = check_for_duplicates(combinations)
	if len(combinations) < 1:
		sys.exit('No combinations fulfilled prerequisites')
	outer = []
	for ii in combinations:
		print('combination '+str(ii.number))
		do_it = True
		out_out = []
		count = 0
		while do_it == True:
			count = count+1
			print(count)
			newtons(ii)
		out_out_len = len(out_out)
		for m in out_out:
			title = '####'+str(out_out_len)+'####'+'>Combination nr '+str(ii.number)+' Start '+str(ii.seq_start)+' Length '+str(ii.seq_length)+' Position '+ ii.location+' Degeneration '+ str(m.degeneration) +'\n'
			seq = ''.join(m.primers_bases)+'\n'
			outer.append(title)
			outer.append(seq)
		
	fasta_file_name = 'primers_code_'+ptype+'.txt'
	fasta_file = open(fasta_file_name,'w')
	fasta_file.write(''.join(outer))
	fasta_file.close()				
				
		
