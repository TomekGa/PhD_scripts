#Converts text file with aligned sequences (.aln extention) into standard fasta

#RUN
#python3 file_to_convert.aln

file_name = sys.argv[1] #ADJUST
data = open(file_name,'r').readlines()
#print(data)

def filter_len(x):
		if len(x) > 1:
			return True
		else: return False

output_list = []
for i in range(len(data)):
	if i == 0 or (data[i].startswith('\n') and i != (len(data)-1)):
		counter = 0
		addition = 1
		if i == 0:
			ii = i
		else:
			ii = i+1
		for j in [k for k in range(ii,len(data))]:
			given_line = data[j].strip('\n')
			given_line = given_line.split(' ')
			given_line = list(filter(filter_len,given_line))
			#print(given_line)
			seq = given_line[1]
			#print(seq)
			if i == 0:
				output_list.append('>'+given_line[0])
				output_list.append(seq)
			else:
				output_list[counter+addition] = output_list[counter+addition]+seq
				counter = counter+1
				addition = addition+1
			if data[j+1].startswith('\n'):
				break

for i in range(len(output_list)):
	output_list[i] = output_list[i]+'\n'

new_file_name = file_name.split('.')[0]+'_converted.fasta'
changed_file = open(new_file_name,'w')
changed_file.write(''.join(output_list))
changed_file.close()

print(output_list)
