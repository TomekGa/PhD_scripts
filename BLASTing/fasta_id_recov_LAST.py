#Code B.4 fasta_id_recov.py
"""
Removes sequences having header not equal to ids given in id file
can even remove even any sequence
"""
usage = "Use:\nfasta_id_recov.py file_name.fasta idFile"
import sys
#import Bio

flanking = 50

def newFunc(allData):
    #function to remove new line character from string list
    seqsList= []
    
    lis = []
    string = ""
    idList = []
    for i in range(0,len(allData)):
        if allData[i][0]==">":
            idList.append(allData[i].strip())
            if string != "":
                seqsList.append(string)
                string = ""
        else:
            string = string.strip() + allData[i].strip()  
    
    seqsList.append(string)

    return [idList,seqsList]
#################################
# def seq_split(seq):
    # i = 0
    # output = ""
    # while i < len(seq):
        # output = output + seq[i:i + 60]+"\n"
        # i = i + 60
    # return output
#################################
def take_care_of_tabular(x,in_title,in_seq):
    checking = 0
    checking_query = 0
    points = []
    for i in range(len(x)):
        if x[i].startswith('#'):
            checking = checking + 1
        if checking == 0 and checking_query == number_of_query:
            splited = x[i].split('\t')
            if splited[1] != 'N/A':
                given_title = '>'+splited[0]+' '+splited[1]
            else:
                given_title = '>'+splited[0]
            if given_title == in_title:
                points.append(int(splited[8]))
                points.append(int(splited[9]))
                #in_seq = in_seq[start:end]
        if checking == 5:
            checking = 0
            checking_query = checking_query + 1
    start = min(points)+1-flanking
    if start < 1:
        start = 1
    #print(start)
    end = max(points)+1+flanking
    if end > (len(in_seq)-1):
        end = len(in_seq)-1
    #print(end)
    in_seq = in_seq[start:end]
    return in_seq
#################################
def get_ids(id_fname):
    i = 0
    idList = []
    f = open(id_fname,"r")
    data = f.read().split("\n")
    f.close()
    reversing = str(data[0])
    #print(reversing)
    data = data[1:]
    for i in data:
        if i != '':
            if i.startswith('>'):
                if i != data[0]:
                    idList.append(id_new.strip()[1:])
                id_new = i
            else:
                id_new = id_new + i
    idList.append(id_new.strip()[1:])
    return [idList,reversing]
#################################
def write_fasta(ofile, id_list, id_template, seqs_template):
    f2 = open(ofile,"w")
    tabular_file = filename.split('.fasta')[0]+'_tabular.txt' 
    tabular_d = open(tabular_file,'r')
    tabular_data = tabular_d.readlines()
    tabular_d.close()
    for i in range(len(id_list)):
        title = '>'+id_list[i]
        index = id_template.index(title)
        #print(index)
        #print(temp)
        #print(id_template[index])
        seq = seqs_template[index].strip()
        seq = take_care_of_tabular(tabular_data,title,seq)
        #print(reversing)
        if reversing[i] == '1':
            #print(seq)
            seq = list(seq)
            for j in range(len(seq)):
                #print(seq[j])
                if seq[j] == 'A':
                    seq[j] = 'T'
                elif seq[j] == 'T':
                    seq[j] = 'A'
                elif seq[j] == 'C':
                    seq[j] = 'G'
                else:
                    seq[j] = 'C'
            seq = ''.join(seq)
            seq = ''.join(reversed(seq))
            #print(seq)
        #print(seq)
        title = title+'\n'
        seq = seq+'\n'
        f2.write(title + seq)    
    f2.close()
    print("file written...")
    return i+1
#################################
if len(sys.argv)<3:
    print("Error in fetching data from input file...")
    print(usage)
else:
    filename =sys.argv[1]
    id_file = sys.argv[2]
    number_of_query = int(sys.argv[3])
    outfileName = filename[:filename.find(".fasta")]+"_given_IDs.fasta"
    f1 = open(filename,"r")
    xx = f1.readlines()
    f1.close()
    x = newFunc(xx)
    idList = x[0]
    seqsList = x[1]
    print(len(x)/2," sequences are now in list")
    ids_outcome = get_ids(id_file)
    ids = ids_outcome[0]
    #print(ids)
    reversing = ids_outcome[1]
    print(len(ids)," IDS are now in list")
    c = write_fasta(outfileName,ids,idList,seqsList)

    print("out of " , str(len(x)/2)," total sequences retained:", c)
    
    
