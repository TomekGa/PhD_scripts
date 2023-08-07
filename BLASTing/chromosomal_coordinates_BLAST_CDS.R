# Script to get chromosomal coordinates for the hits of BLAST perforem on CDS

#REQUIRES:
# - LINUX environment
# - dplyr {package}
# - original fasta file with CDS sequences (headers should include chromosomal coordinates)
# - pairwise outcome of BLAST search

#OUTCOME:
# - BED file with chromosomal coordinates for BLAST hits 

#arguments
# - 1. [file] original fasta file with CDS
# - 2. [file] pairwise outcome of BLAST
# - 3. [logical] if TRUE, just helps to choose value for arguments 4 & 5
# e.g.

# 1:>ENSPMRT00000000002.1
# cds
# primary_assembly:PodMur_1.0:1:851673:869792:-1
# gene:ENSPMRG00000000002.1
# gene_biotype:protein_coding
# transcript_biotype:protein_coding
# gene_symbol:LRFN5
# description:leucine-rich
# repeat
#   and
# fibronectin
# type-III
# domain-containing
# protein
# 5-like
# [Source:NCBI
#   gene;Acc:114596289]

# - 4. [number] line number including chromosomal coordinates (3 for the example above)
# - 5. [number] line number of including name (1 for the example above)

#RUNNING:
# Rscript chromosomal_coordinates_BLAST_CDS.R CDS.fasta BLAST_pairwise.txt TRUE
# Rscript chromosomal_coordinates_BLAST_CDS.R CDS.fasta BLAST_pairwise.txt FALSE 3 1

args <- commandArgs(trailingOnly = T)
input_file <- args[1]
bed_file <- paste(strsplit(input_file,split = '\\.')[[1]][1],'.bed',sep = '')
index_file <- args[2]
check <- args[3]
headings_file <- system(paste("grep -n -e '^>'",input_file),intern = T)
splited <- strsplit(headings_file[1],split = ' ')[[1]]

#info about genes regions
if(check == T){
  cat(paste(splited,'\n'))
  stop('Choose proper line and rerun with additional arguments')
}

#checking proper division
location_value <- as.integer(args[4])
cat('#### division test ####')
new_split <- strsplit(splited[location_value],split = ':')[[1]][3:6]
cat(paste('\nchrom','chromStart','chromEnd','strand\n',sep = '\t'))
cat(paste(new_split,collapse = '\t'))
name_value <- as.integer(args[5])
new_split <- strsplit(splited[name_value],split = '>')[[1]][2]
cat(paste('\nTranscript name',new_split,sep = '\n'))


#creating BED file
bed_list <- list()
for(i in 1:length(headings_file)){
  splited <- strsplit(headings_file[i],split = ' ')[[1]]
  split_1 <- as.character(strsplit(splited[location_value],split = ':')[[1]][3:5])
  #split_1[1] <- paste('chr',split_1[1],sep='')
  split_2 <- as.character(strsplit(splited[name_value],split = '>')[[1]][2])
  split_3 <- as.integer(strsplit(splited[location_value],split = ':')[[1]][6])
  if(split_3 == 1){
    split_3 <- '+'
  } else {
    split_3 <- '-'
  }
  bin_vec <- c(split_1,split_2,split_3)
  bed_list[[i]] <- bin_vec
}
bed_data <- as.data.frame(do.call(rbind,bed_list))
colnames(bed_data) <- c('chrom','chromStart','chromEnd','name','strand')
#write.table(bed_data,file = bed_file,sep = '\t',col.names = T)

#gene names gathered from BLAST outcome
index_data <- readLines(index_file)
start_seq <- which(startsWith(x = index_data,prefix = 'Sequences producing'))

name <- c()
query_order <- c()
for(i in 1:length(start_seq)){
  checker <- T
  ii <- start_seq[i]+1
  while(checker == T){
    ii <- ii+1
    if(index_data[ii] == ''){
      checker <- F
    } else {
      given_line <- strsplit(index_data[ii],split = " ")[[1]][1]
      name <- c(name,given_line)
      query_order <- c(query_order,i)
    }
  }
}
blast_indexes <- as.data.frame(cbind(name,query_order))
#write.table(blast_indexes,file = "blast_indexes.txt",sep = '\t')

# filtering original table
library(dplyr)
proper_data <- inner_join(x = bed_data,y = blast_indexes,by = 'name')
colnames(proper_data) <- c('chrom','chromStart','chromEnd','name','strand')

new_bed_file <- paste(strsplit(input_file,split = '\\.')[[1]][1],'_filtered.bed',sep = '')
write.table(proper_data,file = new_bed_file,sep = '\t',quote = F,row.names = F)

#### OPTIONAL FURTHER STEPS ######
#Checking coverage

## mapped bam file needs to be converted into bed file
# bedtools bamtobed -i <bam file> > <output file>
## both files should be sorted
# sort -k1,1 -k2,2n <input bed file> > <output bed file>
## coverage checking
# bedtools coverage -a <transcriptomes bed file> -b <mapped bed file> -sorted > <outcome>
#OR
#samtools mpileup -u <bam file> -I <transcriptomes bed file> -f <assembly.fasta> -o <outcome.txt>