#get a subset of fasta sequences based on names

input_file <- 'P_muralis_cds.fasta' #original fasta file

transcripts <- read.table("MHC_II.bed",sep = "\t") #table with names

out_to_file <- c()
for(i in 1:nrow(transcripts)){
  cmd <- paste("grep -n -e '>",transcripts[i,4],"' ",input_file, sep = "") #ADJUST - column with names
  trs_header <- system(cmd,intern = T)
  start_num <- strsplit(trs_header,split = ":")[[1]][1]
  end_num <- as.character(as.integer(start_num)+50)
  cmd <- paste("sed -n '",start_num,",",end_num,"p' ",input_file,sep = "")
  found_region <- system(cmd,intern = T)
  
  seque <- c()
  for(j in 2:length(found_region)){
    if(!startsWith(found_region[j],">")){
      seque <- c(seque,found_region[j])
    } else {
      break
    }
  }
  out_to_file <- c(out_to_file,found_region[1],seque)
}

writeLines(out_to_file,"transcripts_seqs.fasta")