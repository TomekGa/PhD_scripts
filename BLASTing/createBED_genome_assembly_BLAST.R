# Script to create BED files
# based on tabular BLAST outcome
# used for multiple hits in the single contigs
# takes out all detected sequences

library(dplyr)

query_num <- 3
query_name <- "MHCII_ex3"
flanking_bp <- 50


#### FUNCTIONS ####
# read BLAST tabular for multiple queries
read_blast_tabular <- function(file_v){
  raw <- readLines(file_v)[-length(readLines(file_v))]
  beginning <- grepl("BLASTN",raw)
  indexes <- c(which(beginning == T),length(raw)+1)
  out_list <- list()
  for(i in 1:sum(beginning)){
    tab_range <- raw[(indexes[i]+5):(indexes[i+1]-1)]
    splited <- strsplit(tab_range,split = "\t")
    df <- as.data.frame(do.call(rbind,splited),stringsAsFactors = F) %>% mutate("query" = i)
    colnames(df) <- c("subject_id","subject_title","%_identity","alignment_length","mismatches",
                      "gap_opens","q_start","q_end","s_start","s_end","e_value","bit_score","query")
    out_list[[i]] <- df
  }
  do.call(rbind,out_list)
}

#Create BED file
create_BED <- function(df,bed_file){
  out <- c()
  for(i in 1:nrow(df)){
    s_1 <- as.integer(df$s_start[i])
    s_2 <- as.integer(df$s_end[i])
    if(s_1 < s_2){direction <- "+"}else{direction <- "-"}
    out[i] <- paste(df$subject_id[i],pmin(s_1,s_2)-flanking_bp,pmax(s_1,s_2)+flanking_bp,
                    query_name,0,direction,sep = "\t")
  }
  writeLines(out,bed_file)
}


#READ BLAST tabular 
tabular <- read_blast_tabular("N_tessellata_new_angsd_tabular.txt") %>% filter(query == query_num)
#MAKE BED
BED_name <- paste("Ntesselata_",query_name,".bed",sep = "")
create_BED(tabular,BED_name)

#OPTIONAL - FURTHER STEPS IN BASH ####
#bedtools getfasta -fi <assembly> -bed <BED> -fo <output> -s #get sequences based on bed file
#linsi --adjustdirection --preservecase <in.fasta> > <out.fasta> #align them using mafft
#seqkit seq -m <minimum_length> <in.fasta> > <out.fasta> #filter based on minimum length
