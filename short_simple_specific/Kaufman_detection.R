#script to check for unconserved MHC amino-acids based on Kaufman 1994
#it helps to determine non-classical MHC alleles

#REQUIRES:
# fasta file with aligned MHC protein sequences
# known position of putatively conserved aa based on Kaufman 1994

#OUTPUT:
#table indicating the number of unexpected aa for each allele

library(seqinr)

### CHANGE HERE ####
infile <- "5_Podarcis_MHCII_020821_fifth_om_alligned_aa.fasta" #alignment fasta
position <- c(52,72,73) #putatively conserved bp positions in your alignment
class_v <- "II" #MHC class (function has build-in Kaufmann's aa for exons2 in both MHC classes)
##############

alleles_badKaufman <- function(infile,positions,class_v){
  require(seqinr)
  require(dplyr)
  
  if(class_v == "I"){
    desired <- list("order1" = c("Y","S"),"order2" = c("Y","L","H"),"order3" = c("Y","R","A","H"))
  } else if(class_v == "II"){
    desired <- list("order1" = c("W","V"),"order2" = c("Y","T","H"),"order3" = c("N","H"))
  }
  
  fas <- read.fasta(infile,forceDNAtolower = F) %>% do.call(what = rbind)
  
  outcome <- list()
  for(i in 1:length(position)){
    temp <- fas[,position[i]]
    indexes <- which(!(temp %in% desired[[i]]))
    outcome[[i]] <- names(temp)[indexes]
  }
  outcome
  
  all_alleles_unlinked <- unique(unlist(p))
  bad_seqs <- c()
  for(i in all_alleles_unlinked){
    comparisons <- c()
    for(j in 1:length(outcome)){
      comparisons[length(comparisons)+1] <- i %in% outcome[[j]]
    }
    bad_seqs[length(bad_seqs)+1] <- sum(comparisons)
  }
  as.data.frame(cbind(all_alleles_unlinked,bad_seqs))
}

p <- alleles_badKaufman(infile,position,class_v)
write.csv(p,"5_Podarcis_MHCII_badBasedKaufman.csv",row.names = F,fileEncoding = "UTF-8")