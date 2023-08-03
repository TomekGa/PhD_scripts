#Analysis of alleles which are similar and might be the effect of genotyping error
#This script was used to determine the proper clustering threshold of AmpliSAS for Emys system.

library(dplyr)

#REQUIRES
#...modified.csv file which is an output of Amplisas_output_analysis.R script

#OUTPUTS
# table showing the results of clustering
# table with combined clustering and original ...modified.csv table
# table with clusters and individuals suspected of posessing alleles resulting from genotyping errors

dt <- read.csv("1_Emys_MHCII_160521_first_nm_modified.csv",stringsAsFactors = F,check.names = F) #ADJUST
dt <- arrange(dt,desc(MAX_FREQ)) %>% filter(LENGTH == 197) #as the clusters are based on direct comparison of nucleotides, all sequences should have equal length

seqs <- dt$SEQUENCE

compare <- function(x,ref){
  x <- strsplit(x,split = "")[[1]]
  ref <- strsplit(ref,split = "")[[1]]
  score <- sum(x == ref)
  differ <- length(ref)-score
  if(differ > 1){
    F
  } else {T}
} # find if seqs differ by more than 1 bp

#clustering sequences which differ by at most 1 bp
out_list <- list(seqs[1])
for(i in 2:length(seqs)){
  cycles <- 0
  for(j in 1:length(out_list)){
    compared <- compare(seqs[i],out_list[[j]][1])
    if(compared == F){
      cycles <- cycles + 1
      if(cycles == length(out_list)){
        out_list[[length(out_list)+1]] <- seqs[i]
      }
    } else {
      out_list[[j]] <- c(out_list[[j]],seqs[i])
      break
    }
  }
}

#make df showing the number of sequences in each cluster
df_list <- list()
for(i in 1:length(out_list)){
  temp <- data.frame("SEQUENCE" = out_list[[i]],"CLUSTER" = length(out_list[[i]]),"ID" = i)
  df_list[[length(df_list)+1]] <- temp
}
df_cluster <- do.call(rbind,df_list)
df_cluster <- arrange(df_cluster,desc(CLUSTER),ID)

write.csv(df_cluster,file = "clusters_seqs.csv",fileEncoding = "UTF-8")

#joining with the original data
dt_c <- left_join(df_cluster,dt,by = "SEQUENCE")
write.csv(dt_c,file = "clusters_table.csv",fileEncoding = "UTF-8")

#identifying clusters with more than 1 sequences and showing individuals posessing less frequent variants in a cluster
clusters <- unique(dt_c$ID[dt_c$CLUSTER > 1])
clust_list <- list()
for(i in clusters){
  used <- filter(dt_c,ID == i) %>% arrange(desc(MAX_FREQ))
  blank_cols_logit <- unname(apply(apply(used,2,is.na),2,sum) == nrow(used))
  used <- used[,!blank_cols_logit]
  #some_dirt <- sum((apply(apply(used,2,is.na),2,sum) == (nrow(used)-1)) & (!is.na(used[1,])))-11
  without_b <- sum(is.na(used[1,]))
  odds <- paste(colnames(used)[is.na(used[1,])],collapse = " ")
  out <- c("ID" = i,"SIZE" = unique(dt_c$CLUSTER[dt_c$ID == i]),"ALL_INDS" = ncol(used)-11,
           "WITHOUT_BEST" = without_b,"ODDS" = odds)
  clust_list[[length(clust_list)+1]] <- out
}
clust_df <- do.call(rbind,clust_list)
write.csv(clust_df,file = "clusters_odds.csv",fileEncoding = "UTF-8")















