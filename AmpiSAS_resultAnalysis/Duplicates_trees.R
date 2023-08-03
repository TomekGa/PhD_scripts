#comparing MHC genotypes of individuals genotyped independently in two runs

#REQUIRED
# ...modified.csv and ...appendix.csv files generated from independent execution of "Amplisas_output_analysis.R" for each run

#OUTPUT
# trees depicting repeatability
# table with the estimated discrepancy 

library(tidyverse)
library(phangorn)
library(ggtree)

source("~/Functions_to_be_loaded.R")

dt <- c("5_Emys_MHCI_080321_fifth_om_modified.csv","5z_Emys_MHCI_080321_fifth_nm_modified.csv")
appx <- c("5_Emys_MHCI_080321_fifth_om_appendix.csv","5z_Emys_MHCI_080321_fifth_nm_appendix.csv")

################################################
# TOTAL - all alleles found in examined runs
# RESEQ - alleles found only in resequenced individuals
################################################

total <- read_and_change(dt,appx,dup_symbols = c("d","f")) #read all files into a single data.frame
# ^ dup_symbols - letters added at the end of ID which indicates replicated individuals within a single run

resequenced <- take_out_resequenced(total) #filter out resequenced individuals

re_info_reseq <- reseq_info(resequenced) #creates the table with the information when was a given allele resequenced 
re_info_total <- reseq_info(total)

dupDR <- discrepancy(resequenced,"discrepancy.csv") #calculates the rate of discrepancy

unique_seqs_reseq <- distinct(resequenced,NAME,.keep_all = T)
unique_seqs_total <- distinct(total,NAME,.keep_all = T)

data_reseq <- save_to_fasta(unique_seqs_reseq,"sequences_dupcheck_reseq.fasta")
data_total <- save_to_fasta(unique_seqs_total,"sequences_dupcheck_total.fasta")

alligned_reseq <- allign_mafft(data_reseq,"allignment_dupcheck_reseq.fasta")
alligned_total <- allign_mafft(data_total,"allignment_dupcheck_total.fasta")
#alligned_reseq <- ape::read.FASTA("allignment_dupcheck_reseq.fasta") ### backdoor
#alligned_total <- ape::read.FASTA("allignment_dupcheck_total.fasta") ### backdoor

tree_reseq <- build_tree(alligned_reseq,"duplicates_tree_reseq.new",multi = F,boot = F)
tree_total <- build_tree(alligned_total,"duplicates_tree_total.new",multi = F,boot = F)

duplicate_gg_tree(tree_reseq,resequenced,re_info_reseq,"duplicates_tree_reseq.png",
                  x_lim = 0.8,boot = F) #draw a tree showing number of replicates and PAF frequency in each
duplicate_gg_tree(tree_total,total,re_info_total,"duplicates_tree_total.png",
                  x_lim = 0.8,boot = F,height_im = NA)
#^ x_lim - modify depending on relative distance between MHC alleles