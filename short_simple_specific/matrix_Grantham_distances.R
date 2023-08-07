#calculates matrix of Grantham distances between protein sequences

library(ape)
library(seqinr)

#reads fasta nt alignment of full codons into DNAbin
nt_alg <- read.FASTA("Example_alignment_full_codons.fas")

#translates into protein in AAbin format
#IMPORTANT - gaps are encoded as "X", not "-"
#you can use this object for calculation of distance
prot_algn <- trans(nt_alg)
#or you can convert it to alignment format used by seqinr
prot_algn_seqinr <- ape::as.alignment(as.character(prot_algn))

#produces object of class dist I want
prot_dist <- dist.gene(as.matrix(prot_algn), method = "pairwise", pairwise.deletion = TRUE)

#g <- read.FASTA("Example_alignment_full_codons_protein.fas",type = "AA")
#g <- as.character(g)

#problematic <- readRDS("AAbin_do_Granthama.rds")

#### funkcja ####

calculate_grantham <- function(protein_allignment){
  # input - object of class AAbin
  # grantham table
  gram_table <- read.csv("https://dl.dropbox.com/s/0awueqeyodsvveu/Gratham_Distance_table.csv?dl=1",
                         header = T, row.names = 1,encoding = "UTF-8")
  #outcome matrix
  out_mat <- matrix(0,nrow = length(protein_allignment),ncol = length(protein_allignment))
  pa_names <- names(protein_allignment)
  colnames(out_mat) <- pa_names
  row.names(out_mat) <- pa_names
  
  # distance between 2 sequences
  gram_dist <- function(seq1,seq2,gr_tab){
    out_dist <- c()
    for(k in 1:length(seq1)){
      one <- seq1[k]
      two <- seq2[k]
      if((one == "X" || two == "X")){
        out_dist <- c(out_dist,NA)
      }
      else if(one == "*" || two == "*"){
        next
      }
      else if(identical(one,two) == T){
        out_dist <- c(out_dist,0)
      }
      else {
        if(two == "S" || is.na(gr_tab[one,two])==T){
          out_dist <- c(out_dist,gr_tab[two,one]) 
        } else {out_dist <- c(out_dist,gr_tab[one,two])}
      }
    }
    out_dist <- na.omit(out_dist)
    sum(out_dist)/length(out_dist)
  }
  
  # calculation of distance for each pair
  for(i in c(1:length(protein_allignment))){
    for(j in c(1:length(protein_allignment))){
      if(pa_names[i] != pa_names[j]){
        out_mat[i,j] <- gram_dist(toupper(as.character(protein_allignment[[i]])),
                                  toupper(as.character(protein_allignment[[j]])),
                                  gram_table)
      }
    }
  }
  as.dist(out_mat)
}

#checking
x <- calculate_grantham(prot_algn)
