#Draws circular tree of DNA sequences - filter sequences, align them, build tree, draw tree
#it also clusters MHC sequences based on the predifined expression categories (HEX and LEX)
#it also marks alleles with known STOP codons or frame shifts

#REQUIRES
# sequences in fasta format OR aligned sequences in fasta format OR built tree in newick format
# csv table with ALLELE and expression CATEGORY columns
# csv table with ALLELE and ISSUE (stop or shift)

source("~/Functions_to_be_loaded.R")

library(ggtree)
library(ggplot2)
library(dplyr)
library(seqinr)
library(ape)
library(tidyverse)
library(tidytree)

setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/Triturus/expression_analysis/MHC_I") #ADJUST

#### FUNCTIONS ####
#dt must be a tibble
# get numbers of all nodes that 
# 1. its bootstrap value is higher or equal to threshold
# 2. there is at least one tip from category_v (no others)
# 3. there are only tips from category_v
get_nodes_by_tips <- function(dt,category_v,category_col,threshold){
  tips <- dt$label[!is.na(dt[,category_col]) & dt[,category_col] == category_v]
  out_nodes <- c()
  for(i in unique(dt$parent)){
    tips_plot <- groupClade(dt,i) %>% filter(group == 1)
    if((as.numeric(dt$label[dt$node == i]) >= threshold) | sum(!(tips_plot[[category_col]][tips_plot$isTip == T] %in% category_v)) == 0){
      if((sum(!is.na(tips_plot[,category_col])) > 0) &(sum(pull(tips_plot,category_col) %in% c(category_v,NA)) == nrow(tips_plot))){
        out_nodes <- c(out_nodes,tips_plot$node,i)
      }
    }
  }
  out_nodes
}

#### OPTIONAL - REMOVE SPECIFIC OR DUPLICATED SEQUENCES ####
# there can be some while you simply merge files from several runs
# org <- read.fasta("withFSandSTOP_allMHCII_Triturus_Amb_org.fasta",forceDNAtolower = F) #ADJUST
# org <- org[!(names(org) %in% c("Tri_MHCIex2_1131","Tri_MHCIex2_1129","Tri_MHCIex2_1117"))] #filter specific sequences
# org <- unique(org) #remove duplicates
# nms <- sapply(org, function(x){attr(x,"name")})
# names(org) <- nms
# write.fasta(org,names = names(org),file.out = "withFSandSTOP_allMHCII_Triturus_Amb_why.fasta") #ADJUST

#### OPTIONAL - ALLIGN SEQUENCES ####
#analyzed_file <- "withFSandSTOP_allMHCII_Triturus_Amb_why.fasta" # fasta - tree base
#aligned <- align_mafft(analyzed_file,"Triturus_all_MHCII_Amb_alligned_why.fasta",rev_comp = F) #align - requires locally installed mafft
#aligned <- read.dna("5_Podarcis_MHCI_020821_fifth_om_alligned.fasta",format = "fasta") #backdoor
# if an error occures, the sequences are probably to distinct and some of them must be removed

#### OPTIONAL - BUILD TREE ####
#internally in R
# tree_nm <- paste(strsplit(analyzed_file,split = "\\.")[[1]][1],"_tree.nwk",sep = "") #name of tree
# #tree_org <- build_tree(aligned,tree_nm,multi = T,boot = T) #phangorn, unrooted, NA means 0 bootstrap support (probably)
# tree_org <- build_tree2(aligned,tree_nm,root_v = "Amti-DAB*0801") #ape, rooted, NA means 0 bootstrap support (probably)
# tree_org <- read.tree("withFSandSTOP_allMHCII_Triturus_Amb.nwk") #also backdoor

#MEGA
#1. Use MEGA 7 (Linux version of MEGA X crashes while performing bootstrap)
#2. It required installation of 'libcanberra-gtk-module' and 'libcanberra-gtk3-module' (Linux packages)
#3. Generate instruction file with MegaProto (.mao file)
#4. megacc -a <.mao> -d <aligned.fasta> -o <output_name>
tree_org <- read.tree("withFSandSTOP_allMHCI_Triturus_Novi_tree_mega.nwk") #backdoor
tree_org$node.label <- as.numeric(tree_org$node.label)*100 #bootstrap values into percentage

#FOR BOTH
tree_org$node.label[is.na(as.numeric(tree_org$node.label))] <- 0 #change NAs into zeros
tree_org <- root(tree_org,"Novi") #ADJUST - root tree

#### PREPARE TABLE ####
tree_tab <- as_tibble(tree_org) %>% mutate("isTip" = !(node %in% parent))
#combine with categories
categories <- read.csv("expression_categories_MHCI.csv",stringsAsFactors = F) %>% na.omit() %>% #ADJUST - imports categories
  select(ALLELE,CATEGORY_MEAN)
tree_tab <- left_join(tree_tab,categories,by = c("label" = "ALLELE"))
#combine with frame and stops
fs_stops <- read.csv("bad_seqs_info.csv",stringsAsFactors = F)
tree_tab <- left_join(tree_tab,fs_stops,by = c("label" = "ALLELE"))
tree_tab$CATEGORY_MEAN[tree_tab$label %in% fs_stops$ALLELE] <- "LEX" #treat all sequences with issues as LEX

Hexes <- get_nodes_by_tips(dt = tree_tab,category_v = "HEX",threshold = 90,category_col = "CATEGORY_MEAN") #detect clades of HEX
Lexes <- get_nodes_by_tips(dt = tree_tab,category_v = "LEX",threshold = 90,"CATEGORY_MEAN") #detect clades of LEX

Branch_col <- rep("nic",times = nrow(tree_tab))
Branch_col[tree_tab$node %in% Hexes] <- "HEX"
Branch_col[tree_tab$node %in% Lexes] <- "LEX"

tree_tab <- add_column(tree_tab, Branch_col) #column with branch colors
write.csv(tree_tab,"tree_tab.csv",fileEncoding = "UTF-8",row.names = F)


tree_tab <- read.csv("tree_tab.csv",stringsAsFactors = F) #import saved file to avoid compatibility problems with ggtree package
tree_tab$Branch_col[tree_tab$label == "Novi"] <- "nic"

#### DRAWING TREE ####
analyzed_file <- "withFSandSTOP_allMHCI_Triturus_Novi.fasta" # fasta - tree base
tree_image <- paste(strsplit(analyzed_file,split = "\\.")[[1]][1],"_tree_branches.png",sep = "")

plot <- 
  ggtree(tr = tree_org,layout = "circular",size = 2,aes(col = Branch_col),branch.length = "none") %<+% 
  tree_tab +
  scale_color_manual(values = c("nic" = "black","HEX" = "green","LEX" = "red"))+
  theme(legend.position = "none",plot.margin = unit(c(0,0,0,0),"cm"))+
  geom_tiplab2(aes(label = label,subset = label %in% "Novi"),hjust = 0,size = 15)+ #add label for the root
  #geom_tiplab2(aes(label = tree_tab[,"label"],x =x + 0.05*max(x)),hjust = 0,size = 8,colour = "black")+ #add labels of sequences
  geom_point2(aes(subset = !is.na(ISSUE),x =x + 0.03*max(x)),size = 6,col = "black")+ #black dots for all stop codons and frame shifts
  geom_point2(aes(subset = ((isTip == F) & (as.numeric(tree_tab[,"label"]) >= 90))),size = 4,col = "orange") #orange dots for bootstrap values higher tha n 90%

ggsave(tree_image,plot,device = "png",scale = 1,limitsize = F,height = 45,width = 45)
