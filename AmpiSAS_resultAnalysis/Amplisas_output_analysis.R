#### AMPLISAS OUTPUT ANALYSIS
#### This script is semi-automatic and each step must be adjusted to the specificity of a system
#### Tomek Gaczorek
#### tomek.gaczorek@gmail.com

# input: AmpliSAS output - excel file
# outputs: 
# - binary matrix with genotypes
# - csv files with filtered Amplisas output (number of reads and coverage separately)
# - log file including descriptive statistics
# - fasta files: original, filtered and alligned alleles (both nucleotide and proteins)
# - tree of detected alleles (newick and png)
# - estimates of genotyping repeatability

#packages - additionally ape, msa and stringr (not loaded generally as they cover some other functions)
library(ggtree)
library(ggplot2)
library(dplyr)
library(seqinr)
library(phangorn)

source("/home/tomek/Dropbox/MHC_projekt/Shared_MHC_project/Amplisas/R_code/Functions_to_be_loaded.R")
# ^ ADJUST!!! - path to the script with created functions de novo
setwd("~/Dropbox/MHC_projekt/Shared_MHC_project/Podarcis/Amplisas_results/MHC_I/Podarcis_MHCI_171120_first_nm")
# ^ ADJUST!!! - folder with the general result from the amplisas
analysed_file <- "2_second_Natrix_MHCII_080423.xlsx"
# ^ ADJUST!!! - the general result from the amplisas

## BODY #### - GO STEP BY STEP ####

#create log file
log_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_log.txt",sep = "")
writeLines(c(analysed_file,""),log_name)

#reading
data_list <- read_and_adjust(analysed_file) #old sequencing method, assuming the idividual IDs are proper
#^it also assumes the existance of 'MEAN_FREQ', 'MAX_FREQ' and 'MIN_FREQ' collumns (6-8 collumn) - if absent add them and fulfill with zeros

data <- data_list[[1]] # actual data
#View(data) #check if everything is ok
appendix <- data_list[[2]] # table with information about coverage

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-6,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,7:ncol(data)],2,function(x){sum(!is.na(x))}))),2),
                                        "\nMean coverage per individual: ",round(mean(as.numeric(appendix[1,2:ncol(appendix)])),2)),log_name)

#sort by length
data <- arrange(data,LENGTH)
most_frequent(data$LENGTH) #printing most frequent value
#View(data) #check the proper length, usually the most frequent one
desired_length <- most_frequent(data$LENGTH)
#desired_length <- 224

#filtering by length
data <- filter_by_length(data,desired_length,2,rem_shifts = F,log_name) 
# ^ digit refers to the acceptable deviation from the proper length (number of codons)
# ^ rem_shifts - do you want to remove frame shifts at this step?

#filtering by amplicon depth (total coverage of an individual)
quantile(as.numeric(appendix[1,2:ncol(appendix)]),probs = seq(0,1,0.05)) #distribution of total coverage
#View(appendix)
data <- filter_by_amplicon_depth(data,appendix,log_name,minimum_depth = 1000)
# ^ usually more than 300; discarded should not constitute more than 5 %

#limit data to species - OPTIONAL and rarely used
#desired_species <- c("cristatus","dobrogicus","macedonicus",'ivanbureschi',"anatolicus")
#data <- limit_to_species(data,"../../Triturus_ids.csv",desired_species,log_name)
# ^ note you need to have a file including species and individual ID - look for the original function to adjust to your file

#check if there are some reverse complements
#data <- check_reversal(data,analysed_file) #without reference
data <- check_reversal_alt(data,analysed_file) #it uses reference file with alleles gathered from previous runs in a system

# contamination checking (delation can change the PAFs) - OPTIONAL
#contamination_checking(data,"../../Triturus_ids.csv",dup_symbols = "d")
# ^it generates html table which might help to see some indications of contamination
# ^note you would need additional file with species and individual IDs

#additionally removed individuals (OPTIONAL, you must have good reason to delete individual e.g. contamination)
# data <- remove_inds(data,c(14791,14627,13157,14624,14663,14784,"14784d"),log_name,
#                     mess = "INDIVIDUALS REMOVED BASED ON TOO HIGH ALLELES NUMBER OR SEQUENCING PROBLEMS")

# adding alleles frequencies - they are recalculated because of filtering before
list_of_things <- adding_freqs_info(data,appendix)
data <- list_of_things[[1]]
freqs_data <- list_of_things[[2]]
write.csv(freqs_data,"freqs_data.csv",row.names = F,fileEncoding = "UTF-8")

# plotting max frequencies - rough way to infer PAF threshold
#plot_name <-  paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_freqPlot.png",sep = "")
#plot_max_freq(data,plot_name,w = 40,t = .1) # w refers to width of picture, the highest the better resolution,
# ^ t - threshold of plotted points (lower than)
data <- filter_by_tag_switch(data,freqs_data,threshold = .25,log_name) #fitering by PAF threshold
#^we often tried with multiple PAF thresholds and compared results afterwards

#modifying names
#system abbreviation + marker + number
#subsequent numbers are provided
#if order > 1, numbering takes into account previous runs (order refer to run)
data <- modify_names_new(data,abr = "Nat",gene = "MHCIIex2",order = 2)

#are there any Ns in sequences?
Are_Ns(data)

#optional - remove some sequences (e.g. if you found Ns)
data <- remove_some_seqs(data,alleles = c("Nat_MHCIIex2_499","Nat_MHCIIex2_240","Nat_MHCIIex2_358","Nat_MHCIIex2_536"),log = log_name,pat_or_vec = "VEC")

#writing filtered sequences to fasta
fasta_name <- save_to_fasta(data,paste("withFrameShifts",analysed_file,sep = "_"))

# OPTIONAL - writing modified file (still with frame shifts - unless filtered above)
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified_FS.csv",sep = "")
write.csv(data,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)
#writing appendix
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_appendix_FS.csv",sep = "")
write.csv(appendix,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

#OPTIONAL filtering based on frame shifts
frame_shifts <- detect_frame_shifts(data,desired_length) #checking frame shifts based on length
data <- remove_some_seqs(data,alleles = frame_shifts,log = log_name,pat_or_vec = "VEC")

#repeatability (OPTIONAL, you must have individuals sequenced twice;
#by convention we add 'd' at the end of ID to distinguish duplicated one)
repeatability(data,log_name)

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-9,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,10:ncol(data)],2,function(x){sum(!is.na(x))}))),2)),log_name)

#writing modified file (after whole filtering, including frame shifts)
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified.csv",sep = "")
write.csv(data,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)
#writing appendix
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_appendix.csv",sep = "")
write.csv(appendix,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

#writing sequences to fasta
fasta_name <- save_to_fasta(data,analysed_file)

#fasta_bin <- seqinr::read.fasta("4_Emys_MHCII_160521_forth_nm.fasta") ### backdoor

allignment_name_file <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned.fasta",sep = "")
#alligned <- allign(fasta_bin,allignment_name_file) #alligned sequences with ClustalW (from msa package)
alligned <- align_mafft(fasta_name,allignment_name_file,rev_comp = F) #allign with locally installed mafft - check function to adjust
# ^path_mafft argument to provide path to installed mafft
# ^built to work on Linux but you can easily make your own alignment and just provide it to the next step (just adjust "allignment_name_file" variable) 
# ^ MUST BE CHECKED MANUALLY if alignment is correct

fasta_bin <- seqinr::read.fasta(allignment_name_file,forceDNAtolower = F)
frame_translate <- 2 #ADJUST proper frame (from 0 to 2)
with_stop_codons <- detect_stop_codons(fasta_bin,log_name,frame_val = frame_translate,stop_verify = F)

frame_shifts <- frame_shifts[!(frame_shifts %in% with_stop_codons)]
# ^if you decided to keep frame shifts from this point all sequences with stop codons are not treated as frame shifts

#building tree
#alligned <- ape::read.FASTA("1_Lissotriton_MHCI_first_080321_om_alligned.fasta") ### backdoor
tree_nm <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree.new",sep = "")
# not sure how it would behave under Windows as the argument mc_cores within bootstrap.pml function 
# (look into build_tree()) should be used only under UNIX based systems. 
# Without this limit my computer was stacked. In case of it maybe turn off 
# parralelization and just wait longer (multi argument)
alligned <- ape::read.FASTA(allignment_name_file)
tree <- build_tree(alligned,tree_nm,multi = T,boot = F) # the heaviest part of the code but if you do not need bootstrap it is super quick

#showing proper tree
tree_image <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree.png",sep = "")
create_gg_tree(tree,with_stop_codons,frame_shi = frame_shifts,
               image_file = tree_image,x_lim = 1,height_im = 25)
# ^ x_lim - size of the x axis - often needs to be adjusted depending on the branches length
# ^ height_im - height of the imag  e - needs to be increased when there is many alleles

### remove non-classical
stops <- with_stop_codons # seqs with stop codons

#to_remove <- c("Tri_MHCIIex2_276","Tri_MHCIIex2_271") # this if additional alleles are removed based on tree
# ^ it can happen that e.g. one allele is clustered with stop codons, you can consider removing it
to_remove <- c() #if nothing additional to remove

#one more data modification
data2 <- remove_some_seqs(data,alleles = c(stops,to_remove,frame_shifts),log = log_name,pat_or_vec = "VEC") #removing stop codons, frame shifts and additional choses alleles
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified_withoutSTOPS.csv",sep = "")
write.csv(data2,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

repeatability(data2,log_name) #final assessment of genotyping repeatability

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-9,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,10:ncol(data)],2,function(x){sum(!is.na(x))}))),2)),log_name)

removed <- alligned[!(labels(alligned) %in% c(stops,to_remove,frame_shifts))] #filtering allinged sequences 

#modified tree
# tree_nm <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree_wihoutSTOPS.new",sep = "")
# tree <- build_tree(removed,tree_nm,multi = T,boot = F) 
# tree_image <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree_withoutSTOPS.png",sep = "")
# create_gg_tree(tree,with_stop_codons,frame_shi = frame_shifts,
#                image_file = tree_image,x_lim = 0.3,height_im = 14)

allignment_func <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_functional.fasta",sep = "")
write.FASTA(x = removed,file = allignment_func) #final alignment of functional alleles

#translate to protein
allignment_aa <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_aa.fasta",sep = "")
protein_alignment(unalligned_DNA = fasta_bin,not_needed = c(),
                  framing = frame_translate,out_file = allignment_aa)
allignment_func_aa <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_functional_aa.fasta",sep = "")
protein_alignment(unalligned_DNA = fasta_bin,not_needed = c(stops,to_remove,frame_shifts),
                  framing = frame_translate,out_file = allignment_func_aa) #protein alignment

#change to matrix with genotypes
#data <- read.csv("Emys_MHCI_ex2_combined_modified.csv",header = T,check.names = F) ### backdoor
gen_file <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_genotypes.csv",sep = "")
genotypes <- change_to_genotypes(data,removed,gen_file) #table with genotypes

