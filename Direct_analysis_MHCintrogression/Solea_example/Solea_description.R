#This is exemplary R script used to perform all necessary analyses for Solea sp.
#Despite small changes related to the nature of systems, each system is analysed in a similar way. 

library(dplyr)
library(tidyverse)

setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/General_analysis/Replicable_code/Solea_example")
#
### RAW DATA ####
genotype <- read.csv("Solea_genotypes_bothMHC.csv") %>% filter(!grepl("d$",x = INDIVIDUAL_ID)) %>% filter(GENOTYPE == 1) 
#^ recieved by running the Amplisas on raw fastq data and then executing "../../AmpiSAS_resultAnalysis/Amplisas_output_analysis.R"

localities <- read.csv("Solea_localities.csv") %>% mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID))
#^ recieved from collaborators or retreived from literature

genetic_info <- read.csv("Solea_geneticInfo.csv") %>% mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID))
#^ recieved from collaborators or retreived from literature
#
### ALLO - PARA ####
#classifying an individual into allo and parapatry
para_and_symp <- c("Annaba","Lake-Mellach","Tabarka","Bizerte","Tunis2","Tunis","Chergui") #it was decided manually based on the map

dt <- left_join(genotype,localities,by = "INDIVIDUAL_ID") %>% 
  left_join(genetic_info,by = "INDIVIDUAL_ID") %>%
  group_by(LOCALITY) %>% mutate("MEAN_ADMIXTURE" = mean(AEGYPTIACA)/sum(mean(AEGYPTIACA),mean(SENEGALENSIS))) %>% ungroup() %>%
  filter(!is.na(SPECIES) & SPECIES != "mixed" & CLEAR == T & (MEAN_ADMIXTURE > .95 | MEAN_ADMIXTURE < .05)) %>%
  select(INDIVIDUAL_ID,LOCALITY,SPECIES,HYBRID_ZONE) %>% distinct() %>%
  group_by(LOCALITY) %>% 
  mutate("TYPE" = case_when(!(LOCALITY %in% para_and_symp) ~ "allo",T ~ "para")) %>% ungroup() %>%
  select(INDIVIDUAL_ID,HYBRID_ZONE,TYPE) %>% distinct() %>% arrange(TYPE)

write.csv(dt,"non_essential/Solea_alloPara.csv",row.names = F,fileEncoding = "UTF-8")
#
### PERMUTATIONS ####
#creates input file for permutation tests - python script used for analyses can be found in "../../MHC_allele_sharing_with_permutations_TG_SES.py"
allo_para <- read.csv("non_essential/Solea_alloPara.csv",fileEncoding = "UTF-8") %>% mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID))
klasa <- "I"

dt <- left_join(genotype,localities,by = "INDIVIDUAL_ID") %>%
  left_join(select(genetic_info,-HYBRID_ZONE),by = "INDIVIDUAL_ID") %>%
  left_join(allo_para,by = "INDIVIDUAL_ID") %>%
  filter(!is.na(TYPE)) %>%
  mutate("HYBRID_ZONE" = stringr::str_replace_all(HYBRID_ZONE,"_","")) %>%
  mutate("InfoType" = paste(SPECIES,HYBRID_ZONE,TYPE,LOCALITY,sep = "_")) %>%
  filter(CLASS == klasa) %>%
  select(INDIVIDUAL_ID,InfoType,ALLELE,GENOTYPE) %>%
  pivot_wider(names_from = ALLELE,values_from = GENOTYPE,values_fill = 0)

write.table(x = dt,file = paste0("./Permutations/templates/Solea_MHC",klasa,"_permutations_template.tsv"),row.names = F,fileEncoding = "UTF-8",sep = "\t",quote = F)
#
#### ALELLES FREQUENCIES IN SPECIES FOR CLINES ####
allo <- read.csv("non_essential/Solea_alloPara.csv",fileEncoding = "UTF-8") %>% mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID)) %>%
  filter(TYPE == "allo")

df1 <- genotype %>%
  left_join(genetic_info,by = "INDIVIDUAL_ID") %>% select(-HYBRID_ZONE) %>% #add genomic data
  left_join(localities,by = "INDIVIDUAL_ID") %>%
  left_join(allo,by = "INDIVIDUAL_ID") %>%
  filter(SPECIES != 'hybrid' & GENOTYPE == 1 & CLEAR == T & TYPE == "allo") %>%
  select(INDIVIDUAL_ID,CLASS,ALLELE,SPECIES,HYBRID_ZONE) %>% #remove unnecessary columns
  group_by(CLASS,SPECIES,HYBRID_ZONE) %>% mutate("n_ind" = n_distinct(INDIVIDUAL_ID)) %>% ungroup() %>% #get number of individuals per species & CLASS
  group_by(CLASS,SPECIES,ALLELE,HYBRID_ZONE) %>% summarise("n_ind_all" = n_distinct(INDIVIDUAL_ID),"n_ind" = unique(n_ind)) %>% ungroup() %>% #get number of individuals possesing a given alelle
  mutate("ALELLE_FREQ" = n_ind_all/n_ind) %>% #calculate frequencies
  select(-n_ind,-n_ind_all) #remove unnecessary columns

write.csv(df1,"./hzar/Solea_alelle_freqs.csv",row.names = F,fileEncoding = "UTF-8")
#
#### TABLE FOR HZAR  - H INDEX ####
#calculated hybrid index
h_max_likelihood <- function(h,present,s,r){ 
  #H-W frequency of allele (observed values reflect phenotype)
  s_null <- sqrt(1-s)
  r_null <- sqrt(1-r)
  s_full <- 1-(sqrt(1-s))
  r_full <- 1-(sqrt(1-r))
  
  # Observed values reflect actual alleles frequencies
  # s_null <- 1-s
  # r_null <- 1-r
  # r_full <- r
  # s_full <- s
  
  #probabilities
  p1 <- (r_full*h)+((1-h)*s_full)
  p2 <- (r_null*h)+((1-h)*s_null)
  
  # <!> There is a problem here <!>
  # -Inf is returned if a allele is absent in both species or has frequency 1 in both species
  vec <- c()
  for(i in 1:length(present)){
    if(present[i] == 1){
      vec <- c(vec,log(p1[i]^2 + (2*(p1[i]*p2[i]))))
    } else {vec <- c(vec,log(p2[i]^2))}
  }
  vec <- vec[!is.infinite(vec)] #<!> solution
  sum(vec)
} #calculates h-index, it must be run with optimization

alelles_freqs <- read.csv("./hzar/Solea_alelle_freqs.csv")

# <*> calculate h-index for all individual for a given MHC class <*>
klasa <- c("I","II")

#preparation of basic table
df_I <- genotype %>% filter(CLASS == "I") %>%
  select(-CLASS) %>%
  pivot_wider(values_from = "GENOTYPE",names_from = "ALLELE",values_fill = 0)

df_II <- genotype %>% filter(CLASS == "II") %>%
  select(-CLASS) %>%
  pivot_wider(values_from = "GENOTYPE",names_from = "ALLELE",values_fill = 0)

df <- left_join(df_I,df_II,by = "INDIVIDUAL_ID") %>%
  pivot_longer(cols = -1,names_to = "ALLELE",values_to = "GENOTYPE") %>%
  left_join(distinct(select(genotype,ALLELE,CLASS)),by = "ALLELE") %>% na.omit() %>%
  left_join(localities,by = "INDIVIDUAL_ID") %>% 
  left_join(genetic_info,by = "INDIVIDUAL_ID") %>%
  select(INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE,LOCALITY,SPECIES) %>% #add species and locality
  na.omit() %>% filter(SPECIES != "mixed") #remove not needed species

out_list <- list()

hz_list <- list("aegyptiaca" = c("aegyptiaca_senegalensis"),
                "senegalensis" = c("aegyptiaca_senegalensis")) #all hybrid zones for analyzed species

for(l in klasa){
  for(k in 1:length(hz_list)){
    sp <- names(hz_list)[k]
    individuals <- unique(df$INDIVIDUAL_ID[df$SPECIES == sp & df$CLASS == l])
    for(i in 1:length(individuals)){ #looping over all individuals
      ind <- individuals[i] #Single ID
      temp <- filter(df,INDIVIDUAL_ID == ind & CLASS == l) %>% arrange(ALLELE) #table for one individual
      
      for(j in hz_list[[sp]]){ #looping over all hybrid zones an individual belongs to
        
        focal_sp <- str_split(j,"_")[[1]][1] #subtable for the focal species
        temp_focal <- filter(alelles_freqs,SPECIES == focal_sp & HYBRID_ZONE == j) %>% select(ALLELE,ALELLE_FREQ) %>% rename("focal" = "ALELLE_FREQ")
        
        other_sp <- str_split(j,"_")[[1]][2] #subtable for the other species
        temp_other <- filter(alelles_freqs,SPECIES == other_sp & HYBRID_ZONE == j) %>% select(ALLELE,ALELLE_FREQ) %>% rename("other" = "ALELLE_FREQ")
        
        temp_j <- temp %>% #joining info about alleles frequencies in 2 hybridizing species
          left_join(temp_focal,by = "ALLELE") %>%
          left_join(temp_other,by = "ALLELE") %>%
          replace_na(replace = list(focal = 0,other = 0))
        
        #h-index optimization
        opt <- optimize(h_max_likelihood,s = temp_j$focal,r = temp_j$other,present = temp_j$GENOTYPE,interval = c(0,1),maximum = T)
        opt <- round(opt$maximum,3)
        
        outcome <- c(ind,sp,l,j,opt) #collecting data
        out_list[[length(out_list)+1]] <- outcome
      }
    }
  }
}

out_table <- as.data.frame(do.call(rbind,out_list))
colnames(out_table) <- c("INDIVIDUAL_ID","SPECIES","CLASS","HYBRID_ZONE","H_INDEX")
write.csv(out_table,"./hzar/Solea_Hindex.csv",fileEncoding = "UTF-8",row.names = F)
#
### STRUCTURE - PREPARATION ####
#preparing input file for STRUCTURE (discarded alternative to h-index)
klasa <- c("I","II")
hz <- c("aegyptiaca_senegalensis")
dir.create("./Structure/inputs")

for(i in klasa){
  for(m in hz){
    dt <- genotype %>% left_join(localities,by = "INDIVIDUAL_ID") %>%
      left_join(genetic_info,by = "INDIVIDUAL_ID") %>%
      filter(CLASS == i & SPECIES %in% strsplit(m,"_")[[1]]) %>%
      select(INDIVIDUAL_ID,ALLELE,GENOTYPE,LOCALITY) %>%
      left_join(data.frame("LOCALITY" = unique(.[,"LOCALITY"]),"POP" = c(1:length(unique(.[,"LOCALITY"])))),by = "LOCALITY") %>%
      relocate(POP,LOCALITY,.after = 1) %>%
      pivot_wider(names_from = ALLELE,values_from = GENOTYPE,values_fill = 0) %>%
      slice(rep(1:n(),each = 2)) %>% #duplicate rows
      add_row(.before = 1)
    
    dt[1,-c(1:3)] <- 0
    
    out_file_name <- paste("./Structure/inputs/Solea_MHC_",i,"_",m,".tsv",sep = "")
    #out_file <- file(out_file_name,"wb")
    write.table(dt,file = out_file_name,quote = F,sep = "\t",na = "",col.names = F,row.names = F,fileEncoding = "UTF-8")
    #close(out_file)
  }
}
#
### STRUCTURE - COMPILATION ####
# <require> directories with structure's results for several replicates
# <do> binds together results from both MHC classes for a given hybrid zone and make sure Q-axes correspond between classes
# <outcome> single csv file with INDIVIDUAL_ID,Q1,Q2, CLASS & SPECIES
choose_max_likelihood <- function(file_list){
  line_num <- readLines(file_list[1]) %>% grepl(x = .,pattern = "-------------") %>% which() %>% .[5]+1
  lapply(file_list, readLines) %>% sapply(.,"[",line_num) %>% sapply(.,str_split,pattern = "= ") %>%
    sapply(.,"[",2) %>% unname() %>% as.numeric() %>% which.max()
} #structure's replicate with maximum likelihood

read.structure <- function(output_file, corr = T){
  x <- readLines(output_file)
  if(corr == T){
    indexes_start <- which(grepl("-------------",x))[5]+12
  } else {
    indexes_start <- which(grepl("-------------",x))[5]+10
  }
  indexes_stop <- which(grepl("^$",x)) %>% .[. > indexes_start] %>% .[1]-1
  x <- x[indexes_start:indexes_stop]
  dt <- stringr::str_trim(x) %>% str_split(.," +") %>% do.call(rbind,.) %>% as.data.frame() %>%
    select(-V1,-V3,-V4,-V5) %>% rename("INDIVIDUAL_ID" = V2)
  colnames(dt)[-1] <- paste("Q",c(1:(ncol(dt)-1)),sep = "")
  dt <- mutate(dt,Q1 = as.numeric(Q1)) 
  dt
} #read structure file

decide_Qaxis <- function(dt){ #if TRUE keep Q1
  des <- group_by(dt,SPECIES) %>% summarise("Q1" = mean(as.numeric(Q1)),"Q2" = mean(as.numeric(Q2))) %>% 
    arrange(SPECIES)
  if(des$Q1[1] > des$Q2[1]){
    temp <- dt %>% select(-Q1) %>% rename("Q" = "Q2")
  } else {
    temp <- dt %>% select(-Q2) %>% rename("Q" = "Q1")
  }
  temp
}

cor_var <- F #is correlated?

klasa <- c("I","II")
hz <- c("aegyptiaca_senegalensis")

out_list <- list()
for(k in klasa){
  for(z in hz){
    dt <- list.files("./Structure/outputs/",full.names = T) %>% .[grepl(paste0("_",k,"_",z),.)] %>% #MHC-I
      .[choose_max_likelihood(.)] %>%
      read.structure(corr = cor_var) %>%
      mutate("CLASS" = k,"HYBRID_ZONE" = z) %>% 
      left_join(select(genetic_info,INDIVIDUAL_ID,SPECIES),by = "INDIVIDUAL_ID") %>%
      decide_Qaxis(.) %>% select(-SPECIES)
    out_list[[length(out_list)+1]] <- dt
  }
}
dt <- as.data.frame(do.call(rbind,out_list))
write.csv(dt,"./Structure/Solea_structureQ.csv",row.names = F,fileEncoding = "UTF-8")
#
### DISTANCE CALCULATION ALONG TRANSECT ####
transects <- read.csv("non_essential/Solea_TransectLocs.csv",fileEncoding = "UTF-8") #information to which transect a given locality belong (manual)

points_transects <- read.csv("./hzar/distance_calculation/Solea_points_transect.csv",fileEncoding = "UTF-8") %>%
  rename("TRANSECT_LINE" = "TRANSECT") #order of points equally spaced along transects (from map - QGIS) 

dist <- read.csv("./hzar/distance_calculation/Solea_distanceMatrix.csv",fileEncoding = "UTF-8") %>% #distance of localities to points along transects (QGIS)
  rename("DISTANCE_LINE" = "DISTANCE") %>%
  left_join(points_transects,by = c("TargetID" = "FID")) %>%
  mutate("DISTANCE" = DISTANCE/1000) %>% 
  group_by(LOCALITY,TRANSECT_LINE) %>%
  filter(DISTANCE_LINE == min(DISTANCE_LINE)) %>% ungroup() %>%
  left_join(transects,by = "LOCALITY") %>% 
  na.omit() %>% filter(TRANSECT_LINE == TRANSECT) %>%
  select(LOCALITY,DISTANCE,TRANSECT) %>%
  arrange(TRANSECT,DISTANCE)

dist <- dist %>% group_by(TRANSECT) %>%
  mutate("DISTANCE" = DISTANCE-min(DISTANCE))

write.csv(dist,"./hzar/Solea_calculated_distance.csv",row.names = F,fileEncoding = "UTF-8")
#
#### TABLE FOR HZAR ####
#creating an input file for geographic clines (check START_HERE_GeneralAnalyses_3rdPaper.R)
#in folder "./geographic_clines/independent_code/ you will find R script, exemplary input file and output

dist <- read.csv("./hzar/Solea_calculated_distance.csv",fileEncoding = "UTF-8")
#Genome_wide - table
genome_wide <- left_join(genetic_info,localities, by = "INDIVIDUAL_ID") %>%
  left_join(dist,by = "LOCALITY") %>% na.omit() %>%
  group_by(LOCALITY,TRANSECT) %>% summarise("genomic" = mean(SENEGALENSIS)*100,"nSamples_genomic" = n())
#write.csv(dt,"../hzar/Lissotriton_GenomeWide.csv",fileEncoding = "UTF-8",row.names = F)

#H_index - table
h_index <- read.csv("./hzar/Solea_Hindex.csv",fileEncoding = "UTF-8") %>%
  mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID)) %>%
  left_join(localities,by = "INDIVIDUAL_ID") %>%
  left_join(dist,by = "LOCALITY") %>% na.omit() %>%
  group_by(LOCALITY,CLASS,TRANSECT) %>%
  summarise("H_INDEX" = mean(H_INDEX)*100,"nSamples_H" = n()) %>%
  pivot_wider(names_from = CLASS,values_from = c(H_INDEX,nSamples_H))

#structure - table
structure <- read.csv("./Structure/Solea_structureQ.csv",fileEncoding = "UTF-8") %>%
  mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID)) %>%
  left_join(localities,by = "INDIVIDUAL_ID") %>%
  left_join(dist,by = "LOCALITY") %>% na.omit() %>%
  group_by(LOCALITY,CLASS,TRANSECT) %>%
  summarise("Q_INDEX" = mean(Q)*100,"nSamples_Q" = n()) %>%
  pivot_wider(names_from = CLASS,values_from = c(Q_INDEX,nSamples_Q))

#x <- left_join(genotype,localities,by = "INDIVIDUAL_ID")

hzar <- h_index %>%
  full_join(structure,by = c("LOCALITY","TRANSECT")) %>%
  full_join(genome_wide,by = c("LOCALITY","TRANSECT")) %>%
  left_join(dist,by = c("LOCALITY","TRANSECT")) %>%
  rename("HYBRID_ZONE" = "TRANSECT") %>%
  mutate("HYBRID_ZONE" = paste0("aegyptiaca_senegalensis",HYBRID_ZONE)) %>%
  select(HYBRID_ZONE,LOCALITY,DISTANCE,H_INDEX_I,nSamples_H_I,H_INDEX_II,nSamples_H_II,Q_INDEX_I,nSamples_Q_I,
         Q_INDEX_II,nSamples_Q_II,genomic,nSamples_genomic)

write.csv(hzar,"./hzar/Solea_hzar_input.csv",row.names = F,fileEncoding = "UTF-8")
#
### GENOMIC CLINES ####
#creating an input for genomic clines estimation
#you can find the code in "./genomic_clines/independent_code/Genomic_clines_code_parallel.R" together with exemplary input file and output files
dt <- read.csv("./hzar/Solea_Hindex.csv",fileEncoding = "UTF-8") %>% #prepare initial dataset
  mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID)) %>%
  left_join(select(genetic_info,-SPECIES,-CLEAR,-MARKER,-HYBRID_ZONE),by = "INDIVIDUAL_ID") %>%
  separate(HYBRID_ZONE,c("SP1","REFERENCE_SP"),sep = "_",remove = F) %>%
  select(-SP1) %>%
  pivot_longer(cols = c("AEGYPTIACA","SENEGALENSIS"),names_to = "GENETIC_SP",values_to = "GENOME_WIDE") %>%
  filter(REFERENCE_SP == tolower(GENETIC_SP)) %>%
  left_join(select(localities,INDIVIDUAL_ID,LOCALITY),by = "INDIVIDUAL_ID") %>% na.omit()
write.csv(dt,"./GenomicClines/GenomicClines_pointsTab.csv",fileEncoding = "UTF-8",row.names = F)
write.csv(dt,"./GenomicClines/input.csv",fileEncoding = "UTF-8",row.names = F)
#