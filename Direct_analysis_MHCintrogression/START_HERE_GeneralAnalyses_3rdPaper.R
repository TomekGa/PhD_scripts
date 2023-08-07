#this script assumes you have all data which must be calculated independently for each system
#all exemplary data are provided and you can easily follow this code
#However, if you wish to see how it was done, follow Solea example in "Solea_example" folder

library(dplyr)
library(tidyverse)

setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/General_analysis/Replicable_code/")
#
### COPYING ALL RELEVANT FILES ####
# taxs <- paste0("../../",taxons,"/GenomicClines/Genomic_clines_plotting_wIntervals.csv")
# outs <- paste0("./genomic_clines/",taxons,"_Genomic_clines_plotting_wIntervals.csv")
# sapply(1:length(taxs),FUN = function(x){file.copy(from = taxs[x],to = outs[x])})
# #
### DESCRIPTIVE STATS ####

#functions
read_and_rbind_gen <- function(taxons){
  #taxs <- paste0("../",taxons,"/Exit_files/",taxons,"_genotypes_bothMHC.csv")
  taxs <- paste0("./genotypes/",taxons,"_genotypes_bothMHC.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") %>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

read_and_rbind_GI <- function(taxons){
  #taxs <- paste0("../",taxons,"/Exit_files/",taxons,"_geneticInfo.csv")
  taxs <- paste0("./genetic_info/",taxons,"_geneticInfo.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") %>% lapply(.,select,INDIVIDUAL_ID,SPECIES,CLEAR,MARKER)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

read_and_rbind_loc <- function(taxons){
  taxs <- paste0("./localities/",taxons,"_localities.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") %>% lapply(.,select,INDIVIDUAL_ID,LOCALITY)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

#vector of taxa
taxons <- c("Salmo","Solea","Bombina","Lissotriton","Triturus","Anguis","Emys","Natrix","Podarcis","Sphyrapicus","Myodes")

#files import
genotypes <- read_and_rbind_gen(taxons) %>% filter(GENOTYPE == 1 & !grepl("d$",INDIVIDUAL_ID))
genetic_info <- read_and_rbind_GI(taxons)
localities <- read_and_rbind_loc(taxons)

# number of inds & localities per species (>50% ancestry)
temp <- left_join(genotypes,localities,by = c("INDIVIDUAL_ID","TAXON")) %>%
  left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>%
  select(INDIVIDUAL_ID,LOCALITY,SPECIES,TAXON) %>% distinct() %>%
  na.omit() %>% group_by(TAXON,SPECIES) %>% summarise("N_IND" = n_distinct(INDIVIDUAL_ID),"N_LOC" = n_distinct(LOCALITY))

# number of genotyped individuals for a given genome-wide marker
temp <- genetic_info %>% filter(INDIVIDUAL_ID %in% genotypes$INDIVIDUAL_ID) %>% count(TAXON,MARKER)

# % of admixed
temp <- genetic_info %>% filter(INDIVIDUAL_ID %in% genotypes$INDIVIDUAL_ID) %>% count(TAXON,CLEAR) %>%
  group_by(TAXON) %>% mutate("SUM_N" = sum(n)) %>% filter(CLEAR == F) %>% mutate("%" = n/SUM_N*100)

#number of alleles per species + plot
taxons <- c("Salmo","Solea","Bom.","Lis.","Triturus","Anguis","Emys","Natrix", "Podarcis","Sph.","Myo.") #abbreviations for display

temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  group_by(TAXON,SPECIES,INDIVIDUAL_ID,CLASS) %>% summarise("N_ALL" = n_distinct(ALLELE)) %>%
  mutate("TAXON" = case_when(TAXON == "Bombina" ~ "Bom.",TAXON == "Lissotriton" ~ "Lis.",TAXON == "Sphyrapicus" ~ "Sph.",
                             TAXON == "Myodes" ~"Myo.",T ~ TAXON)) %>%
  mutate("CLASS" = paste("MHC",CLASS),"TAXON" = factor(TAXON,levels = taxons))

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

p <- ggplot(temp,aes(x = SPECIES,y = N_ALL),linewidth = 1.5) + 
  #geom_boxplot(aes(x = SPECIES,y = N_ALL))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",linewidth = .8,width = 0.8)+
  facet_grid(CLASS~TAXON,scales = "free_x",space = "free")+
  labs(y = "Number of alleles per individual")+
  theme_bw()+
  #facet_wrap(vars(CLASS))+
  theme(axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1,face = "italic"),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        strip.text.x = element_text(size = 25,face = "italic"),strip.text.y = element_text(size = 25),axis.text = element_text(size = 25),axis.title.y = element_text(size = 30),
        strip.background = element_rect(fill = "white"))

ggsave("plots/N_alleles_ref.png",plot = p,device = "png",width = 20,height = 12)

#n_alleles species - summary
temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid") %>%
  filter(CLEAR == T) %>%
  group_by(TAXON,SPECIES,CLASS) %>% summarise("N_ALL" = n_distinct(ALLELE)) %>%
  arrange(TAXON,SPECIES,CLASS)

#n_alleles taxa - summary
temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  group_by(TAXON,CLASS) %>% summarise("N_ALL" = n_distinct(ALLELE)) %>%
  arrange(TAXON,CLASS)

#n_alleles per ind and species - summary
temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  group_by(TAXON,CLASS,SPECIES,INDIVIDUAL_ID) %>% mutate("N_ALL" = n_distinct(ALLELE)) %>% ungroup() %>%
  group_by(TAXON,CLASS,SPECIES) %>% summarise("MEAN_ALL" = round(mean(N_ALL),2)) %>% ungroup() %>%
  arrange(TAXON,CLASS)

#n_alleles per ind and taxa - summary
temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  group_by(TAXON,CLASS,SPECIES,INDIVIDUAL_ID) %>% mutate("N_ALL" = n_distinct(ALLELE)) %>% ungroup() %>%
  group_by(TAXON,CLASS) %>% summarise("MEAN_ALL" = round(mean(N_ALL),2)) %>% ungroup() %>%
  arrange(TAXON,CLASS)

#allelic richness
allelic_richness <- function(x,size_v = 4,cycles = 1000){
  out <- c()
  for(i in 1:cycles){
    inds <- sample(unique(x$INDIVIDUAL_ID),size_v,replace = F)
    out[length(out)+1] <- x %>% filter(INDIVIDUAL_ID %in% inds & GENOTYPE == 1) %>%
      summarise("Allelic_richness" = n_distinct(ALLELE)) %>% .$Allelic_richness
  }
  mean(out)
}

temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  nest_by(TAXON,SPECIES,CLASS) %>% group_by(TAXON,SPECIES,CLASS) %>%
  mutate("ALL_RICHNESS" = unlist(purrr::map(data,allelic_richness))) %>% select(-data)
write.csv(temp,"temp_tables/AllelicRichness_perSpecies.csv",row.names = F,fileEncoding = "UTF-8")

temp <- genotypes %>% left_join(genetic_info,by = c("INDIVIDUAL_ID","TAXON")) %>% na.omit() %>%
  filter(SPECIES != "mixed" & SPECIES != "hybrid" & SPECIES != "thyroideus") %>%
  filter(CLEAR == T) %>%
  nest_by(TAXON,CLASS) %>% group_by(TAXON,CLASS) %>%
  mutate("ALL_RICHNESS" = unlist(purrr::map(data,allelic_richness))) %>% select(-data)
write.csv(temp,"temp_tables/AllelicRichness_perTaxon.csv",row.names = F,fileEncoding = "UTF-8")
#
### DRAW SPECIES TREE ####
library(ape)
library(ggtree)
library(ggplot2)

tree_org <- read.tree("permutations/MHC_introgression_timetree.nwk")
tree_org$tip.label <- stringr::str_replace_all(tree_org$tip.label,"_"," ")

nms <- c("Myodes rutilus","Myodes glareolus",
         "Sphyrapicus nuchalis","Sphyrapicus ruber","Sphyrapicus varius",
         "Emys trinacris","Emys orbicularis",
         "Podarcis guadarramae","Podarcis bocagei","Podarcis virescens",     
         "Podarcis vaucheri","Podarcis hispanicus","Podarcis carbonelli",    
         "Podarcis liolepis",
         "Anguis fragilis","Anguis colchica",
         "Natrix natrix","Natrix helvetica","Natrix astreptophora",
         "Triturus cristatus","Triturus macedonicus",  
         "Triturus anatolicus","Triturus ivanbureschi","Triturus pygmaeus",      
         "Triturus marmoratus","Lissotriton vulgaris","Lissotriton montandoni", 
         "Bombina variegata","Bombina bombina","Solea aegyptiaca",       
         "Solea senegalensis","Salmo trutta","Salmo salar")

tree_org <- ape::rotateConstr(tree_org,nms) #changes the order of tips

### PLOTTING ###
p <- ggtree(tr = tree_org,ladderize = F,size = 5)+ #it draws tree and moves leaves more space around x-axis but the axis is too long
  scale_x_continuous(breaks = seq(0,450,50),limits = c(0,600),expand = c(0,0))+
  geom_tiplab(size = 20)+
  theme_tree2(plot.margin = margin(50, 50, 50, 50))+
  ylim(0,33)+
  theme(axis.line.x = element_line(linewidth = 5),axis.text.x.bottom = element_text(size = 50,vjust = .1))

p <- ggtree(tr = tree_org,ladderize = F,size = 5)+ #this is better but there is no much space around x-axis
  geom_tiplab(size = 20,as_ylab = T)+
  theme(plot.margin = margin(50, 50, 50, 50))+
  theme(axis.line.x = element_line(linewidth = 5),axis.text.x.bottom = element_text(size = 80,vjust = .1),
        axis.text.y.right = element_text(size = 80,face = "italic"),axis.title.x.bottom = element_text(size = 80,margin = margin(t = 30)))

p <- revts(p)+ #reverts time
  scale_x_continuous(breaks = seq(-450,0,50),expand = c(0,0),labels = abs(seq(-450,0,50)),name = "Time [My]")

ggsave("plots/All_species_Tree.png",p,device = "png",scale = 1,limitsize = F,height = 40,width = 60)
#
### PLOTS - PERMUTATIONS ####
taxons <- c("Solea","Bombina","Lissotriton","Triturus","Anguis","Emys","Natrix", "Podarcis","Sphyrapicus","Myodes")
# ^ Salmo was excluded because of lack of parapatric individuals and significant sharing in sympatry

read_and_rbind_perm <- function(taxons){
  taxs <- paste0("./permutations/",taxons,"_permutation_output_manual.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

# DRAW PERMUTATION PLOT
taxons1 <- c("Sol.","Bom.","Lis.","Triturus","Ang.","Emys","Natrix", "Podarcis","Sphyrapicus","Myo.")

dt <- read_and_rbind_perm(taxons) %>%
  filter(VARIANT == "all") %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = " - ")) %>%
  pivot_longer(cols = c("FAR","CLOSE"),names_to = "WHERE",values_to = "TYPE") %>%
  pivot_longer(cols = c("PERC_FAR","PERC_CLOSE"),names_to = "PERC",values_to = "PERC_VALUE") %>%
  mutate("KEEP" = case_when(str_detect(PERC,WHERE) ~ T, T ~ F)) %>%
  filter(KEEP == T) %>%
  mutate("CLASS" = paste("MHC",CLASS)) %>%
  mutate("TYPE" = factor(TYPE,levels = c("para","allo"))) %>%
  mutate("PERC_VALUE" = case_when(PERC_VALUE == 0 ~ 0.1,T ~ PERC_VALUE)) %>%
  select(-WHERE,-KEEP,-PERC) %>%
  mutate("TAXON" = case_when(TAXON == "Bombina" ~ "Bom.",TAXON == "Lissotriton" ~ "Lis.",TAXON == "Myodes" ~"Myo.",
                             TAXON == "Anguis" ~ "Ang.",TAXON == "Solea" ~ "Sol.",T ~ TAXON)) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons1)) %>%
  rowwise() %>%
  mutate("HYBRID_ZONE" = case_when(grepl("montandoni",HYBRID_ZONE) ~ paste("montandoni - vulgaris (",stringr::str_sub(HYBRID_ZONE,-1),")",sep = ""),T ~ HYBRID_ZONE)) %>%
  ungroup()


signif_dt <- dt %>%
  group_by(HYBRID_ZONE,CLASS,P_VAL,TAXON) %>%
  summarise("max_shared" = max(PERC_VALUE)) %>%
  ungroup() %>% arrange(TAXON,CLASS,HYBRID_ZONE) %>%
  group_by(TAXON,CLASS) %>% mutate("order" = 1:n()) %>%
  mutate("signif" = P_VAL <= 0.05) %>%
  mutate("annotation" = case_when(signif == T & P_VAL <= 0.05 & P_VAL > 0.01 ~ "*", signif == T & P_VAL <= 0.01 & P_VAL > 0.001 ~ "**",signif == T & P_VAL <= 0.001 ~ "***"),
         "x" = order-.25,"xend" = order+.25,"y" = max_shared+1) %>%
  na.omit()

p <- ggplot(dt) + 
  geom_col(aes(x = HYBRID_ZONE,y = PERC_VALUE,fill = TYPE),position = position_dodge2(padding = 0))+
  facet_grid(CLASS~TAXON,scales = "free_x",space = "free")+
  labs(y = "Percent of shared MHC alleles")+
  theme_bw()+
  ylim(0,60)+
  #facet_wrap(vars(CLASS))+
  theme(axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1,face = "italic"),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        strip.text.x = element_text(size = 25,face = "italic"),strip.text.y = element_text(size = 25),
        axis.text = element_text(size = 25),axis.title.y = element_text(size = 30),
        strip.background = element_rect(fill = "white"),legend.title = element_text(size = 30),legend.text = element_text(size = 25),
        legend.key.size = unit(1,"cm"),strip.clip = "off")+
  scale_fill_hue(labels = c("parapatric","allopatric"),name = "Between")+
  ggsignif::geom_signif(stat = "identity",data = signif_dt,aes(x=x,xend=xend, y=y, yend=y, annotation=annotation,textsize = 8,group = HYBRID_ZONE))

ggsave("plots/Permutations.png",plot = p,device = "png",width = 22,height = 15)


#DRAW SES plot
dt <- read_and_rbind_perm(taxons) %>%
  filter(VARIANT == "all") %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = " - ")) %>%
  select(TAXON,HYBRID_ZONE,CLASS,SES) %>%
  mutate("CLASS" = paste("MHC",CLASS)) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons)) %>%
  mutate("SES" = case_when(SES == 0 ~ 0.1,T ~ SES)) %>%
  rowwise() %>%
  mutate("HYBRID_ZONE" = case_when(grepl("montandoni",HYBRID_ZONE) ~ paste("montandoni - vulgaris (",stringr::str_sub(HYBRID_ZONE,-1),")",sep = ""),T ~ HYBRID_ZONE)) %>%
  ungroup()

p <- ggplot(dt) + 
  geom_col(aes(x = HYBRID_ZONE,y = SES))+
  facet_grid(CLASS~TAXON,scales = "free_x",space = "free")+
  labs(y = "Standardized effect size")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        strip.text.x = element_text(size = 15,face = "italic"),strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 15),axis.title.y = element_text(size = 20),
        strip.background = element_rect(fill = "white"))

library(grid)
gt <- ggplot_gtable(ggplot_build(p))

#gtable::gtable_show_layout(gt)

#Function to make panels wider but keep the same width of bars
#panels_grobs - grobs' numbers corresponding to panels
#width_index - index of general grob widths gt$widths[?] #check by gtable::gtable_show_layout(gt)
#increament how much should we increase panel's width
# I made changes to posision because there is only a single bar here
increase_panels <- function(gt_list, panels_grobs = c(),width_index,increment = 1.5){
  gt <- gt_list
  panel_width <- as.numeric(gt$widths[width_index]) #collect original panel's width
  gt$widths[width_index] <- increment*gt$widths[width_index] #modify panels' width
  for(i in panels_grobs){
    wid_geom <- gt$grobs[[i]]$children[[3]]$width %>% as.numeric() %>% .[1] #collect original ration corresponding to bar's width
    
    wid_expected <- (panel_width*wid_geom)/(panel_width*increment) #calculate new ration of bar's width
    pos_expected <- unit(.5-(wid_expected/2),units = "native") #new position
    
    gt$grobs[[i]]$children[[3]]$width <- unit(c(wid_expected,wid_expected),units = "native") #change bars' widths
    gt$grobs[[i]]$children[[3]]$x <- pos_expected #change bars' positions
    
  }
  gt
}

gt <- increase_panels(gt_list = gt,panels_grobs = c(4,5),width_index = 7,increment = 1.2)
gt <- increase_panels(gt_list = gt,panels_grobs = c(10,11),width_index = 13,increment = 1.1)
gt <- increase_panels(gt_list = gt,panels_grobs = c(20,21),width_index = 23,increment = 1.2)

png("plots/SES.png",width = 18,height = 10,units = "in",res = 300)
grid.draw(gt)
dev.off()
#

### PERMUTATIONS - COMPARING P-VALUES - STATISTICAL MODEL ####
library(ape)
library(ggtree)
library(nlme)
taxons <- c("Solea","Bombina","Lissotriton","Triturus","Anguis","Emys","Natrix", "Podarcis","Sphyrapicus","Myodes")
read_and_rbind_perm <- function(taxons){
  taxs <- paste0("./permutations/",taxons,"_permutation_output_manual.csv")
  taxs[grepl("Lisso",taxs)] <- "./permutations/Lissotriton_permutation_output_manual_all.csv" #lissotriton averaged over transects
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

tree_org <- ape::read.tree("permutations/MHC_introgression_timetree.nwk")

# Time of divergence within species pair
# x <- cophenetic.phylo(tree_org) %>% as.data.frame() %>%
#   mutate("SP1" = row.names(.)) %>% relocate(SP1,.before = 1) %>%
#   pivot_longer(-SP1,names_to = "SP2",values_to = "DISTANCE") %>%
#   mutate("SP1" = sapply(sapply(SP1,strsplit,split = "_"),"[",2)) %>%
#   mutate("SP2" = sapply(sapply(SP2,strsplit,split = "_"),"[",2))
# write.csv(x,"temp_tables/DistanceTips_FromTree.csv",row.names = F,fileEncoding = "UTF-8")
dist <- read.csv("permutations/DistanceTips_FromTree.csv",fileEncoding = "UTF-8") #distance between tips

#plot.phylo(tree,show.node.label = T,show.tip.label = T)
#x <- as_tibble(tree_org)

#prepare data
dt <- read_and_rbind_perm(taxons) %>%
  filter(VARIANT == "all") %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = "_")) %>%
  select(TAXON,HYBRID_ZONE,CLASS,SES,SP1,SP2) %>%
  mutate("CLASS" = paste("MHC",CLASS)) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons)) %>%
  left_join(dist,c("SP1" = "SP1","SP2" = "SP2")) %>%
  mutate("DISTANCE" = if_else(is.na(DISTANCE) & SP1 == "o. occidentalis",0,DISTANCE)) %>%
  arrange(CLASS,HYBRID_ZONE)

#diagnose node leading to a given species pair
shared_node <- mrca(tree_org) %>% #common node for species pairs
  as.data.frame() %>%
  mutate("SP1" = row.names(.)) %>% relocate(SP1,.before = 1) %>%
  pivot_longer(-SP1,names_to = "SP2",values_to = "NODE") %>%
  mutate("SP1" = sapply(sapply(SP1,strsplit,split = "_"),"[",2)) %>%
  mutate("SP2" = sapply(sapply(SP2,strsplit,split = "_"),"[",2)) %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = "_")) %>%
  select(-SP1,-SP2) %>%
  mutate("NODE" = as.character(NODE)) %>%
  add_row("NODE" = "25","HYBRID_ZONE" = "o. occidentalis_orbicularis") #adding Emys

#calculate distance between nodes
dist_nodes <- dist.nodes(tree_org) %>% as.data.frame() %>%
  mutate("NODE1" = row.names(.)) %>% relocate(NODE1,.before = 1) %>%
  pivot_longer(-NODE1,names_to = "NODE2",values_to = "DISTANCE") %>%
  left_join(shared_node,by = c("NODE1" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE1" = "HYBRID_ZONE") %>%
  left_join(shared_node,by = c("NODE2" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE2" = "HYBRID_ZONE")

# Correlation matrix
cor_mat <- data.frame("HYBRID_ZONE1" = sort(rep(unique(dt$HYBRID_ZONE),times = length(unique(dt$HYBRID_ZONE)))),
                      "HYBRID_ZONE2" = rep(unique(dt$HYBRID_ZONE),times = length(unique(dt$HYBRID_ZONE)))) %>%
  left_join(dist_nodes,by = c("HYBRID_ZONE1","HYBRID_ZONE2"),multiple = "all") %>%
  select(-NODE1,-NODE2) %>% 
  pivot_wider(names_from = "HYBRID_ZONE2",values_from = "DISTANCE",names_sort = F) %>%
  select(-HYBRID_ZONE1)

rownames(cor_mat) <- colnames(cor_mat)
cor_mat <- as.matrix(cor_mat)

max_value <- max(cor_mat) #changing matrix of distances into correlation matrix
cor_mat <- (max_value-cor_mat)/max_value

#keep the same correlation matrix within each MHC class but assume no correlation between MHC classes
null_mat <- cor_mat
null_mat[,] <- 0
cor_mat <- cbind(rbind(cor_mat,null_mat),rbind(null_mat,cor_mat))


# C <- Matrix::nearPD(cor_mat)$mat #make matrix positive-definite
# C <- corSymm(C[lower.tri(C)], fixed = T) #wrapper to incorporate into gls model
#dt <- filter(dt,CLASS == "MHC II")

formula_v <- as.formula("SES ~ 0+CLASS") #change formula here

#test different lambda values
lambda_range <- seq(0,1,0.0001)
ll <- c()
for(lambda in lambda_range){
  if(lambda == 0){
    model_temp <- gls(formula_v,dt)
  } else {
    C_new <-  Matrix::nearPD(cor_mat*lambda)$mat
    C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
    model_temp <- gls(formula_v,dt,correlation = C_new)
  }
  ll[length(ll)+1] <- model_temp$logLik
}
chosen_lambda <- lambda_range[which.max(ll)]
 
#model without correlation
model0 <- gls(formula_v,dt)
summary(model0)

#model with correlation
C_new <-  Matrix::nearPD(cor_mat*chosen_lambda)$mat # modify matrix and make it positive-definite
C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
model1 <- gls(formula_v,dt,correlation = C_new)
summary(model1)

#plotting just CLASS (when formula is SES ~ 0+CLASS)
conf <-  as.data.frame(confint(model1))
dt_means <- data.frame("CLASS" = c("MHC I","MHC II"),"MEAN" = model1$coefficients)
dt_error <- data.frame("CLASS" = c("MHC I","MHC II"),"YMIN" = conf$`2.5 %`,"YMAX" = conf$`97.5 %`)

p <- ggplot()+
  geom_point(aes(x = CLASS,y = SES),data = dt,position = position_jitter(width = 0.15),size = 3)+
  geom_point(aes(x = CLASS,y = MEAN),data = dt_means,col = "red",size = 3)+
  theme_bw()+
  geom_errorbar(aes(x = CLASS,ymin = YMIN,ymax = YMAX),data = dt_error,color = "red",width = 0.2,linewidth = 2)+
  theme(axis.text = element_text(size = 20),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20))

ggsave(filename = "plots/SES_model_CLASS.png",device = "png",width = 5,height = 9)
#

### PARAMETERS OF GEOGRAPHIC CLINES ####
#prior to this step, I fitted geographic clines for each system independently
#in folder "../hzar_wrapper/ you will find R script, exemplary input file and output

taxons <- c("Bombina","Lissotriton","Triturus","Anguis","Emys","Sphyrapicus","Natrix","Solea")

read_and_rbind_perm <- function(taxons){
  taxs <- paste0("./geographic_clines/",taxons,"_LL_parameters.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

read_and_rbind_perm_n <- function(taxons){
  taxs <- paste0("./geographic_clines/",taxons,"_hzar_input.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

#get the number of localities
n_tab <- read_and_rbind_perm_n(taxons) %>%
  select(HYBRID_ZONE,"HI"=nSamples_H_I,"HII"=nSamples_H_II,"QI"=nSamples_Q_I,"Q2"=nSamples_Q_II,"genomic"=nSamples_genomic,TAXON) %>%
  pivot_longer(cols = 2:6,names_to = "MARKER",values_to = "OBSERVATIONS") %>% na.omit() %>%
  group_by(HYBRID_ZONE,TAXON,MARKER) %>% summarise("N" = n())

marker <- "H" #select marker (H - hybrid index, Q - structure)

dt <- read_and_rbind_perm(taxons) %>%
  left_join(n_tab,by = c("HYBRID_ZONE","TAXON","MARKER")) %>%
  filter(PARAMETERS %in% c("center","width")) %>% 
  filter(grepl(marker,MARKER) | MARKER == "genomic") %>%
  mutate("LEVEL" = case_when(LEVEL == "high" ~ "high2",T ~ LEVEL)) %>%
  mutate("LEVEL" = case_when(LEVEL == "low" ~ "high",LEVEL == "high2" ~ "low",T ~ LEVEL)) %>%
  mutate("MARKER" = case_when(grepl("II",MARKER) ~ "MHC II",grepl("I$",MARKER)  ~ "MHC I",T ~ MARKER)) %>%
  filter(!HYBRID_ZONE %in% c("astreptophora_helveticaAH","helvetica_natrixHN3"))

write.csv(dt,"temp_tables/Geographic_clines_original.csv",row.names = F,fileEncoding = "UTF-8")
#
### PARAMETERS OF GEOGRAPHIC CLINES - CALCULATING HEDGES ####
dt <- read.csv("temp_tables/Geographic_clines_original.csv",fileEncoding = "UTF-8") %>%
  pivot_wider(names_from = LEVEL,values_from = VALUES) %>%
  mutate("SE" = (high-low)/3.92) %>%
  mutate("SD" = SE*sqrt(N)) %>%
  mutate("SP1" = sapply(strsplit(HYBRID_ZONE,"_"),"[",1),"SP2" = sapply(strsplit(HYBRID_ZONE,"_"),"[",2)) %>%
  mutate(across(c(SP1,SP2),.fns = function(x){stringr::str_replace_all(x,pattern = "[:upper:].*$",replacement = "")}))

#correct signs of centre
allelic_richness <- read.csv("temp_tables/AllelicRichness_perSpecies.csv") %>%
  mutate("CLASS" = case_when(grepl("II",CLASS) ~ "MHC II",grepl("I$",CLASS)  ~ "MHC I",T ~ CLASS)) %>%
  select(-TAXON)

dt <- left_join(dt,allelic_richness,by = c("MARKER" = "CLASS","SP1" = "SPECIES")) %>%
  rename("AR1" = "ALL_RICHNESS") %>%
  left_join(allelic_richness,by = c("MARKER" = "CLASS","SP2" = "SPECIES")) %>% rename("AR2" = "ALL_RICHNESS") 

#Hedges G
out_list <- list()
for(i in unique(dt$HYBRID_ZONE)){
  temp1 <- filter(dt,HYBRID_ZONE == i)
  for(k in c("width","center")){
    temp <- filter(temp1,PARAMETERS == k)
    gen_value <- temp$main[temp$MARKER == "genomic"]
    gen_var <- temp$SD[temp$MARKER == "genomic"]^2
    gen_n <- temp$N[temp$MARKER == "genomic"]
    
    for(j in unique(temp$MARKER)){
      if(j != "genomic"){
        MHC_value <- temp$main[temp$MARKER == j]
        MHC_var <- temp$SD[temp$MARKER == j]^2
        MHC_n <- temp$N[temp$MARKER == j]
        
        J_estimate <- 1-(3/((4*(MHC_n+gen_n-2))-1))
        
        
        if(k == "width"){
          d <- ((MHC_value-gen_value)/sqrt((((MHC_n-1)*MHC_var)+((gen_n-1)*gen_var))/(MHC_n+gen_n-2)))*J_estimate
        } else if(k == "center"){
          diff_v <- MHC_value-gen_value
          if(temp$AR1[temp$MARKER == j] < temp$AR2[temp$MARKER == j]){
            diff_v <- -diff_v
          }
          d <- (diff_v/sqrt((((MHC_n-1)*MHC_var)+((gen_n-1)*gen_var))/(MHC_n+gen_n-2)))*J_estimate
        }
        v <- ((MHC_n+gen_n)/(MHC_n*gen_n))+((d^2)/(2*(MHC_n+gen_n)))
        out <- c(unique(temp$TAXON),i,j,J_estimate,d,v,k)
        out_list[[length(out_list)+1]] <- out
      }
    }
  }
}
out_table <- as.data.frame(do.call(rbind,out_list))
colnames(out_table) <- c("TAXON","HYBRID_ZONE","MARKER","J","d","v","PARAMETERS")
out_table <- mutate(out_table,across(.cols = c("J","d","v"),.fns = as.numeric))
write.csv(out_table,"temp_tables/Geographic_clines_Hedges_g.csv",fileEncoding = "UTF-8",row.names = F)

#plot
dt <- read.csv("temp_tables/Geographic_clines_Hedges_g.csv") %>%
  mutate("PARAMETERS" = if_else(PARAMETERS == "center","centre",PARAMETERS)) %>%
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegataK" ~ "Bombina bombina - B. variegata (K)",
                                   HYBRID_ZONE == "bombina_variegataJ" ~ "Bombina bombina - B. variegata (J)",
                                   HYBRID_ZONE == "montandoni_vulgarisS" ~ "Lissotriton montandoni - L. vulgaris (S)",
                                   HYBRID_ZONE == "montandoni_vulgarisT" ~ "Lissotriton montandoni - L. vulgaris (T)",
                                   HYBRID_ZONE == "montandoni_vulgarisR2" ~ "Lissotriton montandoni - L. vulgaris (R2)",
                                   HYBRID_ZONE == "montandoni_vulgarisR" ~ "Lissotriton montandoni - L. vulgaris (R)",
                                   HYBRID_ZONE == "montandoni_vulgarisL" ~ "Lissotriton montandoni - L. vulgaris (L)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilisA" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "occidentalis_orbicularisA" ~ "Emys o. occidentalis - E. orbicularis",
                                   HYBRID_ZONE == "nuchalis_ruberNR" ~ "Sphyrapicus nuchalis - S. ruber",
                                   #HYBRID_ZONE == "helvetica_natrixHN3" ~ "Natrix helvetica - N. natrix (HN3)",
                                   HYBRID_ZONE == "helvetica_natrixHN2" ~ "Natrix helvetica - N. natrix (HN2)",
                                   HYBRID_ZONE == "helvetica_natrixHN4" ~ "Natrix helvetica - N. natrix (HN1)",
                                   #HYBRID_ZONE == "astreptophora_helveticaAH" ~ "Natrix astreptophora - N. helvetica",
                                   HYBRID_ZONE == "aegyptiaca_senegalensisA" ~ "Solea aegyptiaca - S. senegalensis"))

dt$HYBRID_ZONE <- factor(dt$HYBRID_ZONE,levels = rev(c("Solea aegyptiaca - S. senegalensis",
                                                   "Bombina bombina - B. variegata (J)","Bombina bombina - B. variegata (K)",
                                                   "Lissotriton montandoni - L. vulgaris (L)","Lissotriton montandoni - L. vulgaris (R)",
                                                   "Lissotriton montandoni - L. vulgaris (R2)",
                                                   "Lissotriton montandoni - L. vulgaris (S)","Lissotriton montandoni - L. vulgaris (T)",
                                                   "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                   "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                   "Anguis colchica - A. fragilis","Emys o. occidentalis - E. orbicularis",
                                                   #"Natrix astreptophora - N. helvetica",
                                                   "Natrix helvetica - N. natrix (HN1)",
                                                   "Natrix helvetica - N. natrix (HN2)",
                                                   #"Natrix helvetica - N. natrix (HN3)",
                                                   "Sphyrapicus nuchalis - S. ruber")))

dt$PARAMETERS <- factor(dt$PARAMETERS,levels = c("width","centre"))

p <- ggplot(dt) + 
  geom_errorbarh(aes(xmin = d-(1.96*sqrt(v)),xmax = d+(1.96*sqrt(v)),y = HYBRID_ZONE,color = MARKER),
                 position = position_dodge(width = 0.4),height = 0.5) +
  geom_point(aes(x = d,y = HYBRID_ZONE,color = MARKER),position = position_dodge(width = 0.4))+
  facet_grid(cols =  vars(PARAMETERS))+
  theme_bw()+labs(x = "Effect size (Hedges g)")+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size = 20),panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),panel.grid.minor = element_blank(),legend.title = element_blank(),
        legend.text = element_text(size = 15),strip.background = element_rect(fill = "white"),strip.text = element_text(size = 15))+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3"))+
  scale_y_discrete(expand = expansion(mult = c(0.05,0.15)))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")

ggsave("plots/Effect_sizes_geoographicClines.png",device = "png",width = 15,height = 9,plot = p)
#
### GEOGRAPHIC CLINES - STATISTICAL MODEL ####
library(ape)
library(ggtree)
library(nlme)
taxons <- c("Solea","Bombina","Lissotriton","Triturus","Anguis","Emys","Natrix","Sphyrapicus")

tree_org <- ape::read.tree("permutations/MHC_introgression_timetree.nwk")

# Time of divergence within species pair
dist <- read.csv("permutations/DistanceTips_FromTree.csv",fileEncoding = "UTF-8") #distance between tips

#prepare data
dt1 <- read.csv("temp_tables/Geographic_clines_Hedges_g.csv") %>%
  mutate("HYBRID_ZONE" = stringr::str_replace_all(HYBRID_ZONE,pattern = "[:upper:].*$",replacement = "")) %>%
  group_by(TAXON,HYBRID_ZONE,MARKER,PARAMETERS) %>%
  summarise("d" = mean(d),"v" = mean(v)) %>%
  mutate("SP1" = sapply(strsplit(HYBRID_ZONE,"_"),"[",1),"SP2" = sapply(strsplit(HYBRID_ZONE,"_"),"[",2)) %>%
  select(TAXON,HYBRID_ZONE,"CLASS" = MARKER,PARAMETERS,d,v,SP1,SP2) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons)) %>%
  left_join(dist,c("SP1" = "SP1","SP2" = "SP2")) %>%
  mutate("DISTANCE" = if_else(is.na(DISTANCE) & SP1 == "occidentalis",0,DISTANCE)) %>%
  arrange(CLASS,HYBRID_ZONE) %>%
  mutate("CLASS" = stringr::str_replace_all(CLASS,"-"," ")) #%>%
  #filter(!HYBRID_ZONE %in% c("astreptophora_helvetica","helveticaNH3_natrixHN3"))

#diagnose node leading to a given species pair

shared_node <- mrca(tree_org) %>% #common node for species pairs
  as.data.frame() %>%
  mutate("SP1" = row.names(.)) %>% relocate(SP1,.before = 1) %>%
  pivot_longer(-SP1,names_to = "SP2",values_to = "NODE") %>%
  mutate("SP1" = sapply(sapply(SP1,strsplit,split = "_"),"[",2)) %>%
  mutate("SP2" = sapply(sapply(SP2,strsplit,split = "_"),"[",2)) %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = "_")) %>%
  select(-SP1,-SP2) %>%
  mutate("NODE" = as.character(NODE)) %>%
  add_row("NODE" = "25","HYBRID_ZONE" = "occidentalis_orbicularis") # adding Emys
  
#calculate distance between nodes
dist_nodes <- dist.nodes(tree_org) %>% as.data.frame() %>%
  mutate("NODE1" = row.names(.)) %>% relocate(NODE1,.before = 1) %>%
  pivot_longer(-NODE1,names_to = "NODE2",values_to = "DISTANCE") %>%
  left_join(shared_node,by = c("NODE1" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE1" = "HYBRID_ZONE") %>%
  left_join(shared_node,by = c("NODE2" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE2" = "HYBRID_ZONE")

# Correlation matrix
cor_mat <- data.frame("HYBRID_ZONE1" = sort(rep(unique(dt1$HYBRID_ZONE),times = length(unique(dt1$HYBRID_ZONE)))),
                      "HYBRID_ZONE2" = rep(unique(dt1$HYBRID_ZONE),times = length(unique(dt1$HYBRID_ZONE)))) %>%
  left_join(dist_nodes,by = c("HYBRID_ZONE1","HYBRID_ZONE2"),multiple = "all") %>%
  select(-NODE1,-NODE2) %>% 
  pivot_wider(names_from = "HYBRID_ZONE2",values_from = "DISTANCE",names_sort = F) %>%
  select(-HYBRID_ZONE1)

rownames(cor_mat) <- colnames(cor_mat)
cor_mat <- as.matrix(cor_mat)

max_value <- max(cor_mat) #changing matrix of distances into correlation matrix
cor_mat <- (max_value-cor_mat)/max_value

#keep the same correlation matrix within each MHC class but assume no correlation between MHC classes
null_mat <- cor_mat
null_mat[,] <- 0
cor_mat <- cbind(rbind(cor_mat,null_mat),rbind(null_mat,cor_mat))


#C <- Matrix::nearPD(cor_mat)$mat #make matrix positive-definite
#C <- corSymm(C[lower.tri(C)], fixed = T) #wrapper to incorporate into gls model
#dt <- filter(dt,CLASS == "MHC II")
#m0 <- gls(d ~ 0+CLASS,dt,weights = ~(1/v),correlation = C)

parameter <- "width" #either width or centre

dt <- filter(dt1,PARAMETERS == parameter)

formula_v <- as.formula("d ~ 0+CLASS")

#test different lambda values
lambda_range <- seq(0,1,0.0001)
ll <- c()
for(lambda in lambda_range){
  if(lambda == 0){
    model_temp <- gls(formula_v,dt)
  } else {
    C_new <-  Matrix::nearPD(cor_mat*lambda)$mat
    C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
    model_temp <- gls(formula_v,dt,correlation = C_new,weights = ~(1/v))
  }
  ll[length(ll)+1] <- model_temp$logLik
}
chosen_lambda <- lambda_range[which.max(ll)]

#model without correlation
model0 <- gls(formula_v,dt)
summary(model0)

#model with correlation
if(chosen_lambda != 0){
C_new <-  Matrix::nearPD(cor_mat*chosen_lambda)$mat # modify matrix and make it positive-definite
C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
model1 <- gls(formula_v,dt,correlation = C_new,weights = ~(1/v))
summary(model1)
} else {
  model1 <- gls(formula_v,dt,weights = ~(1/v))
  summary(model1)
}

#plotting just CLASS
conf <-  as.data.frame(nlme::intervals(model1)$coef)
dt_means <- data.frame("CLASS" = c("MHC I","MHC II"),"MEAN" = model1$coefficients)
dt_error <- data.frame("CLASS" = c("MHC I","MHC II"),"YMIN" = conf$lower,"YMAX" = conf$upper)

p <- ggplot()+
  geom_point(aes(x = CLASS,y = d),data = dt,position = position_jitter(width = 0.15),size = 3)+
  geom_point(aes(x = CLASS,y = MEAN),data = dt_means,col = "red",size = 3)+
  theme_bw()+
  geom_errorbar(aes(x = CLASS,ymin = YMIN,ymax = YMAX),data = dt_error,color = "red",width = 0.2,linewidth = 2)+
  theme(axis.text = element_text(size = 20),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20))+
  labs(y = "Hedges g")+
  scale_y_continuous(limits = c(-0.9,1.4),breaks = c(-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1,1.2),expand = c(0,0))

ggsave(filename = "plots/Geographic_Hedges_width_model_CLASS.png",device = "png",width = 5,height = 9)


# plotting full model including divergence time
dist_seq <- seq(0,20,0.01)
dt_predict <- data.frame("CLASS" = sort(rep(c("MHC I","MHC II"),length(dist_seq))),"DISTANCE" = c(dist_seq,dist_seq))
dt1 <- dt1 %>% mutate("PARAMETERS" = if_else(PARAMETERS == "center","centre",PARAMETERS))


formula_v <- as.formula("d ~ CLASS*DISTANCE")
model_width <- gls(formula_v,filter(dt1,PARAMETERS == "width"),weights = ~(1/v)) #at this point we know there is no correlation
model_center <- gls(formula_v,filter(dt1,PARAMETERS == "centre"),weights = ~(1/v))

predicted_width <- predict(model_width,newdata = dt_predict)
predicted_center <- predict(model_center,newdata = dt_predict)

dt_predict <- cbind(rbind(dt_predict,dt_predict),"d" = c(predicted_width,predicted_center),
                    "PARAMETERS" = c(rep("width",2*length(dist_seq)),rep("centre",2*length(dist_seq))))

p <- ggplot()+
  geom_point(aes(x = DISTANCE,y = d,color = CLASS),data = dt1,size = 3)+
  geom_line(aes(x = DISTANCE,y = d,color = CLASS),data = dt_predict)+
  theme_bw()+
  labs(y = "Effect size (Hedges g)",x = "Time of divergence [My]")+
  facet_grid(cols = vars(PARAMETERS))+
  theme(axis.title = element_text(size = 20),panel.grid = element_blank(),
        axis.text = element_text(size = 15),legend.title = element_blank(),
        legend.text = element_text(size = 15),strip.background = element_rect(fill = "white"),strip.text = element_text(size = 15))+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3"))

ggsave(filename = "plots/Geographic_clines_Hedges_modelsDivergence.png",device = "png",width = 9,height = 9)
#
### PLOTTING GEOGRAPHIC CLINES ####
taxons <- c("Solea","Bombina","Lissotriton","Triturus","Anguis","Emys","Natrix","Sphyrapicus")

combine_all_parameters <- function(taxons){
  taxs <- paste0("./geographic_clines/",taxons,"_parameters.RDS")
  out <- lapply(taxs,readRDS)
  unlist(out,recursive = F)
}

df <- combine_all_parameters(taxons) %>% .[!names(.) %in% c("astreptophora_helveticaAH","helvetica_natrixHN3")] 
# ^ those were excluded as the null model (horizontal line) was optimal for at least one MHC class

hzs <- names(df)
markers <- c("HI","HII","genomic")

out_list_points <- list()
out_list_lines <- list()
for(i in hzs){
  for(j in markers){
    temp <- df[[i]][[j]]
    dist <- temp$dist_xs
    lines_tab <- temp$CI_data %>%
      rename("dist" = "x") %>%
      mutate("obsFreq" = temp$function_cline(dist),"HYBRID_ZONE" = i,"marker" = j)
    points_tab  <- temp$obs %>% mutate("marker" = j,"HYBRID_ZONE" = i)
    out_list_points[[length(out_list_points)+1]] <- points_tab
    out_list_lines[[length(out_list_lines)+1]] <- lines_tab
  }
}

points_tab <- as.data.frame(do.call(rbind,out_list_points)) %>%
  rename("CLASS" = "marker") %>%
  mutate("CLASS" = case_when(CLASS == "HI" ~ "MHC I",CLASS == "HII" ~ "MHC II",CLASS == "genomic" ~ "genome-wide")) %>% 
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegataK" ~ "Bombina bombina - B. variegata (K)",
                                   HYBRID_ZONE == "bombina_variegataJ" ~ "Bombina bombina - B. variegata (J)",
                                   HYBRID_ZONE == "montandoni_vulgarisS" ~ "Lissotriton montandoni - L. vulgaris (S)",
                                   HYBRID_ZONE == "montandoni_vulgarisT" ~ "Lissotriton montandoni - L. vulgaris (T)",
                                   HYBRID_ZONE == "montandoni_vulgarisR2" ~ "Lissotriton montandoni - L. vulgaris (R2)",
                                   HYBRID_ZONE == "montandoni_vulgarisR" ~ "Lissotriton montandoni - L. vulgaris (R)",
                                   HYBRID_ZONE == "montandoni_vulgarisL" ~ "Lissotriton montandoni - L. vulgaris (L)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilis" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "occidentalis_orbicularisA" ~ "Emys o. occidentalis - E. orbicularis",
                                   HYBRID_ZONE == "nuchalis_ruberNR" ~ "Sphyrapicus nuchalis - S. ruber",
                                   #HYBRID_ZONE == "helvetica_natrixHN3" ~ "Natrix helvetica - N. natrix (HN3)",
                                   HYBRID_ZONE == "helvetica_natrixHN2" ~ "Natrix helvetica - N. natrix (HN2)",
                                   HYBRID_ZONE == "helvetica_natrixHN4" ~ "Natrix helvetica - N. natrix (HN1)",
                                   #HYBRID_ZONE == "astreptophora_helveticaAH" ~ "Natrix astreptophora - N. helvetica",
                                   HYBRID_ZONE == "aegyptiaca_senegalensis" ~ "Solea aegyptiaca - S. senegalensis"))

points_tab$HYBRID_ZONE <-  factor(points_tab$HYBRID_ZONE,levels = c("Solea aegyptiaca - S. senegalensis",
                                                                        "Bombina bombina - B. variegata (J)","Bombina bombina - B. variegata (K)",
                                                                        "Lissotriton montandoni - L. vulgaris (L)","Lissotriton montandoni - L. vulgaris (R)",
                                                                        "Lissotriton montandoni - L. vulgaris (R2)",
                                                                        "Lissotriton montandoni - L. vulgaris (S)","Lissotriton montandoni - L. vulgaris (T)",
                                                                        "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                                        "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                                        "Anguis colchica - A. fragilis","Emys o. occidentalis - E. orbicularis",
                                                                        #"Natrix astreptophora - N. helvetica",
                                                                        "Natrix helvetica - N. natrix (HN1)",
                                                                        "Natrix helvetica - N. natrix (HN2)",
                                                                        #"Natrix helvetica - N. natrix (HN3)",
                                                                        "Sphyrapicus nuchalis - S. ruber"))

lines_tab <- as.data.frame(do.call(rbind,out_list_lines)) %>%
  rename("CLASS" = "marker") %>%
  mutate("CLASS" = case_when(CLASS == "HI" ~ "MHC I",CLASS == "HII" ~ "MHC II",CLASS == "genomic" ~ "genome-wide")) %>% 
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegataK" ~ "Bombina bombina - B. variegata (K)",
                                   HYBRID_ZONE == "bombina_variegataJ" ~ "Bombina bombina - B. variegata (J)",
                                   HYBRID_ZONE == "montandoni_vulgarisS" ~ "Lissotriton montandoni - L. vulgaris (S)",
                                   HYBRID_ZONE == "montandoni_vulgarisT" ~ "Lissotriton montandoni - L. vulgaris (T)",
                                   HYBRID_ZONE == "montandoni_vulgarisR2" ~ "Lissotriton montandoni - L. vulgaris (R2)",
                                   HYBRID_ZONE == "montandoni_vulgarisR" ~ "Lissotriton montandoni - L. vulgaris (R)",
                                   HYBRID_ZONE == "montandoni_vulgarisL" ~ "Lissotriton montandoni - L. vulgaris (L)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilis" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "occidentalis_orbicularisA" ~ "Emys o. occidentalis - E. orbicularis",
                                   HYBRID_ZONE == "nuchalis_ruberNR" ~ "Sphyrapicus nuchalis - S. ruber",
                                   #HYBRID_ZONE == "helvetica_natrixHN3" ~ "Natrix helvetica - N. natrix (HN3)",
                                   HYBRID_ZONE == "helvetica_natrixHN2" ~ "Natrix helvetica - N. natrix (HN2)",
                                   HYBRID_ZONE == "helvetica_natrixHN4" ~ "Natrix helvetica - N. natrix (HN1)",
                                   #HYBRID_ZONE == "astreptophora_helveticaAH" ~ "Natrix astreptophora - N. helvetica",
                                   HYBRID_ZONE == "aegyptiaca_senegalensis" ~ "Solea aegyptiaca - S. senegalensis"))



lines_tab$HYBRID_ZONE <-  factor(lines_tab$HYBRID_ZONE,levels = c("Solea aegyptiaca - S. senegalensis",
                                                                       "Bombina bombina - B. variegata (J)","Bombina bombina - B. variegata (K)",
                                                                       "Lissotriton montandoni - L. vulgaris (L)","Lissotriton montandoni - L. vulgaris (R)",
                                                                       "Lissotriton montandoni - L. vulgaris (R2)",
                                                                       "Lissotriton montandoni - L. vulgaris (S)","Lissotriton montandoni - L. vulgaris (T)",
                                                                       "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                                       "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                                       "Anguis colchica - A. fragilis","Emys o. occidentalis - E. orbicularis",
                                                                       #"Natrix astreptophora - N. helvetica",
                                                                       "Natrix helvetica - N. natrix (HN1)",
                                                                       "Natrix helvetica - N. natrix (HN2)",
                                                                       #"Natrix helvetica - N. natrix (HN3)",
                                                                       "Sphyrapicus nuchalis - S. ruber"))
#solving for Triturus so it also starts at 0
lines_tab <- lines_tab %>% group_by(HYBRID_ZONE) %>%
  mutate("dist" = dist-min(dist))
points_tab <- points_tab %>%
  group_by(HYBRID_ZONE) %>%
  mutate("dist" = dist-min(dist))

my_plot <- ggplot(lines_tab)+ theme_bw() +
  labs(x = "Distance [km]",y = "Ancestry estimate")+
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 25),
        legend.position = c(0.92,0.10),legend.title = element_blank(),legend.key.size = unit(2,units = "cm"),
        legend.text = element_text(size = 30),strip.text = ggtext::element_markdown(size = 18),
        strip.background = element_rect(fill = "white"),plot.margin = unit(c(1,1,1,1),"cm"),panel.spacing = unit(0.5,"cm")) +
  geom_point(data = points_tab,aes(x = dist,y = obsFreq,color = CLASS),alpha = 0.3,size = 3)+
  geom_ribbon(aes(x = dist,ymin = yMin,ymax = yMax,fill = CLASS),alpha = 0.3)+
  geom_line(aes(x = dist,y = obsFreq,col = CLASS),size = 1,alpha = 0.3)+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3","genome-wide" = "red"))+
  scale_fill_manual(values = c("MHC I" = "black","MHC II" = "green3","genome-wide" = "red"))+
  facet_wrap(vars(HYBRID_ZONE),scales = "free_x")

ggsave("plots/Drawn_geographic_clines.png",my_plot,device = "png",width = 27,height = 18)
#
### GENOMIC CLINES - STATISTICAL MODEL ####
#all files used at this stage were generated independently for each system
#you can find the code in "../genomicClines/Genomic_clines_code_parallel.R" together with exemplary input file and output files

library(ape)
library(ggtree)
library(nlme)
taxons <- c("Bombina","Lissotriton","Triturus","Anguis","Natrix","Sphyrapicus")

tree_org <- ape::read.tree("permutations/MHC_introgression_timetree.nwk")

dist <- read.csv("permutations/DistanceTips_FromTree.csv",fileEncoding = "UTF-8") #distance between tips

#prepare data
read_and_rbind_perm <- function(taxons){
  taxs <- paste0("./genomic_clines/",taxons,"_GenomicClines_parameters.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

#correct signs of centre
allelic_richness <- read.csv("temp_tables/AllelicRichness_perSpecies.csv") %>%
  mutate("CLASS" = case_when(grepl("II",CLASS) ~ "MHC II",grepl("I$",CLASS)  ~ "MHC I",T ~ CLASS)) %>%
  select(-TAXON)

dt1 <- read_and_rbind_perm(taxons) %>%
  mutate("SE" = (UPPER_2LL-LOWER_2LL)/3.92) %>%
  mutate("VAR" = SE^2) %>%
  mutate("HYBRID_ZONE" = stringr::str_replace_all(HYBRID_ZONE,pattern = "[:upper:].*$",replacement = "")) %>%
  group_by(TAXON,HYBRID_ZONE,CLASS,PARAMETERS) %>%
  summarise("d" = mean(VALUE),"v" = mean(VAR)) %>%
  mutate("SP1" = sapply(strsplit(HYBRID_ZONE,"_"),"[",1),"SP2" = sapply(strsplit(HYBRID_ZONE,"_"),"[",2)) %>%
  select(TAXON,HYBRID_ZONE,CLASS,PARAMETERS,d,v,SP1,SP2) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons)) %>%
  left_join(dist,c("SP1" = "SP1","SP2" = "SP2")) %>%
  arrange(CLASS,HYBRID_ZONE) %>%
  mutate("CLASS" = paste("MHC",CLASS)) %>%
  left_join(allelic_richness, by = c("CLASS" = "CLASS","SP1" = "SPECIES")) %>%
  rename("AR1" = "ALL_RICHNESS") %>%
  left_join(allelic_richness, by = c("CLASS" = "CLASS","SP2" = "SPECIES")) %>%
  rename("AR2" = "ALL_RICHNESS") %>%
  mutate("d" = case_when(PARAMETERS == "ALPHA" & AR1 > AR2 ~ -d,T ~ d))

write.csv(dt1,"temp_tables/GenomicClines_tableModels.csv",fileEncoding = "UTF-8",row.names = F)

param <- "BETA" #either alpha or beta parameter
dt <- read.csv("temp_tables/GenomicClines_tableModels.csv") %>% filter(PARAMETERS == param)

shared_node <- mrca(tree_org) %>% #common node for species pairs
  as.data.frame() %>%
  mutate("SP1" = row.names(.)) %>% relocate(SP1,.before = 1) %>%
  pivot_longer(-SP1,names_to = "SP2",values_to = "NODE") %>%
  mutate("SP1" = sapply(sapply(SP1,strsplit,split = "_"),"[",2)) %>%
  mutate("SP2" = sapply(sapply(SP2,strsplit,split = "_"),"[",2)) %>%
  mutate("HYBRID_ZONE" = paste(SP1,SP2,sep = "_")) %>%
  select(-SP1,-SP2) %>%
  mutate("NODE" = as.character(NODE)) %>%
  add_row("NODE" = "25","HYBRID_ZONE" = "occidentalis_orbicularis") #%>% #adding Emys

#calculate distance between nodes
dist_nodes <- dist.nodes(tree_org) %>% as.data.frame() %>%
  mutate("NODE1" = row.names(.)) %>% relocate(NODE1,.before = 1) %>%
  pivot_longer(-NODE1,names_to = "NODE2",values_to = "DISTANCE") %>%
  left_join(shared_node,by = c("NODE1" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE1" = "HYBRID_ZONE") %>%
  left_join(shared_node,by = c("NODE2" = "NODE"),multiple = "all") %>%
  rename("HYBRID_ZONE2" = "HYBRID_ZONE")

# Correlation matrix
cor_mat <- data.frame("HYBRID_ZONE1" = sort(rep(unique(dt$HYBRID_ZONE),times = length(unique(dt$HYBRID_ZONE)))),
                      "HYBRID_ZONE2" = rep(unique(dt$HYBRID_ZONE),times = length(unique(dt$HYBRID_ZONE)))) %>%
  left_join(dist_nodes,by = c("HYBRID_ZONE1","HYBRID_ZONE2"),multiple = "all") %>%
  select(-NODE1,-NODE2) %>% 
  pivot_wider(names_from = "HYBRID_ZONE2",values_from = "DISTANCE",names_sort = F) %>%
  select(-HYBRID_ZONE1)

rownames(cor_mat) <- colnames(cor_mat)
cor_mat <- as.matrix(cor_mat)

max_value <- max(cor_mat) #changing matrix of distances into correlation matrix
cor_mat <- (max_value-cor_mat)/max_value

#keep the same correlation matrix within each MHC class but assume no correlation between MHC classes
null_mat <- cor_mat
null_mat[,] <- 0
cor_mat <- cbind(rbind(cor_mat,null_mat),rbind(null_mat,cor_mat))

formula_v <- as.formula("d ~ 0+CLASS")

#test different lambda values
lambda_range <- seq(0,1,0.0001)
ll <- c()
for(lambda in lambda_range){
  if(lambda == 0){
    model_temp <- gls(formula_v,dt)
  } else {
    C_new <-  Matrix::nearPD(cor_mat*lambda)$mat
    C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
    model_temp <- gls(formula_v,dt,correlation = C_new,weights = ~(1/v))
  }
  ll[length(ll)+1] <- model_temp$logLik
}
chosen_lambda <- lambda_range[which.max(ll)]

#model without correlation
model0 <- gls(formula_v,dt)
summary(model0)

#model with correlation
if(chosen_lambda != 0){
  C_new <-  Matrix::nearPD(cor_mat*chosen_lambda)$mat # modify matrix and make it positive-definite
  C_new <- corSymm(C_new[lower.tri(C_new)], fixed = T)
  model1 <- gls(formula_v,dt,correlation = C_new,weights = ~(1/v))
  summary(model1)
} else {
  model1 <- gls(formula_v,dt,weights = ~(1/v))
  summary(model1)
}

#plotting just CLASS (when formula is d ~ 0+CLASS)
conf <-  as.data.frame(nlme::intervals(model1)$coef)
dt_means <- data.frame("CLASS" = c("MHC I","MHC II"),"MEAN" = model1$coefficients)
dt_error <- data.frame("CLASS" = c("MHC I","MHC II"),"YMIN" = conf$lower,"YMAX" = conf$upper)

p <- ggplot()+
  geom_point(aes(x = CLASS,y = d),data = dt,position = position_jitter(width = 0.15),size = 3)+
  geom_point(aes(x = CLASS,y = MEAN),data = dt_means,col = "red",size = 3)+
  theme_bw()+
  geom_errorbar(aes(x = CLASS,ymin = YMIN,ymax = YMAX),data = dt_error,color = "red",width = 0.2,linewidth = 2)+
  theme(axis.text = element_text(size = 20),panel.grid.major.x = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20))+
  labs(y = "β")+ylim(-1.05,1.05)#α
  #scale_y_continuous(limits = c(-0.4,0.7),breaks = c(-.4,-.2,0,.2,.4,.6),expand = c(0,0))

ggsave(filename = "plots/Genomic_clines_Beta_model_CLASS.png",device = "png",width = 5,height = 9)


# plotting full model including divergence
dist_seq <- seq(0,20,0.01)
dt_predict <- data.frame("CLASS" = sort(rep(c("MHC I","MHC II"),length(dist_seq))),"DISTANCE" = c(dist_seq,dist_seq))

dt <- read.csv("temp_tables/GenomicClines_tableModels.csv") %>%
  filter(PARAMETERS %in% c("BETA","ALPHA")) %>%
  mutate("PARAMETERS" = case_when(PARAMETERS == "BETA" ~ "β",T ~ "α"))

formula_v <- as.formula("d ~ CLASS*DISTANCE")

model_beta <- gls(formula_v,filter(dt,PARAMETERS == "β"),weights = ~(1/v))

model_alpha <- gls(formula_v,filter(dt,PARAMETERS == "α"),weights = ~(1/v))


predicted_beta <- predict(model_beta,newdata = dt_predict)
predicted_alpha <- predict(model_alpha,newdata = dt_predict)

dt_predict <- cbind(rbind(dt_predict,dt_predict),"d" = c(predicted_beta,predicted_alpha),
                    "PARAMETERS" = c(rep("β",2*length(dist_seq)),rep("α",2*length(dist_seq))))

p <- ggplot()+
  geom_point(aes(x = DISTANCE,y = d,color = CLASS),data = dt,size = 3)+
  geom_line(aes(x = DISTANCE,y = d,color = CLASS),data = dt_predict)+
  theme_bw()+
  labs(y = "Value",x = "Time of divergence [My]")+
  facet_grid(cols = vars(PARAMETERS))+
  theme(axis.title = element_text(size = 20),panel.grid = element_blank(),
        axis.text = element_text(size = 15),legend.title = element_blank(),
        legend.text = element_text(size = 15),strip.background = element_rect(fill = "white"),strip.text = element_text(size = 15))+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3"))

ggsave(filename = "plots/Genomic_modelsDivergence.png",device = "png",width = 9,height = 9)
#
### GENOMIC CLINES - OVERALL PLOT ####
taxons <- c("Bombina","Lissotriton","Triturus","Anguis","Natrix","Sphyrapicus")

read_and_rbind_perm <- function(taxons){
  taxs <-  paste0("./genomic_clines/",taxons,"_GenomicClines_parameters.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

allelic_richness <- read.csv("temp_tables/AllelicRichness_perSpecies.csv") %>%
  mutate("CLASS" = case_when(grepl("II",CLASS) ~ "MHC II",grepl("I$",CLASS)  ~ "MHC I",T ~ CLASS)) %>%
  select(-TAXON)

dt <- read_and_rbind_perm(taxons) %>%
  mutate("SE" = (UPPER_2LL-LOWER_2LL)/3.92) %>%
  mutate("v" = SE^2,"d" = VALUE) %>%
  mutate("SP1" = sapply(strsplit(HYBRID_ZONE,"_"),"[",1),"SP2" = sapply(strsplit(HYBRID_ZONE,"_"),"[",2)) %>%
  mutate("SP1" = stringr::str_replace_all(SP1,pattern = "[:upper:].*$",replacement = ""),
         "SP2" = stringr::str_replace_all(SP2,pattern = "[:upper:].*$",replacement = "")) %>%
  select(TAXON,HYBRID_ZONE,CLASS,PARAMETERS,d,v,SP1,SP2,UPPER_2LL,LOWER_2LL,SE) %>%
  mutate("TAXON" = factor(TAXON,levels = taxons)) %>%
  arrange(CLASS,HYBRID_ZONE) %>%
  mutate("CLASS" = paste("MHC",CLASS)) %>%
  left_join(allelic_richness, by = c("CLASS" = "CLASS","SP1" = "SPECIES")) %>%
  rename("AR1" = "ALL_RICHNESS") %>%
  left_join(allelic_richness, by = c("CLASS" = "CLASS","SP2" = "SPECIES")) %>%
  rename("AR2" = "ALL_RICHNESS") %>%
  mutate("d" = case_when(PARAMETERS == "ALPHA" & AR1 > AR2 ~ -d,T ~ d),
         "UPPER_2LL" = case_when(PARAMETERS == "ALPHA" & AR1 > AR2 ~ -UPPER_2LL,T ~ UPPER_2LL),
         "LOWER_2LL" = case_when(PARAMETERS == "ALPHA" & AR1 > AR2 ~ -LOWER_2LL,T ~ LOWER_2LL)) %>%
  filter(PARAMETERS %in% c("BETA","ALPHA")) %>%
  mutate("PARAMETERS" = case_when(PARAMETERS == "BETA" ~ "β",T ~ "α")) %>%
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegata" ~ "Bombina bombina - B. variegata",
                                   HYBRID_ZONE == "montandoni_vulgarisIN" ~ "Lissotriton montandoni - L. vulgaris (IN)",
                                   HYBRID_ZONE == "montandoni_vulgarisOUT" ~ "Lissotriton montandoni - L. vulgaris (OUT)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilis" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "nuchalis_ruber" ~ "Sphyrapicus nuchalis - S. ruber",
                                   HYBRID_ZONE == "ruber_varius" ~ "Sphyrapicus ruber - S. varius",
                                   HYBRID_ZONE == "nuchalis_varius" ~ "Sphyrapicus nuchalis - S. varius",
                                   HYBRID_ZONE == "helvetica_natrix" ~ "Natrix helvetica - N. natrix",
                                   HYBRID_ZONE == "astreptophora_helvetica" ~ "Natrix astreptophora - N. helvetica"))

dt$HYBRID_ZONE <- factor(dt$HYBRID_ZONE,levels = rev(c("Bombina bombina - B. variegata",
                                                       "Lissotriton montandoni - L. vulgaris (IN)",
                                                       "Lissotriton montandoni - L. vulgaris (OUT)",
                                                       "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                       "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                       "Anguis colchica - A. fragilis",
                                                       "Natrix astreptophora - N. helvetica",
                                                       "Natrix helvetica - N. natrix",
                                                       "Sphyrapicus nuchalis - S. ruber",
                                                       "Sphyrapicus nuchalis - S. varius",
                                                       "Sphyrapicus ruber - S. varius")))

dt$PARAMETERS <- factor(dt$PARAMETERS,levels = c("β","α"))

p <- ggplot(dt) + 
  geom_errorbarh(aes(xmin = LOWER_2LL,xmax = UPPER_2LL,y = HYBRID_ZONE,color = CLASS),
  #geom_errorbarh(aes(xmin = d-(SE*1.96),xmax = d+(SE*1.96),y = HYBRID_ZONE,color = CLASS),
                 position = position_dodge(width = 0.4),height = 0.5) +
  geom_point(aes(x = d,y = HYBRID_ZONE,color = CLASS),position = position_dodge(width = 0.4))+
  facet_grid(cols =  vars(PARAMETERS))+
  theme_bw()+labs(x = "Value")+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size = 20),panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),panel.grid.minor = element_blank(),legend.title = element_blank(),
        legend.text = element_text(size = 15),strip.background = element_rect(fill = "white"),strip.text = element_text(size = 15))+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3"))+
  scale_y_discrete(expand = expansion(mult = c(0.05,0.15)))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")

ggsave("plots/Parameters_genomicClines_LL.png",device = "png",width = 15,height = 9,plot = p)
#
### ALL GENOMIC CLINES TOGETHER - PLOT ####
taxons <- c("Bombina","Lissotriton","Triturus","Anguis","Natrix","Sphyrapicus")

read_and_rbind_points <- function(taxons){
  taxs <- paste0("./genomic_clines/",taxons,"_GenomicClines_pointsTab.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}
read_and_rbind_intervals <- function(taxons){
  taxs <- paste0("./genomic_clines/",taxons,"_Genomic_clines_plotting_wIntervals.csv")
  dts <- lapply(taxs,read.csv,fileEncoding = "UTF-8") #%>% lapply(.,select,INDIVIDUAL_ID,CLASS,ALLELE,GENOTYPE)
  dts <- lapply(1:length(dts), function(x){mutate(dts[[x]],"TAXON" = taxons[x])})
  as.data.frame(do.call(rbind,dts))
}

dt_points <- read_and_rbind_points(taxons) %>%
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegata" ~ "Bombina bombina - B. variegata",
                                   HYBRID_ZONE == "montandoni_vulgarisIN" ~ "Lissotriton montandoni - L. vulgaris (IN)",
                                   HYBRID_ZONE == "montandoni_vulgarisOUT" ~ "Lissotriton montandoni - L. vulgaris (OUT)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilis" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "nuchalis_ruber" ~ "Sphyrapicus nuchalis - S. ruber",
                                   HYBRID_ZONE == "ruber_varius" ~ "Sphyrapicus ruber - S. varius",
                                   HYBRID_ZONE == "nuchalis_varius" ~ "Sphyrapicus nuchalis - S. varius",
                                   HYBRID_ZONE == "helvetica_natrix" ~ "Natrix helvetica - N. natrix",
                                   HYBRID_ZONE == "astreptophora_helvetica" ~ "Natrix astreptophora - N. helvetica")) %>%
  mutate("CLASS" = paste("MHC",CLASS))

dt_points$HYBRID_ZONE <- factor(dt_points$HYBRID_ZONE,levels = c("Bombina bombina - B. variegata",
                                                               "Lissotriton montandoni - L. vulgaris (IN)",
                                                               "Lissotriton montandoni - L. vulgaris (OUT)",
                                                               "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                               "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                               "Anguis colchica - A. fragilis",
                                                               "Natrix astreptophora - N. helvetica",
                                                               "Natrix helvetica - N. natrix",
                                                               "Sphyrapicus nuchalis - S. ruber",
                                                               "Sphyrapicus nuchalis - S. varius",
                                                               "Sphyrapicus ruber - S. varius"))

dt_lines <- read_and_rbind_intervals(taxons) %>% filter(INTERVAL_TYPE == "LL") %>%
  mutate("HYBRID_ZONE" = case_when(HYBRID_ZONE == "bombina_variegata" ~ "Bombina bombina - B. variegata",
                                   HYBRID_ZONE == "montandoni_vulgarisIN" ~ "Lissotriton montandoni - L. vulgaris (IN)",
                                   HYBRID_ZONE == "montandoni_vulgarisOUT" ~ "Lissotriton montandoni - L. vulgaris (OUT)",
                                   HYBRID_ZONE == "ivanbureschi_macedonicus" ~ "Triturus ivanbureschi - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_macedonicus" ~ "Triturus cristatus - T. macedonicus",
                                   HYBRID_ZONE == "cristatus_ivanbureschi" ~ "Triturus cristatus - T. ivanbureschi",
                                   HYBRID_ZONE == "anatolicus_ivanbureschi" ~ "Triturus anatolicus - T. ivanbureschi",
                                   HYBRID_ZONE == "colchica_fragilis" ~ "Anguis colchica - A. fragilis",
                                   HYBRID_ZONE == "nuchalis_ruber" ~ "Sphyrapicus nuchalis - S. ruber",
                                   HYBRID_ZONE == "ruber_varius" ~ "Sphyrapicus ruber - S. varius",
                                   HYBRID_ZONE == "nuchalis_varius" ~ "Sphyrapicus nuchalis - S. varius",
                                   HYBRID_ZONE == "helvetica_natrix" ~ "Natrix helvetica - N. natrix",
                                   HYBRID_ZONE == "astreptophora_helvetica" ~ "Natrix astreptophora - N. helvetica")) %>%
  mutate("CLASS" = paste("MHC",CLASS))

dt_lines$HYBRID_ZONE <- factor(dt_lines$HYBRID_ZONE,levels = c("Bombina bombina - B. variegata",
                                                       "Lissotriton montandoni - L. vulgaris (IN)",
                                                       "Lissotriton montandoni - L. vulgaris (OUT)",
                                                       "Triturus anatolicus - T. ivanbureschi","Triturus cristatus - T. ivanbureschi",
                                                       "Triturus cristatus - T. macedonicus","Triturus ivanbureschi - T. macedonicus",
                                                       "Anguis colchica - A. fragilis",
                                                       "Natrix astreptophora - N. helvetica",
                                                       "Natrix helvetica - N. natrix",
                                                       "Sphyrapicus nuchalis - S. ruber",
                                                       "Sphyrapicus nuchalis - S. varius",
                                                       "Sphyrapicus ruber - S. varius"))

library(ggplot2)

plot1 <- ggplot(dt_lines)+
  geom_line(aes(x = XS,y = YS,col = CLASS))+
  geom_ribbon(aes(x = XS,ymin = MIN_YS,ymax = MAX_YS,fill = CLASS),alpha = 0.3)+
  geom_point(data = dt_points,aes(x = GENOME_WIDE,y = H_INDEX,col = CLASS),alpha = 0.3)+
  scale_color_manual(values = c("MHC I" = "black","MHC II" = "green3"))+
  scale_fill_manual(values = c("MHC I" = "black","MHC II" = "green3"))+
  facet_wrap(vars(HYBRID_ZONE))+
  theme_bw()+labs(x = "GENOME-WIDE ANCESTRY",y = "MHC H-INDEX")+
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 25),
        legend.position = c(0.88,0.10),legend.title = element_blank(),legend.key.size = unit(2,units = "cm"),
        legend.text = element_text(size = 30),strip.text = element_text(size = 18,face = "italic"),
        strip.background = element_rect(fill = "white"),plot.margin = unit(c(1,1,1,1),"cm"),panel.spacing = unit(0.5,"cm"))

ggsave("plots/Drawn_genomic_clines.png",plot1,device = "png",width = 27,height = 18)

#