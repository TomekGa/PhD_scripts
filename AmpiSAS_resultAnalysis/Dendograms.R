#Generates dendogram for allele sequences or individual genotypes based on binary matrix

#REQUIRES
#aligned fasta with allele sequences or binary matrix of MHC genotypes
#csv file with species information for each individual or some other information for each allele

#OUTPUT
#dendogram showing similarity and potential clustering of individuals or MHC alleles

library(ggplot2)
library(ggdendro)
library(dplyr)

#for sequences
alligned <- seqinr::read.alignment("4_Anguis_MHCII_forth_om_lib1_150621_alligned.fasta",format = "fasta") #ADJUST
distance <- seqinr::dist.alignment(alligned)

#for binary matrices
dt <- read.csv("4_Anguis_MHCII_forth_om_lib1_150621_genotypes.csv",row.names = 1) #ADJUST
distance <- dist(dt,method = "binary")

p <- hclust(distance,method = "average")

p <- as.dendrogram(p)

g_data <- dendro_data(p,type = "rectangle")

sps <- read.csv("Anguis_ids.csv") #ADJUST
sps$ID <- as.factor(as.character(sps$ID))

g_data$labels <- left_join(g_data$labels,sps,by = c("label" = "ID"))
g_data$labels <- select(g_data$labels,-y)
g_data$segments <- left_join(g_data$segments,g_data$labels,by = "x")


g <- ggplot(g_data$segments) + 
  geom_segment(aes(x,y, xend = xend,yend = yend))+
  geom_segment(aes(x, y, xend = xend,yend = yend,col = SPECIES))+
  geom_text(aes(x, y=-0, label = label), size = 1,angle = 90)+
  theme_bw()+theme(axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.x = element_blank(),
                   axis.text.x = element_blank())+
  ylab("Distance")

ggsave("Dendogram.png",device = "png",plot = g,width = 20)
