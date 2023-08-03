#Script to plot the logit tranformation of PAFs to choose a better PAF threshold for filtering
#Note it requires the prior analysis of ampliSAS output with "Amplisas_output_analysis.R" script (you can initially set PAF threshold to 0%)

source("~/Functions_to_be_loaded.R")
library(ggplot2)
library(patchwork)

analyzed_file <- "1_first_Myodes_MHCI_110222_om.xlsx" #ADJUST - original AmpliSAS output
mod <- paste(strsplit(analyzed_file,split = "\\.")[[1]][1],"modified.csv",sep = "_") #output from "Amplisas_output_analysis.R" script
appx <- paste(strsplit(analyzed_file,split = "\\.")[[1]][1],"appendix.csv",sep = "_") #output from "Amplisas_output_analysis.R" script

data <- read.csv(mod,check.names = F,stringsAsFactors = F)[,-c(1:3)]
appxes <- read.csv(appx,check.names = F,stringsAsFactors = F)

list_data <- adding_freqs_info(data,appxes) #recalculating allele frequencies
fre <- list_data[[2]]

max_freq <- sort(fre$MAX_FREQ,decreasing = T)

logit_max_freq <- log((max_freq)/(1-(max_freq))) #logit transformation

xs <- c(1:length(logit_max_freq)) #ordering of alleles

p1 <- ggplot()+
  geom_point(aes(x = xs,y = logit_max_freq))+
  theme_bw()+
  labs(y = "Logit max frequency")+
  scale_y_continuous(sec.axis = sec_axis(trans =~(1/(1+exp(-.))*100),
                                         breaks = seq(round(min(max_freq),4),round(max(max_freq),4),0.01)*100,
                                         name = "Max frequency [%]"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 20))

p2 <- ggplot()+
  geom_col(aes(x = xs,y = max_freq*100),fill = "black",color = "black")+
  theme_bw()+
  labs(y= "Max frequency [%]")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 20))

p_combined <- p1+p2

ggsave("logit_PAF.png",p_combined,device = "png",width = 18,height = 9)