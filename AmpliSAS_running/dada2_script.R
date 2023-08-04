####### Variables ###########

#REMEMBER ABOUT FILE WITH IDs (named <species>_ids.csv) - 1st column with IDs

species <- "Anguis"
MHC_class <- "I"
pr_F <- "TYVCRCTYCHWKYDSHMK"
pr_R <- "NYVYDRYRMNCASAVYRR" # not a reverse complement
reads_path <- "/mnt/matrix/projects/backup_miseq/210122_M01530_0066_000000000-D9DM3/Data/Intensities/BaseCalls"

#############################
#Tomek Gaczorek
#tomek.gaczorek@gmail.com
#script using dada2 genotyping software as the alternative to AmpliSAS

# Requirement - cutadapt
library(dada2)

######## FUNCTIONS ###########
copy_reads <- function(ID_file,MHC_class,in_path){
  ids <- read.csv(ID_file)$individual_id
  cmd <- paste("cp ",in_path,"/*",ids,"*-",MHC_class,"_*R*.fastq.gz ./",sep = "")
  sapply(cmd,system)
}

trim_primers <- function(vfiles_f,vfiles_r,primer_f,primer_r,rev_comp = F){
  if(rev_comp == F){
    primer_r <- stringi::stri_reverse(chartr("ACGTMRWSYKVHDBN","TGCAKYWSRMBDHVN",primer_r))
  }
  for(i in 1:length(vfiles_f)){
    cmd <- paste("/home/tomek/miniconda3/bin/cutadapt -e 0 -g",primer_f,"-G",primer_r,"-o out_1.fastq",
                 "-p out_2.fastq",vfiles_f[i],vfiles_r[i]) #local
    # cmd <- paste("cutadapt -e 0 -g",primer_f,"-G",primer_r,"-o out_1.fastq",
    #              "-p out_2.fastq",vfiles_f[i],vfiles_r[i])
    system(cmd)
    file.rename("out_1.fastq",vfiles_f[i])
    file.rename("out_2.fastq",vfiles_r[i])
  }
}

##############################
#proper destination
dir.create(paste(species,"_MHC",MHC_class,sep = ""))
setwd(paste(species,"_MHC",MHC_class,sep = ""))

#copying proper files
copy_reads(paste("../",species,"_ids.csv",sep = ""),MHC_class,reads_path)

#unzip
system("gzip -d *.gz")

#list reads files
fls_f <- sort(list.files(pattern = "_R1_"))
fls_r <- sort(list.files(pattern = "_R2_"))

#trim primers with cutadapt - overrides previous files
trim_primers(fls_f,fls_r,pr_F,pr_R)

#quality profiles - randomly chosen 9
png("./ForwardQprofile.png",type = "cairo",width = 800,height = 700)
plotQualityProfile(sample(fls_f,size = pmin(9,length(fls_f))))
dev.off()
png("./ReverseQprofile.png",type = "cairo",width = 800,height = 700)
plotQualityProfile(sample(fls_r,size = pmin(9,length(fls_r))))
dev.off()

#filtering
sample_names <- sapply(strsplit(basename(fls_f), "_"), `[`, 1)
filtFs <- file.path("./filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path("./filtered", paste0(sample_names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fls_f,filtFs, fls_r, filtRs,maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=T, multithread=TRUE) # On Windows set multithread=FALSE

#learn error rates
errF <- learnErrors(filtFs, multithread=TRUE,randomize = T,MAX_CONSIST = 30)
errR <- learnErrors(filtRs, multithread=TRUE,randomize = T,MAX_CONSIST = 30)

#error plots
png("./ErrorsF.png",type = "cairo",width = 800,height = 700)
plotErrors(errF, nominalQ=TRUE)
dev.off()
png("./ErrorsR.png",type = "cairo",width = 800,height = 700)
plotErrors(errR, nominalQ=TRUE)
dev.off()

#inference (pooling can be added for low frequency alleles)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE,MIN_ABUNDANCE = 1,selfConsist = T,pool = T) #min abundace can be ignored
dadaRs <- dada(filtRs, err=errR, multithread=TRUE,MIN_ABUNDANCE = 1,selfConsist = T,pool = T) #in the number of reads
# ^ selfConsist - more attention to error model generated
# ^ pool = consider all samples together (generates a lot of noice)

#merge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE,minOverlap = 5)
names(mergers) <- sample_names

#ASV
seqtab <- makeSequenceTable(mergers) #sequences in columns, individuals in rows

#remove chimeras
seqtab.nochim <- t(removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE))
sum(seqtab.nochim)/sum(seqtab) #percent of non-chimera sequences

#number of reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(t(seqtab.nochim)))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
write.csv(track,"track.csv",fileEncoding = "UTF-8")

#changing format to sequences in rows and individuals in columns
dt <- cbind("SEQ" = rownames(seqtab.nochim), data.frame(seqtab.nochim, row.names=NULL,check.names = F))
write.csv(dt,"dada2_output.csv",row.names = F,fileEncoding = "UTF-8")

#removal
unlink("./filtered",recursive = T)
file.remove(c(fls_f,fls_r))

