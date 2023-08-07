#STRUCTURE ACCEPTS ONLY LINUX'S LINE ENDINGS

#perform STRUCTURE analysis in parallel

#REQUIRES
#input table for STRUCTURE (see STRUCTURE documentation)
#driving files for structure software (see STRUCTURE documentation)

#OUTPUT
#standard STRUCTURE output

library(foreach)
library(doParallel)
no_cores <- 20
taxon <- "Triturus"
cl <- makeCluster(no_cores)
registerDoParallel(cl)

inputs <- list.files(path = paste0("./",taxon,"/inputs"),include.dirs = F)
ks <- c(2)
reps <- c(1:5)
seeds <- sample.int(10000,length(klasy)*length(ks)*length(reps))
dir.create("Structure_outputs")

writeLines("","log.txt")

foreach(i = klasy,.packages = c("dplyr")) %:% #loop over inputs
  foreach(k = ks,.packages = c("dplyr")) %:% #loop over k-values
  foreach(j = reps,.packages = c("dplyr")) %dopar% { #loop over replicates
    
    system(paste0("dos2unix ",i))
    
    sink("log.txt",append = T)
    print(paste(i,k,j))
    sink()
    
    core <- strsplit(i,split = "\\.|\\/")[[1]][4]
    out_file_name <- paste0("./",taxon,"/Structure_outputs/Structure_",core,"_k",k,"_r",j,".txt")
    in_temp_file <- read.table(i,sep = "\t",fileEncoding = "UTF-8")
    
    ind_num <- (nrow(in_temp_file)-1)/2
    num_loci <- ncol(in_temp_file)-3
    i_index <- which(klasy == i)  
    
    which_seed <- ((i_index-1)*length(reps))+j #use consecutive seeds generated before the loop
    
    cmd <- paste("structure -i",i,"-m Structure_mainparams -e Structure_extraparams -L",num_loci,"-N",ind_num,
                 "-D",seeds[which_seed],"-K",k,"-o",out_file_name)
    system(cmd)
  }
stopCluster(cl)
