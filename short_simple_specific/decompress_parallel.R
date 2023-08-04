#Decompress multiple files in parallel

library(foreach)
library(doParallel)

gzs <- list.files(pattern = ".gz")

no_cores <- length(gzs) #number of cores to use
cl <- makeCluster(no_cores)
registerDoParallel(cl)

foreach(i = gzs) %dopar% {
  system(paste("gzip -d",i))
}

stopCluster(cl)