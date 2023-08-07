#code to fit genomic clines
#fitting genomic clines is super fast but setting their 2LL intervals is slow
#As the parameters are restricted from -1 to 1, it must be done numerically

library(dplyr)
library(tidyr)
#max likelihood functions
dt <- read.csv("input.csv",fileEncoding = "UTF-8") %>%
  mutate("GENOME_WIDE" = case_when(GENOME_WIDE == 0 ~ 0.1,GENOME_WIDE == 1 ~ 0.9, T ~ GENOME_WIDE))

normalPDF_likelihood_optim <- function(vec,p,y){ #p & y are vectors
  sd <- vec[1]
  a <- vec[2]
  b <- vec[3]
  n <- length(p)
  q <- 1-p
  #s <- sqrt(sd)
  first_part <- -(n/2)*log(2*pi)
  second_part <- -(n/2)*log(sd)
  third_part <- -(1/(2*sd))
  forth_part <- sum(((p+((2*p*q)*(a+(b*(p-q))))) - y)^2)
  #forth_part <- sum((y - (a+b*p))^2)
  (first_part+second_part+(third_part*forth_part))
}

#running all fitting inside a loop and return the parameters
hz_names <- unique(dt$HYBRID_ZONE) 
klasy <- unique(dt$CLASS)

#initiate parallelisation
library(foreach)
library(doParallel)
no_cores <- 12
cl <- makeCluster(no_cores)
registerDoParallel(cl)

excluded <- c()

#fitting the optimal cline and tracing the path to the optimal value
for(l in hz_names){
  for(m in klasy){
    print(paste(l,m))
    
    temp <- filter(dt,HYBRID_ZONE == l & CLASS == m)
    
    writeLines("",paste(l,m,"log.txt",sep = "_")) #create a file to collect tracing
    sink(paste(l,m,"log.txt",sep = "_"),append = T) #send to file
    trace(normalPDF_likelihood_optim,exit = quote(print(c(vec,returnValue()))),print = F) #trace - save each iteration
    res_optim_eq <-  try(optim(par = c(0,0,0), #optimalization
                          fn = normalPDF_likelihood_optim,
                          lower = c(0.001,-1,-1),
                          upper = c(1,1,1),
                          method="L-BFGS-B",
                          p = temp$GENOME_WIDE,y = temp$H_INDEX,
                          control = list("fnscale"=-1),
                          hessian = T))
    untrace(normalPDF_likelihood_optim) #untrace
    sink() #close connection
    
    if("try-error" %in% class(res_optim_eq)){
      excluded <- c(excluded,l)
    }
  }
}


writeLines("","log.txt") 
hz_names <- hz_names[!hz_names %in% excluded]

#here we start parallelisation
out_list <- foreach(l = hz_names, .packages = c("dplyr","tidyr")) %:% # loop through hybrid zones
  foreach(m = klasy,.packages = c("dplyr","tidyr")) %dopar% { # loop through loci
    
    
    sink("log.txt",append = T)
    print(paste(l,m))
    sink()
    
    # get some rough 2LL estimates from tracing
    LL2 <- read.table(paste(l,m,"log.txt",sep = "_"),sep = "",header = F,skip = 2,strip.white = T) %>% #read collected file
      .[,-1] %>% rename("SD" = 1,"ALPHA" = 2,"BETA" = 3,"LL" = 4) %>%
      filter(LL >= res_optim_eq$value-2 & BETA <= 1 & BETA >= -1 & ALPHA <= 1 & ALPHA >= -1 & SD <= 1 & SD >= 0.001) %>%
      select(-LL) %>% pivot_longer(cols = everything(),values_to = "VALUE",names_to = "PARAMETERS") %>%
      group_by(PARAMETERS) %>% summarise("LOWER" = min(VALUE),"UPPER" = max(VALUE))

    # add resolution to 2LL intervals (to 0.01 resolution)
    sd_vec <- seq(from = max(0,LL2$LOWER[LL2$PARAMETERS == "SD"]-0.3),to = min(1,LL2$UPPER[LL2$PARAMETERS == "SD"]+0.3),by = 0.01)
    sd_vec[sd_vec == 0] <- 0.001
    a <- seq(from = max(-1,LL2$LOWER[LL2$PARAMETERS == "ALPHA"]-0.3),to = min(1,LL2$UPPER[LL2$PARAMETERS == "ALPHA"]+0.3),by = 0.01)
    b <- seq(from = max(-1,LL2$LOWER[LL2$PARAMETERS == "BETA"]-0.3),to = min(1,LL2$UPPER[LL2$PARAMETERS == "BETA"]+0.3),by = 0.01)
    out_list1 <- list()
    for(i in sd_vec){
      for(j in a){
        for(k in b){
          out_list1[[length(out_list1)+1]] <- c(i,j,k,normalPDF_likelihood_optim(vec = c(i,j,k),p = temp$GENOME_WIDE,y = temp$H_INDEX))
        }
      }
    }
    out1 <- as.data.frame(do.call(rbind,out_list1))
    colnames(out1) <- c("SD","ALPHA","BETA","LL")
    chosen <- out1 %>% filter(LL >= (res_optim_eq$value-2)) %>%
      select(-LL) %>% pivot_longer(cols = everything(),values_to = "VALUE",names_to = "PARAMETERS") %>%
      group_by(PARAMETERS) %>% summarise("UPPER_2LL" = max(VALUE),"LOWER_2LL" = min(VALUE))
    ###

    # second stage (0.001 resolution - restricted range) - highly time consuming
    sd_vec <- seq(from = max(0,chosen$LOWER_2LL[chosen$PARAMETERS == "SD"]-0.1),to = min(1,chosen$UPPER_2LL[chosen$PARAMETERS == "SD"]+0.1),by = 0.001)
    sd_vec[sd_vec == 0] <- 0.001
    a <- seq(from = max(-1,chosen$LOWER_2LL[chosen$PARAMETERS == "ALPHA"]-0.1),to = min(1,chosen$UPPER_2LL[chosen$PARAMETERS == "ALPHA"]+0.1),by = 0.001)
    b <- seq(from = max(-1,chosen$LOWER_2LL[chosen$PARAMETERS == "BETA"]-0.1),to = min(1,chosen$UPPER_2LL[chosen$PARAMETERS == "BETA"]+0.1),by = 0.001)
    #out_list1 <- list()

    sd_min <- Inf
    sd_max <- -Inf
    a_min <- Inf
    a_max <- -Inf
    b_min <- Inf
    b_max <- -Inf

    total <- length(sd_vec)*length(a)*length(b)
    
    sink("log.txt",append = T)
    print(paste("Total",total,l,m))
    sink()
    
    counter <- 0
    for(i in sd_vec){
      for(j in a){
        for(k in b){

          counter <- counter+1
          #perc <- counter/total*100

          sink("log.txt",append = T)
          if(counter %% 10000000 == 0){
            print(paste(counter,l,m))
          }
          sink()

          ll_val <- normalPDF_likelihood_optim(vec = c(i,j,k),p = temp$GENOME_WIDE,y = temp$H_INDEX)
          if(ll_val >= (res_optim_eq$value-2)){
            if(i < sd_min){
              sd_min <- i
            }
            if(i > sd_max){
              sd_max <- i
            }
            if(j < a_min){
              a_min <- j
            }
            if(j > a_max){
              a_max <- j
            }
            if(k < b_min){
              b_min <- k
            }
            if(k > b_max){
              b_max <- k
            }
          }
          #out_list1[[length(out_list1)+1]] <- c(i,j,k,normalPDF_likelihood_optim(vec = c(i,j,k),p = temp$GENOME_WIDE,y = temp$H_INDEX))
        }
      }
    }
    chosen <- data.frame("PARAMETERS" = c("SD","ALPHA","BETA"),"UPPER_2LL" = c(sd_max,a_max,b_max),"LOWER_2LL" = c(sd_min,a_min,b_min))
    ###
    
    #changing 2LL intervals into 95% CI 
    fisher_info<-solve(-res_optim_eq$hessian)
    prop_sigma<-sqrt(diag(fisher_info))
    upper<-res_optim_eq$par+1.96*prop_sigma
    lower<-res_optim_eq$par-1.96*prop_sigma
    interval<-data.frame("HYBRID_ZONE" = l,"CLASS" = m,
                         "VALUE"=res_optim_eq$par, "UPPER_95"=upper, "LOWER_95"=lower,
                         "PARAMETERS" = c("SD","ALPHA","BETA"),"SE"=prop_sigma,"N" = nrow(temp),"LL" = res_optim_eq$value) %>%
      left_join(chosen,by = "PARAMETERS")

    saveRDS(interval,paste(l,m,"interval.RDS",sep = "_"))
    interval
} 

stopCluster(cl)
out <- as.data.frame(do.call(rbind,out_list))
write.csv(out,"GenomicClines_parameters.csv",row.names = F,fileEncoding = "UTF-8")
#

### PLOTTING GENOMIC CLINES - PREPARATION ####
# it scans within 2LL range of parameters and chooses the min and max value for each x-value
parameters <- read.csv("GenomicClines_parameters.csv",fileEncoding = "UTF-8")

cline_drawn <- function(p,a,b){
  q <- 1-p
  (p+(2*p*q)*(a+b*(p-q)))
}

xs <- seq(0,1,0.01)

no_cores <- 12
cl <- makeCluster(no_cores)
registerDoParallel(cl)

out_list <- foreach(i = unique(parameters$HYBRID_ZONE), .packages = c("dplyr","tidyr")) %:% # loop through hybrid zones
  foreach(j = unique(parameters$CLASS),.packages = c("dplyr","tidyr")) %dopar% { # loop through loci
    
    temp <- filter(parameters,HYBRID_ZONE == i & CLASS == j)
    
    intervals <- list()
    for(k in c("95","LL")){
      if(k == "95"){
        temp1 <- temp %>% select(-UPPER_2LL,-LOWER_2LL) %>% rename("UPPER" = "UPPER_95","LOWER" = "LOWER_95")
      } else {
        temp1 <- temp %>% select(-UPPER_95,-LOWER_95) %>% rename("UPPER" = "UPPER_2LL","LOWER" = "LOWER_2LL")
      }
      for(n in xs){
        
        sink("log.txt",append = T)
        if((n*100) %% 10 == 0){
          print(n)
        }
        sink()
        
        min_ys <- Inf
        max_ys <- -Inf
        for(l in seq(temp1$LOWER[temp1$PARAMETERS == "ALPHA"],temp1$UPPER[temp1$PARAMETERS == "ALPHA"],0.001)){
          for(m in seq(temp1$LOWER[temp1$PARAMETERS == "BETA"],temp1$UPPER[temp1$PARAMETERS == "BETA"],0.001)){
            val_ys <- cline_drawn(p = n,a = l,b = m)
            if(val_ys > max_ys){
              max_ys <- val_ys
            }
            if(val_ys < min_ys){
              min_ys <- val_ys
            }
          }
        }
        intervals[[length(intervals)+1]] <- c(k,n,min_ys,max_ys)
      }
    }
    intervals <- as.data.frame(do.call(rbind,intervals))
    colnames(intervals) <- c("INTERVAL_TYPE","XS","MIN_YS","MAX_YS")
    
    data.frame("HYBRID_ZONE" = i,"CLASS" = j,"XS" = as.character(xs),"YS" = cline_drawn(xs,a = temp$VALUE[temp$PARAMETERS == "ALPHA"],
                                                                                        b = temp$VALUE[temp$PARAMETERS == "BETA"])) %>%
      left_join(intervals,by = "XS") %>%
      mutate("ALPHA" = temp1$VALUE[temp$PARAMETERS == "ALPHA"],"BETA" = temp$VALUE[temp1$PARAMETERS == "BETA"])
}
stopCluster()

out <- as.data.frame(do.call(rbind,out_list))
write.csv(out,"Genomic_clines_plotting_wIntervals.csv",row.names = F,fileEncoding = "UTF-8")
#