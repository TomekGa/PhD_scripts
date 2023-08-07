#Tomek Gaczorek
#tomek.gaczorek@gmail.com

#wrapper to fit geographic clines with hzar package using multiple cores
#it also automatically gathers the best model and all needed paramethers

library(hzar)
library(dplyr)
library(foreach)

beg_time <- Sys.time()
### FUNCTIONS ####
substraction_vector <- function(x,minus=2,times){
  out_vec <- c()
  for(i in 1:times){
    if(i == 1){
      out_vec <- x
    } else{
      out_vec <- c(out_vec,out_vec[length(out_vec)]-minus)
    }
  }
  sort(out_vec)
}

list.append <- function(in_list,new){
  in_list[[length(in_list)+1]] <- new
  in_list
}

with_timeout <- function(expr, cpu, elapsed){
  expr <- substitute(expr)
  envir <- parent.frame()
  setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  eval(expr, envir = envir)
}
#########################

#REQUIRES:
args <- c("Triturus_hzar_input.csv",3,3,1) #ADJUST
# 1. path to input file (check exemplary_data folder for the proper format of input table)
# 2. number of loci (be aware that they should be at the end of table followed by 'number of samples')
# 3. number of column with distance values
# 4. number of column with unique values for different hybrid zones

## decide about type of the model here #ADJUST
plotting <- T
scaling_args <- c('none','fixed','free')
tails_args <- c('none','left','right','mirror','both')
models_args <- list()
cwd <- getwd()
## or alternalively here
#models_args <- list(c('none','none'),c('none','left'))

cwd <- getwd()
markers <- c("MHCI","MHCII","genomic")

#OUTPUTS:
# - Logs
# - Plots of observed data
# - Tracing plots (fitting)
# - Plot of all clines fitted
# - Plot of maximum likelihood cline
# - table with 2LL parameters
# - RDS file with gathered parameters
# - RDS files for individual hybrid zones (not included in "exemplary_data" folder due to big size)

#RUNNING
# Rscript hzar_clines_analysis_parallel.R

### ### ### ### ###

## Read file for analysis
in_data <- read.csv(args[1],encoding = "UTF-8")

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528))


#registerDoSEQ()

#### EVIL LOOP ####

# prepare values for looping
hybrid_zones <- sort(as.character(unique(in_data[,as.integer(args[4])])),decreasing = T)
#hybrid_zones <- "ivanbureschi_macedonicus"
loci_col_numbers <- substraction_vector(x = ncol(in_data)-1,times = as.integer(args[2]))
sample_size_col_numbers <- substraction_vector(x = ncol(in_data),times = as.integer(args[2]))

#CAL <- list()
dir.create('obsData')
dir.create('Trace_plots')
dir.create('Compared_clines')
dir.create('ML_clines')
dir.create('CALs')
dir.create('Logs')

library(foreach)
library(doParallel)
#no_cores <- 10
no_cores <- length(hybrid_zones)*as.integer(args[2])
cl <- makeCluster(no_cores)
registerDoParallel(cl)

foreach(i = hybrid_zones, .packages = c("dplyr","hzar")) %:% # loop through hybrid zones
            foreach(j = c(1:as.integer(args[2])),.packages = c("dplyr","hzar")) %dopar% { # loop through loci
              
              writeLines("",paste(cwd,"/Logs/",i,"_",j,".txt",sep = ""))
              sink(paste(cwd,"/Logs/",i,"_",j,".txt",sep = ""),append = T)
              ###create all neccesary lists
              #lista_smierci <- list()
              
              ###load data
              in_data_hz <- filter(in_data,in_data[,as.integer(args[4])] == i)
              if(length(na.omit(in_data_hz[,loci_col_numbers[j]])) != 0){
                
                #create empty lists
                CAL <- list()
                CAL$obs <- list()
                CAL$models <- list()
                CAL$fitRs <- list()
                CAL$runs <- list()
                CAL$analysis <- list()
                CAL$to_check <- list()
                
                ###
                #in_data_hz$distance <- in_data_hz$distance*10
                #in_data_hz$percantage.gen <- round(in_data_hz$percantage.gen,3)
                ### 
                CAL$obs <- hzar.doMolecularData1DPops(in_data_hz[,as.integer(args[3])],
                                                      in_data_hz[,loci_col_numbers[j]]/100,
                                                      in_data_hz[,sample_size_col_numbers[j]])
                
                ###plot observed data
                if(plotting == T){
                  setwd(paste(cwd,'obsData',sep = '/'))
                  bitmap(width=900, height=900, res=200,
                         file=paste('obsData',i,j,'.png',sep = '_'),pointsize=8,units="px")
                  hzar.plot.obsData(CAL$obs)
                  dev.off()
                  setwd(cwd)
                }
                
                ###calculates preliminary parameters based on equations
                if(length(models_args)==0){
                  for(k in scaling_args){
                    for(l in tails_args){
                      id <- paste(k,l,sep = '_')
                      CAL$models[[id]] <- hzar.makeCline1DFreq(CAL$obs,k,l)
                    }
                  }
                } else {
                  for(k in 1:length(models_args)){
                    id <- paste(models_args[[k]],collapse = '_')
                    CAL$models[[id]] <- hzar.makeCline1DFreq(CAL$obs,
                                                             models_args[[k]][1],
                                                             models_args[[k]][2])
                  }
                }
                
                ### focus on given distance from hz
                min_dis <- min(in_data_hz[,as.integer(args[3])])
                max_dis <- max(in_data_hz[,as.integer(args[3])])
                # CAL$models <- lapply(CAL$models,hzar.model.addBoxReq,
                #                      min_dis-30,max_dis+30)
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_0.RDS",sep = "")) #checkpoint 0
                
                # Compile each of the models to prepare for fitting
                #giving 100 chances to fit model
                #important step - initial generation of data
                g_names <- names(CAL$models)
                for(k in 1:length(CAL$models)){
                  print(k)
                  counter <- 0
                  checker <- F
                  while(checker == F){
                    counter <- counter + 1
                    print(paste("counter:",counter))
                    outcome <- hzar.first.fitRequest.old.ML(CAL$models[[k]],obsData = CAL$obs,verbose=FALSE)
                    try_var <- try(chol(outcome$cM))
                    try_doFit <- capture.output(with_timeout(hzar.doFit(outcome),cpu = 3000, elapsed = 3000))
                    if(((class(try_var) != "try-error") && (grepl(" successfully",try_doFit[length(try_doFit)])==T))){
                      CAL$fitRs$init[[g_names[k]]] <- outcome
                      #lista_smierci[[k]] <- outcome
                      checker <- T
                    } else if(counter > 100) {
                      checker <- T
                    }
                  }
                }
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_1.RDS",sep = "")) #checkpoint 1
                #CAL$fitRs$init <- lapply(CAL$models, hzar.first.fitRequest.old.ML, obsData=CAL$obs, verbose=FALSE)
                
                
                CAL$fitRs$init <- lapply(CAL$fitRs$init, function(mdl) {
                  mdl$mcmcParam$chainLength <- chainLength #1e5
                  mdl$mcmcParam$burnin <- chainLength %/% 10
                  mdl})
                
                ### seems that it creates independend seeds for each model
                CAL$fitRs$init <- hzar.multiFitRequest(CAL$fitRs$init,
                                                       rotateSeed=TRUE,skip=50,
                                                       baseChannel=NULL,each=1,
                                                       baseSeed=c(596,528,124,978,544,99))
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_2.RDS",sep = "")) #checkpoint 2
                
                ### run initial chain
                CAL$runs$init <- list()
                
                
                # giving 10 chances to fit model
                g_names <- names(CAL$fitRs$init)
                for(k in 1:length(CAL$fitRs$init)){
                  counter <- 0
                  checker <- F
                  while(checker == F){
                    counter <- counter + 1
                    try_var <- try(hzar.doFit(CAL$fitRs$init[[k]]))
                    if((length(try_var) > 4) | (counter > 10)){
                      CAL$runs$init[[g_names[k]]] <- try_var
                      checker <- T
                    }
                  }
                }
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_3.RDS",sep = "")) #checkpoint 3
                
                ### ploting trace
                if(plotting == T){
                  setwd(paste(cwd,'Trace_plots',sep = '/'))
                  for(k in 1:length(CAL$runs$init)){
                    bitmap(width=900, height=900, res=200,
                           file=paste('Trace',i,j,names(CAL$runs$init)[k],
                                      '.png',sep = '_'),pointsize=8,units = "px");
                    plot(hzar.mcmc.bindLL(CAL$runs$init[[k]]));
                    dev.off()
                  }
                  setwd(cwd)
                }
                
                ## Compile a new set of fit requests using the initial chains 
                CAL$fitRs$chains <- lapply(CAL$runs$init,hzar.next.fitRequest)
                
                ## Replicate each fit request 3 times, keeping the original
                ## seeds while switching to a new seed channel.
                CAL$fitRs$chains <- hzar.multiFitRequest(CAL$fitRs$chains,
                                                         each=3,baseSeed=NULL)
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_4.RDS",sep = "")) #checkpoint 4
                ### randomize initial fit values and run a chain of 3 runs for every fit request
                ### not all values result in a proper model so it is repeated until 
                for(k in 1:length(CAL$fitRs$chains)){
                  names_params <- names(CAL$fitRs$chains[[k]]$modelParam$init)
                  checker <- T
                  while(checker == T){
                    for(l in names_params){
                      if(l == 'center'){
                        CAL$fitRs$chains[[k]]$modelParam$init[[l]] <- runif(1,min_dis-30,max_dis+30)
                      } else if(l == 'width' | startsWith(l,prefix = "delta")){
                        CAL$fitRs$chains[[k]]$modelParam$init[[l]] <- runif(1,0,sum(abs(c(min_dis-30,max_dis+30))))
                      } else {CAL$fitRs$chains[[k]]$modelParam$init[[l]] <- runif(1,0,1)}
                    }
                    print(paste("order",k))
                    try_val <- tryCatch(hzar.chain.doSeq(CAL$fitRs$chains[[k]]),warning = function(x){NA})
                    if(is.na(try_val) != T){
                      CAL$runs$chains[[k]] <- try_val
                      checker <- F
                    }
                  }
                }
                
                ## Go ahead and run a chain of 3 runs for every fit request
                # CAL$runs$chains <-  hzar.doChain.multi(CAL$fitRs$chains,
                #                                        doPar=F,inOrder=T,
                #                                        count=3)
                # for(k in 1:length(CAL$fitRs$chains)){
                #   print(paste("order",k))
                #   try_val <- tryCatch(hzar.chain.doSeq(CAL$fitRs$chains[[k]]),warning = function(x){NA})
                #   CAL$runs$chains[[k]] <- hzar.chain.doSeq(CAL$fitRs$chains[[k]])
                # }
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_5.RDS",sep = "")) #checkpoint 5
                
                ### models convergence - important to check
                ## find it later on in panel analysis$mdl_convergence 
                summary_list <- list()
                for(k in seq(1,length(CAL$runs$chains),by=3)){
                  x <- summary(do.call(mcmc.list, lapply(CAL$runs$chains[k:k+2],function(x){hzar.mcmc.bindLL(x[[3]])})))
                  x_stats <- cbind(x$statistics,rep((length(summary_list)+1),times = nrow(x$statistics)))
                  summary_list <- list.append(summary_list,x_stats)
                }
                summary_list <- do.call(rbind,summary_list)
                colnames(summary_list)[5] <- "model_num"
                summary_list <- as.data.frame(summary_list)
                summary_list$parameter <- row.names(summary_list)
                row.names(summary_list) <- c(1:nrow(summary_list))
                CAL$to_check$mdl_convergence <- summary_list
                
                ## Create a model data group (hzar.dataGroup object) for each
                ## model from the initial runs.
                CAL$analysis$initDGs <- lapply(CAL$runs$init,hzar.dataGroup.add)
                
                ## Create a model data group for the null model (expected allele
                ## frequency independent of distance along cline) to include in analysis.
                CAL$analysis$initDGs$nullModel <- hzar.dataGroup.null(CAL$obs)
                
                ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
                CAL$analysis$oDG <- hzar.make.obsDataGroup(CAL$analysis$initDGs)
                #copy labels
                CAL$analysis$oDG <- hzar.copyModelLabels(CAL$analysis$initDGs,CAL$analysis$oDG)
                
                ## Convert all runs to hzar.dataGroup objects, adding them to
                ## the hzar.obsDataGroup object.
                CAL$analysis$oDG <- hzar.make.obsDataGroup(lapply(CAL$runs$chains,
                                                                  hzar.dataGroup.add),
                                                           CAL$analysis$oDG);
                
                ##compare models graphically
                if(plotting == T){
                  setwd(paste(cwd,'Compared_clines',sep = '/'))
                  bitmap(width=900, height=900, res=200,
                         file=paste('clines',i,j,'.png',sep = '_'),pointsize=8,units = "px")
                  hzar.plot.cline(CAL$analysis$oDG)
                  dev.off()
                  setwd(cwd)
                }
                
                ### AIC table
                CAL$to_check$AICcTable <- hzar.AICc.hzar.obsDataGroup(CAL$analysis$oDG)
                
                ### model with minimum AIC score
                CAL$to_check$lo_AIC_model_name <- rownames(CAL$to_check$AICcTable)[which.min(CAL$to_check$AICcTable$AICc)]
                ### dataGroup object for above model
                CAL$to_check$lo_AIC_DG <- CAL$analysis$oDG$data.groups[[CAL$to_check$lo_AIC_model_name]]
                
                #saveRDS(CAL,file = paste(cwd,"/Logs/",i,"_",j,"_6.RDS",sep = "")) #checkpoint 6
                
                if(CAL$to_check$lo_AIC_model_name != "nullModel"){ #clines for the null model does not make any sense and it doesn't work
                  ## Look at the variation in parameters for the selected model
                  CAL$to_check$lo_AIC_parameters_variation <- 
                    hzar.getLLCutParam(CAL$to_check$lo_AIC_DG,
                                       names(CAL$to_check$lo_AIC_DG$data.param))
                  ## Print the maximum likelihood cline for the selected model
                  CAL$to_check$lo_AIC_ML_cline_parameters <- 
                    hzar.get.ML.cline(CAL$to_check$lo_AIC_DG)
                  
                  ## Plot the maximum likelihood cline for the selected model
                  if(plotting == T){
                    setwd(paste(cwd,'ML_clines',sep = '/'))
                    bitmap(width=900, height=900, res=200,
                           file=paste('ML_cline',i,j,'.png',sep = '_'),pointsize=8,units = "px")
                    hzar.plot.fzCline(CAL$to_check$lo_AIC_DG,
                                      main = CAL$to_check$lo_AIC_model_name)
                    dev.off()
                    setwd(cwd)
                  }
                } else {
                  writeLines("",paste(cwd,'/ML_clines/IamNull',"_",i,"_",j,'.txt',sep = ''))
                }
                setwd(paste(cwd,'CALs',sep = '/'))
                file_name <- paste(i,"locus",j,".RDS",sep = "")
                saveRDS(CAL,file = file_name)
                setwd(cwd)
              }
              sink()
            }
stopCluster(cl)

setwd(paste(cwd,'CALs',sep = '/'))
files_list <- list.files()
loci <- c(1:as.integer(args[2]))
CAL <- list()
for(i in hybrid_zones){
  CAL[[i]] <- list()
  for(j in loci){
    file_name <- paste(i,"locus",j,".RDS",sep = "")
    if(file_name %in% files_list){
      CAL[[i]][[j]] <- readRDS(file = file_name)
      file.remove(file_name)
    } else {
      CAL[[i]][[j]] <- list()
    }
  }
}

setwd(cwd)
saveRDS(CAL,"CAL.RDS")

### Taking out parameters for plotting in a flow ###

hz_names <- names(CAL)
new_list <- list()
LL_list <- list()
print("Merging parameters")
for(i in hz_names){
  new_list[[i]] <- list()
  for(j in 1:length(markers)){
    print(paste(i,j))
    if(length(CAL[[i]][[j]])==0 | CAL[[i]][[j]]$to_check$lo_AIC_model_name == "nullModel"){
      new_list[[i]][[markers[j]]]$params_all <- NA
      new_list[[i]][[markers[j]]]$LL_variance <- NA
      new_list[[i]][[markers[j]]]$function_cline <- NA
      new_list[[i]][[markers[j]]]$obs <- NA
      new_list[[i]][[markers[j]]]$dist_xs <- NA
      new_list[[i]][[markers[j]]]$CI_data <- NA
    } else {
      new_list[[i]][[markers[j]]]$params_all <- CAL[[i]][[j]]$to_check$lo_AIC_ML_cline_parameters$param.all
      new_list[[i]][[markers[j]]]$LL_variance <- CAL[[i]][[j]]$to_check$lo_AIC_parameters_variation
      new_list[[i]][[markers[j]]]$function_cline <- CAL[[i]][[j]]$to_check$lo_AIC_ML_cline_parameters$clineFunc
      new_list[[i]][[markers[j]]]$obs <- CAL[[i]][[j]]$obs$frame
      temp <- hzar.getCredParamRed(CAL[[i]][[j]]$to_check$lo_AIC_DG)
      new_list[[i]][[markers[j]]]$dist_xs <- seq(min(new_list[[i]][[markers[j]]]$obs$dist),max(new_list[[i]][[markers[j]]]$obs$dist),1) #every 1 unit of distance, modify if needed
      new_list[[i]][[markers[j]]]$CI_data <- temp$fzCline(new_list[[i]][[markers[j]]]$dist_xs)
      #LL table
      VALUES <- c(unname(unlist(new_list[[i]][[markers[j]]]$params_all)),unname(unlist(new_list[[i]][[markers[j]]]$LL_variance)))
      PARAMETERS <- c(names(new_list[[i]][[markers[j]]]$params_all),names(new_list[[i]][[markers[j]]]$LL_variance))
      HYBRID_ZONE <- rep(i,times = length(VALUES))
      MARKER <- rep(markers[j],times = length(VALUES))
      out_table <- data.frame(HYBRID_ZONE,MARKER,PARAMETERS,VALUES) %>%
        mutate("LEVEL" = case_when(grepl("2LL",PARAMETERS) == F ~ "main",
                                   grepl("2LLLow",PARAMETERS) == F ~ "low",
                                   grepl("2LLHigh",PARAMETERS) == F ~ "high")) %>%
        mutate("PARAMETERS" = stringr::str_replace(PARAMETERS,"2LLLow|2LLHigh",""))
      LL_list[[length(LL_list)+1]] <- out_table
    }
  }
}

saveRDS(new_list,"parameters.RDS")
write.csv(as.data.frame(do.call(rbind,LL_list)),"LL_parameters.csv",row.names = F,fileEncoding = "UTF-8")

#####
print(Sys.time()-beg_time)


