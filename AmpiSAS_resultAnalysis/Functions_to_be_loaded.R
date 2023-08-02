#Functions needed for the automatic analysis of ampliSAS results
#Tomek Gaczorek
#tomek.gaczorek@gmail.com

###### OPERATORS ##########

`%s%` <- function(x,y){
  #syntax - "x" %s% "y"
  cmd <- paste(x,"[",x,y,"]",sep = "")
  eval(parse(text = cmd))
}

`%sr%` <- function(x,y){
  #syntax - "x" %s% "y"
  cmd <- paste(x,"[",x,y,"]",sep = "")
  out <- eval(parse(text = cmd))
  assign(x,out,envir = .GlobalEnv)
}

`%src%` <- function(x,y,z,env = .GlobalEnv){
  #syntax - '%src%'(x,y,z)
  cmd <- paste(x,"[",x,y,"] <- ",z,sep = "")
  eval(parse(text = cmd),envir = env)
}

#### FUNCTIONS #### 

#translate protein but keep Gaps
translate_wGaps <- function(x,frame_v){
  out <- seqinr::translate(x,frame = frame_v)
  x_indexes <- which(out == "X")
  x_indexes <- x_indexes[!(x_indexes %in% c(1,length(out)))]
  if(length(x_indexes > 0)){
    for(i in x_indexes){
      indexes_DNA <- c((i*3-2+frame_v):(i*3+frame_v))
      if(sum(!(x[indexes_DNA] == "-")) == 0){
        out[i] <- "-"
      }
    }
  }
  out
}

# Easily add row to any chosen place
add_row2 <- function(df,row,place){
  if(place > 1 & place < nrow(df)){
    out <- rbind(df[1:(place-1),],row,df[place:nrow(df),])
  } else if(place == 1){
    out <- rbind(row,df)
  } else if(place == nrow(df)){
    out <- rbind(df,row)
  }
  row.names(out) <- NULL
  out
}

# Used for additional AmpliSAS data analysis
read_and_adjust_II <- function(file_name,split_var,pattern_var){
  out_list <- list()
  raw <- readxl::read_xlsx(file_name,col_names = F)
  head_3 <- raw[1:3,9:ncol(raw)]
  out_list[[1]] <- readxl::read_xlsx(file_name,skip = 3)[,-c(6:8)]
  
  col_nam <- strsplit(colnames(out_list[[1]])[-(1:6)],split_var)
  indexes <- lapply(col_nam,grep,pattern = pattern_var)
  new_col_nam <- character(length = length(col_nam))
  for(i in 1:length(col_nam)){
    new_col_nam[i] <- col_nam[[i]][indexes[[i]]]
  }
  
  colnames(out_list[[1]])[-(1:5)] <- c("NAME",new_col_nam)
  colnames(head_3) <- c("DESCRIPTION",colnames(out_list[[1]])[7:ncol(out_list[[1]])])
  out_list[[2]] <- head_3
  out_list
}

read_and_adjust <- function(file_name){
  out_list <- list()
  raw <- readxl::read_xlsx(file_name,col_names = F)[,-c(6:8)]
  head_3 <- raw[1:3,6:ncol(raw)]
  out_list[[1]] <- readxl::read_xlsx(file_name,skip = 3)[,-c(6:8)]
  colnames(out_list[[1]])[6] <- "NAME"
  colnames(head_3) <- c("DESCRIPTION",colnames(out_list[[1]])[7:ncol(out_list[[1]])])
  out_list[[2]] <- head_3
  out_list
}

repeatability <- function(dt,log){
  duplicated_names <- colnames(dt)[grep("d",colnames(dt))]
  counterpart_names <- sapply(sapply(duplicated_names,strsplit,"d"),"[",1)
  ###
  counterpart_names <- counterpart_names[counterpart_names %in% colnames(dt)]
  duplicated_names <- duplicated_names[duplicated_names %in% paste(counterpart_names,"d",sep = "")]
  ###
  all_names <- sort(c(counterpart_names,duplicated_names))
  temp <- dt[,all_names]
  #print(quantile(na.omit(unlist(temp)),probs = seq(0,1,0.1)))
  #View(temp)
  write.csv(temp,"repeatability.csv",row.names = F,fileEncoding = "UTF-8")
  #threshold <- as.integer(readline(prompt = "Minimal acceptable coverage of an allele: "))
  #temp[temp < threshold] <- NA
  write(paste("REPEATABILITY\n\nID\tID_duplicate\trepeatability\n"),log,append = T)
  out_lines <- c()
  reps <- c()
  for(i in 1:ncol(temp)){
    if(i %% 2 != 0){
      indexes_1 <- which(!is.na(temp[,i]))
      indexes_2 <- which(!is.na(temp[,i+1]))
      repeat_value <- round(sum(indexes_1 %in% indexes_2)/length(unique(c(indexes_1,indexes_2))),2)
      reps[length(reps)+1] <- repeat_value
      out_lines <- c(out_lines,paste(colnames(temp)[i],colnames(temp)[i+1],repeat_value,"\n",sep = "\t"))
    }
  }
  out_lines <- c(out_lines, paste("MEAN",round(mean(reps),2),"\n",sep = "\t"))
  out_lines <- c(out_lines,paste(rep("*",times = 100),collapse = ""))
  write(out_lines,log,append = T)
}

filter_by_length <- function(dt,desired_len,codons_range,rem_shifts = F,log){
  minimum <- min(desired_len) - (3*codons_range)
  maximum <- max(desired_len) + (3*codons_range)
  out <- filter(dt,LENGTH >= minimum & LENGTH <= maximum)
  shifts <- "OFF"
  
  if(rem_shifts == T){
    out_nrow <- nrow(out)
    out <- filter(out,(LENGTH-desired_len) %% 3 == 0)
    shifts <- out_nrow - nrow(out)
  }
  ###checking empty inds
  alleles_inds <- apply(out[,7:ncol(out)],2,function(x){sum(!is.na(x))})
  if(0 %in% alleles_inds){
    index <- which(alleles_inds == 0)
    inds_names <- names(alleles_inds)[index]
    removed_num <- length(index)
    out <- out[,-(index+6)]
  } else{
    inds_names <- c()
    removed_num <- 0
  }
  n_o <- paste(inds_names,collapse = " ")
  ###
  number_filtered <- nrow(dt)-nrow(out)
  ###
  counts <- freqs(dt$LENGTH)
  occur <- paste("(",counts,")",sep = "")
  c_o <- paste(paste(names(counts),occur),collapse = " ")
  ###
  body <- paste("Minimum: ",minimum,"\tMaximum: ",maximum,"\nSequences filtered: ",number_filtered,
                "\nIncluding frame shifts: ",shifts,
                "\nLengths (occurrence): ",c_o,"\nIndividuals filtered: ",removed_num,"\t",n_o)
  write_to_log("FILTERED BY LENGTH",body,log)
  out
}

remove_some_seqs <- function(dt,pat = "",alleles = c(),log,pat_or_vec = "PAT"){
  if(pat_or_vec == "PAT"){
    alleles_names <- get_alleles_names(pat,alleles)
  } else {
    alleles_names <- alleles
  }
  out <- filter(dt,!(NAME %in% alleles_names))
  
  ###checking empty inds
  alleles_inds <- apply(out[,7:ncol(out)],2,function(x){sum(!is.na(x))})
  if(0 %in% alleles_inds){
    index <- which(alleles_inds == 0)
    inds_names <- names(alleles_inds)[index]
    removed_num <- length(index)
    out <- out[,-(index+6)]
  } else{
    inds_names <- c()
    removed_num <- 0
  }
  n_o <- paste(inds_names,collapse = " ")
  ###
  number_filtered <- nrow(dt)-nrow(out)
  ###
  body <- paste("Sequences filtered: ",number_filtered,
                "\n",paste(alleles_names,collapse = " "),
                "\nIndividuals filtered: ",removed_num,"\t",n_o)
  write_to_log("ADDITIONALLY REMOVED SEQS",body,log)
  out
}

filter_by_amplicon_depth <- function(dt,appx,log,minimum_depth = 2000){
  indexes <- which(as.numeric(appx[1,2:ncol(appx)]) < minimum_depth)+1
  del_names <- colnames(appx)[indexes]
  del_cov <- paste("(",appx[1,indexes],")",sep = "")
  del_all <- paste(del_names,del_cov)
  temp <- dt[,!(colnames(dt) %in% del_names)]
  kept_rows <- c()
  bad_counter <- 0
  for(i in 1:nrow(temp)){
    if(sum(!is.na(temp[i,7:ncol(temp)])) > 0){
      kept_rows <- c(kept_rows,i)
    } else{bad_counter <- bad_counter+1}
  }
  distr <- quantile(as.numeric(appx[1,2:ncol(appx)]),probs = seq(0,1,0.05))
  occur <- paste("(",distr,")",sep = "")
  c_o <- paste(paste(names(distr),occur),collapse = " ")
  body <- paste("Minimum depth: ",minimum_depth,"\nIndividuals filtered (coverage): ",length(del_names),
                "\t",paste(del_all,collapse = " "),"\nSequences filtered: ",bad_counter,"\nDistribution (value): ",
                c_o)
  write_to_log("FILTERED BY AMPLICON DEPTH",body,log)
  temp[kept_rows,]
}

filter_by_tag_switch <- function(dt,cp_dt,threshold,log_name){
  data <- dt[dt$MAX_FREQ >= (threshold/100),]
  seqs_removed <- nrow(dt)-nrow(data)
  f_dt <- cp_dt[cp_dt$MAX_FREQ >= (threshold/100),10:ncol(cp_dt)]
  sub <- data[,10:ncol(data)]
  sub[f_dt <= (threshold/100)] <- NA
  data[,10:ncol(data)] <- sub
  ###
  indexes <- apply(data,2,function(x){if(sum(is.na(x))==nrow(data)){0}else{1}})
  names_removed <- names(indexes)[indexes == 0]
  data <- data[,as.logical(indexes)]
  ###
  write_to_log("FILTERED BY PAF",body = paste("Chosen threshold:",threshold,"%\nNumber of sequences removed:",
                                              seqs_removed,"\nNumber of individuals removed:",length(names_removed),
                                              "\t",paste(names_removed,collapse = " ")),log = log_name)
  data
}

filter_by_tag_switch_dual <- function(dt,cp_dt,threshold,dual,log_name){
  data <- dt[dt$MAX_FREQ >= (threshold/100),]
  seqs_removed <- nrow(dt)-nrow(data)
  f_dt <- cp_dt[cp_dt$MAX_FREQ >= (threshold/100),10:ncol(cp_dt)]
  sub <- data[,10:ncol(data)]
  sub[f_dt <= (dual/100)] <- NA
  data[,10:ncol(data)] <- sub
  ###
  indexes <- apply(data,2,function(x){if(sum(is.na(x))==nrow(data)){0}else{1}})
  names_removed <- names(indexes)[indexes == 0]
  data <- data[,as.logical(indexes)]
  ###
  write_to_log("FILTERED BY PAF",body = paste("Chosen threshold:",threshold,"% Dual:",dual,"%\nNumber of sequences removed:",
                                              seqs_removed,"\nNumber of individuals removed:",length(names_removed),
                                              "\t",paste(names_removed,collapse = " ")),log = log_name)
  data
}

adding_freqs_info <- function(dt,appx){
  cp_dt <- dt
  out_list <- list()
  for(i in 1:nrow(dt)){
    freqs_allele <- c()
    for(j in 7:ncol(dt)){
      if(is.na(dt[i,j])==F){
        value <- as.numeric(dt[i,j])/as.numeric(appx[1,(colnames(dt)[j])])
        freqs_allele <- c(freqs_allele,value)
        cp_dt[i,j] <- value
      }
    }
    out_list[[length(out_list)+1]] <- c(min(freqs_allele),mean(freqs_allele),max(freqs_allele))
    #out_list[[length(out_list)+1]] <- round(c(min(freqs_allele),mean(freqs_allele),max(freqs_allele)),3)
  }
  combined <- as.data.frame(do.call(rbind,out_list))
  colnames(combined) <- c("MIN_FREQ","MEAN_FREQ","MAX_FREQ")
  list(cbind(combined,dt),cbind(combined,cp_dt))
}

adding_freqs_info_allele <- function(dt,appx){
  cp_dt <- dt
  out_list <- list()
  for(i in 1:nrow(dt)){
    freqs_allele <- c()
    for(j in 7:ncol(dt)){
      if(is.na(dt[i,j])==F){
        value <- as.numeric(dt[i,j])/as.numeric(appx[2,(colnames(dt)[j])])
        freqs_allele <- c(freqs_allele,value)
        cp_dt[i,j] <- value
      }
    }
    out_list[[length(out_list)+1]] <- round(c(min(freqs_allele),mean(freqs_allele),max(freqs_allele)),3)
  }
  combined <- as.data.frame(do.call(rbind,out_list))
  colnames(combined) <- c("MIN_FREQ","MEAN_FREQ","MAX_FREQ")
  list(cbind(combined,dt),cbind(combined,cp_dt))
}

plot_max_freq <- function(dt,file_name,w = NA,h = NA,t = 0.1){
  max_freqs1 <- sort(dt$MAX_FREQ,decreasing = T)
  max_freqs <- max_freqs1[max_freqs1<t]
  if(length(max_freqs) == 0){
    max_freqs <- max_freqs1
    print("BE AWARE!!! WHOLE RANGE!!!")
  }
  x <- c(1:length(max_freqs))
  p <- ggplot(mapping = aes(x,max_freqs))+theme_bw()+geom_point()+xlab("Order")+ylab("Maximum PAF")+
    ggrepel::geom_text_repel(aes(label = max_freqs),nudge_y = c(0.1,0.2,0.5))
  ggsave(file_name,plot = p,device = "png",width = w,height = h,units = "cm")
  print(paste("Plot saved in: ",getwd(),"/",file_name,sep = ""))
}

modify_names <- function(dt,gene){
  temp <- arrange(dt,desc(SAMPLES))
  abr <- stringr::str_trunc(analysed_file,3,side = "right",ellipsis = "")
  new_names <- paste(abr,gene,sprintf("%03d",c(1:nrow(dt))),sep = "_")
  temp$NAME <- new_names
  temp
}

modify_names_new <- function(dt,abr,gene,order){
  old_ff <- paste(paste(abr,gene,"alleles_names",order-1,sep = "_"),".csv",sep = "")
  out_name <- paste("../",paste(abr,gene,"alleles_names",order,sep = "_"),".csv",sep = "")
  temp <- arrange(dt,desc(SAMPLES))
  if(order == 1){
    new_names <- paste(abr,gene,sprintf("%03d",c(1:nrow(dt))),sep = "_")
    temp$NAME <- new_names
    out <- cbind(temp[,c("NAME","SEQUENCE")],"ORDER" = rep(order,times = nrow(temp)))
    write.csv(out,file = out_name,fileEncoding = "UTF-8",row.names = F)
  } else {
    old_f <- read.csv(paste("../",old_ff,sep = ""),stringsAsFactors = F)
    new_names <- c()
    for(i in 1:nrow(temp)){
      if(temp[i,"SEQUENCE"] %in% old_f$SEQUENCE){
        new_names[i] <- old_f$NAME[old_f$SEQUENCE == temp[i,"SEQUENCE"]]
      } else {
        new_names[i] <- paste(abr,gene,sprintf("%03d",nrow(old_f)+1),sep = "_")
        add_rowe <- data.frame("NAME" = new_names[i],"SEQUENCE" = temp[i,"SEQUENCE"],"ORDER" = order)
        old_f <- rbind(old_f,add_rowe)
      }
    }
    temp$NAME <- new_names
    write.csv(old_f,file = out_name,fileEncoding = "UTF-8",row.names = F)
  }
  temp
}

save_to_fasta <- function(dt,analyzed){
  headers <- dt$NAME
  seqs <- dt$SEQUENCE
  name <- paste(strsplit(analyzed,split = "\\.")[[1]][1],".fasta",sep = "")
  seqinr::write.fasta(as.list(seqs),headers,name,as.string = T)
  name
}

most_frequent <- function(vec){
  uniqs <- unique(vec)
  freqs <- c()
  for(i in uniqs){
    freqs <- c(freqs,sum(vec == i))
  }
  index <- which.max(freqs)
  uniqs[index]
}

detect_stop_codons <- function(seqs,log,frame_val,stop_verify = F){
  translated <- lapply(seqs,seqinr::translate,frame = frame_val)
  stops <- c()
  probable_stops <- c()
  for(i in 1:length(translated)){
    if("*" %in% translated[[i]]){
      index <- which.max(translated[[i]] == "*")
      index_xx <- which(translated[[i]] == "X")
      #if "X" is present after
      if(length(index_xx) == 0 || sum(index_xx < index) == 0){
        stops <- c(stops,names(translated)[i])
      } else { #if "X" is present before
        probable_stops <- c(probable_stops,names(translated)[i])
        if(stop_verify == T){
          print(paste("Possible STOP in:",names(translated)[i]))
        } else {
          stops <- c(stops,names(translated)[i])
        }
      }
    }
  }
  probable_stops <- paste(probable_stops,collapse = " ")
  write_to_log("SEQUENCES WITH STOP CODONS",c(paste(stops,collapse = " "),paste("Frame:",frame_val),
                                              paste("Probable:",probable_stops)),log)
  print("With_stop / all")
  print(paste(length(stops),"/",length(seqs)))
  stops
}

detect_frame_shifts <- function(dt,des_len){
  f_s <- c()
  for(i in 1:nrow(dt)){
    #print((dt$LENGTH[i] - des_len))
    if((dt$LENGTH[i] - des_len) %% 3 != 0){
      f_s <- c(f_s,dt$NAME[i])
    }
  }
  write_to_log(title = "WITH FRAME SHIFTS",body = paste(c(f_s),collapse = " "),log = log_name)
  f_s
}

write_to_log <- function(title,body,log){
  out_log <- paste(title,"\n",paste(body,collapse = "\n"),"\n\n",paste(rep("*",times = 100),collapse = ""),"\n",sep = "")
  write(out_log,log,append = T)
}

allign <- function(seqs,out_file){
  seqs_collapsed <- sapply(seqs,paste,collapse = "")
  allignment <- msa::msaClustalW(seqs_collapsed,type = "dna")
  converted <- msa::msaConvert(allignment,type = 'seqinr::alignment')
  outcome <- ape::as.DNAbin.alignment(converted)
  ape::write.FASTA(outcome,out_file)
  outcome
}

align_mafft <- function(in_file,out_file,rev_comp = T,path_mafft = "/home/tomek/miniconda3/bin/mafft",method_v = "localpair"){
  initial_bin <- ape::read.FASTA(in_file)
  #alligned <- ips::mafft(initial_bin,exec = "/usr/bin/mafft",op = 5,ep = 3,options = c("--adjustdirection"))
  if(rev_comp == T){
  aligned <- ips::mafft(initial_bin,exec = path_mafft,op = 5,ep = 3,options = c("--adjustdirection"),method = method_v)
  } else {
    aligned <- ips::mafft(initial_bin,exec = path_mafft,op = 5,ep = 3,method = method_v)
  }
  ape::write.FASTA(aligned,out_file)
  aligned
}

build_tree <- function(alligned_DNAbin,tree_name,multi = T,boot = T){
  require(phangorn)
  phyD <- phyDat(alligned_DNAbin)
  dm <- dist.ml(phyD,model = "JC69",exclude = "pairwise")
  nj <- NJ(dm)
  fit <- pml(nj,phyD)
  if(boot == T){
    fitJC <- optim.pml(fit, model = "JC", rearrangement = "ratchet")
    bs <- bootstrap.pml(fitJC, bs=500, optNni=TRUE, multicore=multi,mc.cores = 10,control = pml.control(trace = 0))
    consensus_tree <- plotBS(midpoint(fitJC$tree), bs, type="none",p = -1)
    ape::write.tree(consensus_tree,tree_name)
    consensus_tree
  } else {
    fitJC <- optim.pml(fit, model = "JC", rearrangement = "none")
    ape::write.tree(fitJC$tree,tree_name)
    fitJC$tree
  }
}

build_tree2 <- function(alligned_DNAbin,tree_name,root_v){
  get_nj_tree <- function(x){root(phy = nj(dist.dna(x,model = "JC69",pairwise.deletion = T)), outgroup = root_v,resolve.root = TRUE)}
  nj <- get_nj_tree(alligned_DNAbin)
  nj$node.label <- boot.phylo(phy = nj,x = alligned_DNAbin,FUN = get_nj_tree,mc.cores = 5)
  ape::write.tree(nj,tree_name)
  nj
}

create_gg_tree <- function(treee,non_functional = c(),frame_shi = c(),removed = c(),
                           image_file,x_lim = 0.5,height_im = NA){
  tip_labels <- treee$tip.label
  to_col <- rep("functional",times = length(treee$tip.label))
  if(length(frame_shi) > 0){
    indexes_r <- which(treee$tip.label %in% frame_shi)
    to_col[indexes_r] <- "shifted"
  }
  if(length(removed) > 0){
    indexes_r <- which(treee$tip.label %in% removed)
    to_col[indexes_r] <- "removed"
  }
  if(length(non_functional) > 0){
    indexes <- which(treee$tip.label %in% non_functional)
    to_col[indexes] <- "non-functional"
  }
  meta <- as.data.frame(cbind(tip_labels,to_col))
  
  plot <- ggtree(treee) %<+% meta + geom_tiplab(hjust = -0.2,size = 1.5)+ geom_nodelab(geom = "label",hjust = 0.9,size = 1.5) + 
    geom_tippoint(aes(x =x + 0.01*max(x),color = to_col),size = 1.5) + 
    scale_color_manual(values = c("functional" = "green","non-functional" = "red","removed" = "orange","shifted" = "purple")) + 
    theme(legend.title = element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size = 5))+
    xlim(0,x_lim)+geom_treescale(y = -1,linesize = 1,fontsize = 1.5,width = 0.05)
  ggsave(image_file,plot,device = "png",scale = 1,limitsize = F,height = height_im)
  print(paste("Plot saved in: ",getwd(),"/",image_file,sep = ""))
}

create_gg_tree_circu <- function(treee,list_samp = list(),
                                 image_file,x_lim = 0.5,height_im = NA,width_im = NA){
  env_v <- environment()
  '%src%'("treee$node.label","<90","'NA'",env_v)
  colors_v <- c("black","green","red")
  
  plot <- ggtree(treee,layout = "circular",size = 2) +
    geom_nodelab(aes(label = label,subset = (!isTip & label != "NA")),geom = "label") + 
    geom_tiplab2(aes(label = label,subset = label %in% "Novi"),hjust = -0.8,size = 15)+
    theme(legend.position = "none",plot.margin = unit(c(0,0,0,0),"cm"))+
    #negative xlim increases the radius of circular tree
    xlim(-0.1,x_lim)+geom_treescale(y = 0,linesize = 1,fontsize = 1.5,width = 0.05)
  if(length(list_samp)>0){
    for(i in 1:length(list_samp)){
      #subsetting relies on data within 'plot' - parent, node, branch.length, label
      plot <- plot + geom_point2(aes(subset = (isTip & label %in% list_samp[[i]]),x =x + 0.01*max(x)),
                                 col = colors_v[i],size = 6)
      # this part is needed to not let the points layer override previous 
      qs <- paste('(isTip & label %in% list_samp[[',i,']])',sep = "")
      plot$layers[[length(plot$layers)]]$mapping$subset <- rlang::expr_interp(rlang::parse_quo(qs,env = env_v))
    }
  }
  ggsave(image_file,plot,device = "png",scale = 1,limitsize = F,height = height_im,width = width_im)
  print(paste("Plot saved in: ",getwd(),"/",image_file,sep = ""))
}

change_to_genotypes <- function(dt,dnabin_org,genotypes_file){
  data <- dt[,9:ncol(dt)]
  out_list <- list()
  out_list[["individuals"]] <- colnames(data)[2:ncol(data)]
  dnabin <- as.list(dnabin_org)
  for(i in names(dnabin)){
    out_vec <- c()
    row_index <- which(data$NAME == i)
    for(j in 2:ncol(data)){
      if(is.na(data[row_index,j]) == F){
        out_vec <- c(out_vec,1) 
      } else {
        out_vec <- c(out_vec,0)
      }
    }
    out_list[[i]] <- out_vec
  }
  genotypes <- as.data.frame(do.call(cbind,out_list))
  write.csv(genotypes,genotypes_file,quote = F,fileEncoding = "UTF-8",row.names = F)
  genotypes
}

get_alleles_names <- function(pattern,numbers){
  paste(pattern,sprintf("%03d",numbers),sep = "_")
}

protein_alignment <- function(unalligned_DNA,not_needed,framing,out_file){
  translated <- lapply(unalligned_DNA,translate_wGaps,frame = framing)
  translated_removed <- subset(translated,!(names(translated) %in%  not_needed))
  seqs_collapsed <- sapply(translated_removed,paste,collapse = "")
  allignment <- msa::msaClustalW(seqs_collapsed,type = "protein")
  converted <- msa::msaConvert(allignment,type = 'seqinr::alignment')
  seqinr::write.fasta(as.list(converted$seq),converted$nam,out_file)
  print(paste("Alignment saved in: ",getwd(),"/",out_file,sep = ""))
}

freqs <- function(vec){
  labels <- sort(unique(vec))
  freqs_out <- sapply(labels,function(x){sum(vec == x)})
  names(freqs_out) <- labels
  freqs_out
}

remove_inds <- function(dt,to_remove,log,mess = "INDIVIDUALS REMOVED BY USER"){
  indexes <- c(1:ncol(dt))[colnames(dt) %in% to_remove]
  out <- dt[,-indexes]
  ### checking sequences
  seqs <- apply(out[,7:ncol(out)],1,function(x){sum(!is.na(x))})
  if(0 %in% seqs){
    index <- which(seqs == 0)
    removed_num <- length(index)
    out <- out[-index,]
  } else{removed_num <- 0}
  ###
  n_o <- paste(colnames(dt)[indexes],collapse = " ")
  body <- paste("Individuals removed: ",length(indexes),"\t",n_o,"\nSequences removed: ",removed_num,sep = "")
  write_to_log(mess,body,log)
  out
}

########## FOR NICE TABLE #####
read_it <- function(dt){
  my_data <- as.data.frame(t(dt[,-c(1:6)]))
  my_data <- cbind(as.integer(colnames(dt)[-c(1:6)]),my_data)
  colnames(my_data) <- c("ID",dt$NAME)
  rownames(my_data) <- NULL
  my_data
}

add_species <- function(dt,file){
  species_info <- read.table(file,sep = ",",header = T)[,c("individual_id","species")]
  joined <- left_join(dt,species_info,by = c("ID" = "individual_id")) %>% select(ID,species,everything())
  joined
}

count_filled_rows <- function(dt){
  raw <- dt[,3:ncol(dt)]
  alleles_number <- rowSums(!is.na(raw))
  cbind(alleles_number,dt)
}

freqs <- function(vec){
  labels <- sort(unique(vec))
  freqs_out <- sapply(labels,function(x){sum(vec == x)})
  names(freqs_out) <- labels
  freqs_out
}

detect_relatedness_allele <- function(x,species){
  occurs <- as.character(species[!is.na(x)])
  if(length(unique(occurs)) > 1){
    paste(unique(occurs),collapse = "_")
  } else {
    unique(occurs)
  }
}

detect_relatedness_ind <- function(x,all_info){
  occurs <- all_info[!(is.na(x))]
  freqs_out <- freqs(occurs)
  index <- which.max(freqs_out)
  c(round(freqs_out[index]/sum(freqs_out),2),names(freqs_out)[index])
}

detect_others_ind <- function(x,all_info,species_var){
  occurs <- all_info[!(is.na(x))]
  splited <- strsplit(occurs,"_")
  lapply(splited, function(x){logic_initial <- species_var %in% x
  if(length(x)>2){logic_initial <- c(logic_initial,T)
  } else {logic_initial <- c(logic_initial,F)}
  if(length(x)==length(species_var)){return(c(logic_initial,T))
  } else {return(c(logic_initial,F))}
  }) %>% as.data.frame() %>% apply(1,sum)
}

adjust_alleles_col <- function(dt,range,desired_colname,all_allles){
  dt_cp <- dt
  for(i in range[1]:range[2]){
    actual_name <- colnames(dt)[i]
    is_species <- actual_name == desired_colname
    if(is_species == T){
      remember <- i
    } else {
      dt_cp[,i] <- dt[,i] - all_allles
    }
  }
  dt_cp[,-remember]
}

adjust_col <- function(reac_object,dt,col_number,funkcja_x){
  ro_cp <- reac_object
  style <- list()
  for(i in 1:nrow(dt)){
    decision <- funkcja_x(dt[i,col_number],dt[i,]) #vector change_name oraz change
    if(is.null(decision)){
      style[i] <- NULL
    } else{
      added_list <- list(decision[2])
      names(added_list) <- decision[1]
      style[[i]] <- added_list
    }
  }
  ro_cp[["x"]][["tag"]][["attribs"]][["columns"]][[col_number]][["style"]] <- style
  ro_cp
}

contamination_checking <- function(dt,id_file,dup_symbols = " "){
  dt <- dt[,!grepl(paste(dup_symbols,collapse = "|"),colnames(dt))]
  dt <- dt[rowSums(dt[,7:ncol(dt)],na.rm = T) > 0,]
  my_data_old <- read_it(dt) %>% add_species(id_file) %>% count_filled_rows()
  rel_allele <- apply(my_data_old[,-(1:3)],2,detect_relatedness_allele,my_data_old$species)
  
  rel_ind <- t(apply(my_data_old[,-(1:3)],1,detect_relatedness_ind,rel_allele)) %>% as.data.frame()
  colnames(rel_ind) <- c("degree","related_to")
  rel_ind$degree <- as.numeric(as.character(rel_ind$degree))
  rel_ind$related_to <- as.character(rel_ind$related_to)
  my_data_old <- cbind(rel_ind,my_data_old)
  
  species_names <- as.character(sort(unique(my_data_old$species)))
  rel_oth <- t(apply(my_data_old[,-(1:5)],1,detect_others_ind,rel_allele,species_names)) %>% as.data.frame()
  colnames(rel_oth) <- c(paste("alleles",species_names,sep = "_"),"more_than_2","alleles_all_species")
  my_data_old <- cbind(my_data_old[,1:5],rel_oth,my_data_old[,-(1:5)])
  
  View(arrange(my_data_old,desc(alleles_number)))
  t_all_num <- as.integer(readline(prompt = "Threshold for the number of alleles: "))
  chosen_colors <- rainbow(length(unique(my_data_old$species)))
  
  
  for(i in species_names){
    my_data <- filter(my_data_old,species == i)
    des_cn <- paste("alleles",i,sep = "_")
    my_data <- adjust_alleles_col(my_data,c(6,5+length(species_names)),des_cn,my_data$alleles_all_species)
    
    nice_table <- reactable::reactable(my_data,bordered = T,defaultPageSize = nrow(my_data),outlined = T) %>%
      adjust_col(my_data,col_number = 3,function(value,row){
        if(value >= t_all_num){
          c("background","red")
        } else {NULL}
      }) %>%
      adjust_col(my_data,col_number = 2,function(value,row){
        if(value != row[,5]){
          c("background","red")
        } else {NULL}
      }) %>%
      adjust_col(my_data,col_number = 1,function(value,row){
        if(value >= .8){
          c("background","green")
        } else if(value >= 0.7){
          c("background","yellow")
        } else if(value >= 0.6){
          c("background","orange")
        } else {
          c("background","red")
        }
      }) %>%
      adjust_col(my_data,col_number = 5,function(value,row){
        if(value %in% species_names){
          index <- which(species_names == value)
          return(c("background",chosen_colors[index]))
        } else {NULL}
      })
    
    starting <- 7+length(species_names)
    for(j in starting:ncol(my_data)){
      nice_table <- adjust_col(nice_table,my_data,col_number = j,function(value,row){
        if(rel_allele[j-(starting-1)] %in% species_names){
          index <- which(species_names == rel_allele[j-(starting-1)])
          app_col <- chosen_colors[index]
          return(c("background",app_col))
        } else {
          present_species <- strsplit(rel_allele[j-(starting-1)],split = "_")[[1]]
          if(length(present_species) == 2 && is.na(value) == F){
            index <- which(species_names == present_species[present_species != i])
            app_col <- chosen_colors[index]
            return(c("background",app_col))
          } else {NULL}
        }
      })
    }
    chosen_name <- paste("contamination_checking_",i,".html",sep = "")
    chosen_name <- stringr::str_replace(chosen_name,"/","_")
    htmlwidgets::saveWidget(nice_table,chosen_name)
  }
}

######### FOR DUPLICATES #####
set_colnames <- function(dt,col_nm){
  colnames(dt) <- col_nm
  dt
}

read_and_change <- function(dt_vs,appx_vs,dup_symbols){
  outcome <- list()
  for(i in 1:length(dt_vs)){
    dt_v <- read.csv(dt_vs[i],stringsAsFactors = F,check.names = F)[,-c(1:3,5:8)] %>%
      #replace(is.na(.), 0) %>%
      pivot_longer(cols = -c(SEQUENCE,NAME),names_to = "ID",values_to = "n_reads") %>%
      mutate(order = i) %>% na.omit()
    
    cov_table <- group_by(dt_v,ID) %>% summarise(cov_alleles = sum(n_reads))
    appx_v <- t(read.csv(appx_vs[i],stringsAsFactors = F,check.names = F,header = F))[-1,] %>% 
      as.data.frame(row.names = NA) %>% set_colnames(c("ID","cov_total","cov_alleles","alleles_num")) %>%
      mutate_all(as.character) %>% mutate_at(.vars = c(2:4),as.integer)
    cov_table <- left_join(cov_table,appx_v[,1:2],by = "ID")
    
    dt_v <- left_join(dt_v,cov_table,by = "ID") %>% mutate(PAF = n_reads/cov_total)
    
    for(j in dup_symbols){
      indexes <- grep(j,dt_v$ID)
      dt_v$ID[indexes] <- gsub(j,"",dt_v$ID[indexes])
      dt_v$order[indexes] <- paste(i,j,sep = "")
    }
    outcome[[length(outcome)+1]] <- dt_v
  }
  as.data.frame(do.call(rbind,outcome))
}

take_out_resequenced <- function(dt){
  kept <- c()
  for(i in unique(dt$ID)){
    subset <- dt[dt$ID == i,]
    if(length(unique(subset$order))>1){
      kept <- c(kept,i)
    }
  }
  filter(dt,ID %in% kept) %>% group_by(ID,NAME) %>% mutate(num_rep = n()) %>% ungroup()
}

reseq_info <- function(dt){
  out_list <- list()
  reps <- sort(unique(dt$order))
  for(i in unique(dt$NAME)){
    out_list[[length(out_list)+1]] <- c(i,as.integer(reps %in% dt$order[dt$NAME == i]))
  } 
  outcome <- as.data.frame(do.call(rbind,out_list))
  colnames(outcome) <- c("NAME",reps)
  outcome
}

discrepancy <- function(dt,file_nm){
  out_list <- list()
  combinations <- combn(unique(dt$order),2)
  for(j in 1:ncol(combinations)){
    dt_ft <- filter(dt, order %in% combinations[,j]) %>% take_out_resequenced()
    for(i in unique(dt_ft$ID)){
      temp <- filter(dt_ft,ID == i)
      all_alleles <- length(unique(temp$NAME))
      n_comb1 <- length(unique(temp$NAME[temp$order == combinations[1,j]]))
      n_comb2 <- length(unique(temp$NAME[temp$order == combinations[2,j]]))
      shared <- sum(temp$NAME[temp$order == combinations[1,j]] %in% temp$NAME[temp$order == combinations[2,j]])
      DR <- round(1-(shared/all_alleles),2)
      out_list[[length(out_list)+1]] <- c(i,combinations[,j],all_alleles,n_comb1,n_comb2,shared,DR)
    }
  }
  outcome <- as.data.frame(do.call(rbind,out_list)) %>% set_colnames(c("ID","comb1","comb2","n_alleles","n_comb1","n_comb2","n_shared","DR"))
  write.csv(outcome,file_nm,row.names = F,fileEncoding = "UTF-8")
}

duplicate_gg_tree <- function(treee,dt,info,image_file,x_lim = 0.5,boot = T,height_im = NA){
  reps <- sort(unique(dt$order))
  PAF_list <- list()
  tip_labels <- treee$tip.label
  to_col <- rep("unrepeated",times = length(treee$tip.label))
  for(i in 1:length(tip_labels)){
    PAF_list[[length(PAF_list)+1]] <- rep(0,times = length(reps))
    temp <- filter(dt,NAME == tip_labels[i]) %>% group_by(order) %>% summarise(PAFer = mean(PAF)) %>%
      arrange(order)
    PAF_list[[length(PAF_list)]][reps %in% temp$order] <- temp$PAFer
    if(tip_labels[i] %in% info$NAME){
      indexes <- which(info[info$NAME == tip_labels[i],] == 1)
      if(length(indexes) > 1){
        to_col[i] <- paste(colnames(info)[indexes],collapse = "_")
      }
    }
  }
  PAFs <- as.data.frame(do.call(rbind,PAF_list)) %>% set_colnames(reps) %>% round(3) %>% 
    apply(1,paste,collapse = " ")
  meta <- as.data.frame(cbind(tip_labels,to_col,PAFs))
  meta$PAFs <- as.character(meta$PAFs)
  
  if(boot == T){
    plot <- ggtree(treee) %<+% meta + geom_tiplab(hjust = -0.2,size = 1.5)+ 
      geom_nodelab(geom = "label",hjust = 0.9,size = 1.5) + 
      geom_tippoint(aes(x =x + 0.01*max(x),color = to_col),size = 1.5) + 
      scale_color_brewer(palette = "Set1") + 
      theme(legend.title = element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size = 5))+
      xlim(0,x_lim)+geom_treescale(linesize = 1,fontsize = 1.5,width = 0.05)+
      geom_tiplab(aes(label = c(meta$PAFs,rep(NA,times = length(meta$PAFs)-1))),size = 1.5,offset = x_lim/6)
  } else {
    plot <- ggtree(treee) %<+% meta + geom_tiplab(hjust = -0.2,size = 1.5)+ 
      geom_nodelab(geom = "label",hjust = 0.9,size = 1.5) + 
      geom_tippoint(aes(x =x + 0.01*max(x),color = to_col),size = 1.5) + 
      scale_color_brewer(palette = "Set1") + 
      theme(legend.title = element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size = 5))+
      xlim(0,x_lim)+geom_treescale(linesize = 1,fontsize = 1.5,width = 0.05)+
      geom_tiplab(aes(label = c(meta$PAFs,rep(NA,times = length(meta$PAFs)-2))),size = 1.5,offset = x_lim/6)
  }
  
  
  ggsave(image_file,plot,device = "png",scale = 1,limitsize = F,height = height_im)
  print(paste("Plot saved in: ",getwd(),"/",image_file,sep = ""))
}

check_reversal <- function(dt,analysed_file_v){
  p <- save_to_fasta(dt,analysed_file_v)
  ala <- align_mafft(p,"temp.fasta")
  pp <- seqinr::read.fasta("temp.fasta",forceDNAtolower = F)
  ppp <- lapply(pp, function(i){paste(i[i %in% c("A","T","C","G")],collapse = "")})
  if(length(names(ppp)[grepl("_R_",names(ppp))]) > 0){
    print("DETECTED AND REPLACED!!!")
    changed <- sapply(names(ppp)[grepl("_R_",names(ppp))],substring,first = 4)
    for(j in 1:length(changed)){
      dt$SEQUENCE[dt$NAME == changed[j]] <- unlist(ppp[names(ppp) == names(changed)[j]])
    }
    if(length(unique(dt$SEQUENCE)) != nrow(dt)){
      print("AMPLISAS SUCKS!!! YOU HAVE REPEATED SEQUENCES!!!")
    }
  }
  file.remove("temp.fasta")
  dt
}

check_reversal_alt <- function(dt,analysed_file_v,path_mafft = "/home/tomek/miniconda3/bin/mafft"){
  #concordance with previous run
  alt_f <- sort(list.files("../")[grepl(".csv",list.files("../"))])
  alt_f <- paste("../",alt_f[length(alt_f)],sep = "")
  alt_f <- read.csv(alt_f,stringsAsFactors = F)
  alt_f$NAME <- paste("PRE",alt_f$NAME,sep = "_")
  save_to_fasta(alt_f,"pre")
  
  p <- save_to_fasta(dt,analysed_file_v)
  system(paste("cat ",p," >> pre.fasta",sep = ""),wait = T)
  
  ala <- align_mafft("pre.fasta","temp.fasta",path_mafft = path_mafft)
  pp <- seqinr::read.fasta("temp.fasta",forceDNAtolower = F)
  ppp <- lapply(pp, function(i){paste(i[i %in% c("A","T","C","G")],collapse = "")})
  
  #making sure pre is not reversed
  if(length(names(ppp)[grepl("_R_PRE",names(ppp))]) > 0){
    ppp <- lapply(ppp,function(x){system(paste("echo ",x," | tr ACGTMRWSYKVHDBN TGCAKYWSRMBDHVN | rev"),
                                         wait = T, intern = T)}) 
  }
  ppp <- ppp[!grepl("PRE",names(ppp))]
  
  if(length(names(ppp)[grepl("_R_",names(ppp))]) > 0){
    print("DETECTED AND REPLACED!!!")
    changed <- sapply(names(ppp)[grepl("_R_",names(ppp))],substring,first = 4)
    for(j in 1:length(changed)){
      dt$SEQUENCE[dt$NAME == changed[j]] <- unlist(ppp[names(ppp) == names(changed)[j]])
    }
    if(length(unique(dt$SEQUENCE)) != nrow(dt)){
      print("AMPLISAS SUCKS!!! YOU HAVE REPEATED SEQUENCES!!!")
    }
  }
  #file.remove("temp.fasta")
  #file.remove("pre.fasta")
  dt
}

Are_Ns <- function(dt){
  indexes <- grepl("N",dt$SEQUENCE)
  if(sum(indexes) > 0){
    print("Ns detected!!!")
    print(paste(dt$NAME[indexes],collapse = " "))
  }
}

limit_to_species <- function(dt,sp_id_file,desired_sp,log_name,dup_symbol = "d"){
  sp <- read.csv(sp_id_file) %>% filter(species %in% desired_sp)
  imiona <- colnames(dt)[7:ncol(dt)]
  rem <- imiona[!(imiona %in% c(sp$individual_id,paste(sp$individual_id,"d",sep = "")))]
  remove_inds(dt,rem,log_name,mess = paste("INDIVIDUALS REMOVED BASED ON SPECIES\nAccepted species:",
                                           paste(desired_sp,collapse = " ")))
}
