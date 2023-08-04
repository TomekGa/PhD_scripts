#Provides consensus sequence for provided sequences

seqs <- c("GAGGGAGAGAGTCTGTGT","GARGGAGAGATACTGTCT") #ADJUST

# ALTERNATIVE - directly from amplisas
seqs <- "AAWABRCGCCTGGAGTGC
AAATAGCGCAGCGAGTGC"
seqs <- strsplit(seqs,"\n")[[1]]


split_degeneration <- function(nucs,dt){
  nucs_cp <- nucs
  if(sum(nucs %in% dt$code) > 0){
    indexes <- which(nucs %in% dt$code)
    for(i in indexes){
      to_change <- nucs_cp[i]
      to_replace <- strsplit(dt$description[dt$code == to_change],split = "")[[1]]
      nucs_cp <- c(nucs_cp,to_replace)
    }
    return(nucs_cp[-indexes])
  } else {
    return(nucs_cp)
  }
}

get_degeneration <- function(nucs,dt){
  simpler <- unique(nucs)
  dt_list <- sapply(dt$description,strsplit,split = "")
  names(dt_list) <- dt$code
  if(length(simpler)>1){
    for(i in 1:length(dt_list)){
      if(sum(simpler %in% dt_list[[i]]) == length(simpler) & length(simpler) == length(dt_list[[i]])){
        return(names(dt_list)[i])
      }
    }
  } else {
    return(simpler)
  }
}

get_primer <- function(seqs){
  dt_degenerate <- data.frame("code" = c("M","R","W","S","Y","K","V","H","D","B","N"),
                              "description" = c("AC","AG","AT","CG","CT","GT","ACG","ACT","AGT","CGT","ACGT"),
                              stringsAsFactors = F)
  p <- sapply(seqs,strsplit,split = "")
  out <- c()
  for(i in 1:length(p[[1]])){
    nucs <- c()
    for(j in 1:length(p)){
      nucs <- c(nucs,p[[j]][i])
    }
    nucs <- split_degeneration(nucs,dt_degenerate)
    out <- c(out,get_degeneration(nucs,dt_degenerate))
  }
  print(paste(out,collapse = ""))
}

get_primer(seqs)
