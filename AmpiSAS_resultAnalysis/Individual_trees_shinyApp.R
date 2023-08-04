#Shiny app to interactively display the tree of all alleles and mark the alleles of a single individual

#REQUIRES
#the analysed AmpliSAS outputs with Amplisas_output_analysis.R script

library(shiny)
library(ggplot2)
library(ape)
library(dplyr)
library(ggtree)
library(shinyFiles)

path_var <- 'exemplary_data' #ADJUST

#FUNCTIONS ####

create_gg_tree <- function(treee,present = c(),PAFs = c(),x_lim = 0.5,height_im = NA,PAF_thr = 10){
    to_col <- rep("brak",times = length(treee$tip.label))
    
    if(length(present) > 0){
        for(i in 1:length(present)){
            index <- which(treee$tip.label == present[i])
            treee$tip.label[index] <- paste(treee$tip.label[index],round(PAFs[i],2),sep = " ")
            if(PAFs[i] >= PAF_thr){
                to_col[index] <- "above"
            } else {
                to_col[index] <- "below" 
            }
        }
    }
    tip_labels <- treee$tip.label
    meta <- as.data.frame(cbind(tip_labels,to_col))
    
    plot_out <- ggtree(treee) %<+% meta + geom_tiplab(hjust = -0.2,size = 3)+ geom_nodelab(geom = "label",hjust = 0.9,size = 1.5) +
        geom_tippoint(aes(x =x + 0.01*max(x),color = to_col),size = 3) +
        scale_color_manual(values = c("brak" = NA,"above" = "green","below" = "red")) +
        theme(legend.title = element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size = 10))+
        xlim(0,x_lim)+geom_treescale(y = -1,linesize = 1,fontsize = 1.5,width = 0.05)+
        labs(title = paste("PAF ",PAF_thr,"%",sep = ""))
    plot_out
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

# Define UI  ####
ui <- fluidPage(
    
    # Application title,
    sidebarLayout(
        sidebarPanel(
            titlePanel("Individual trees"),
            headerPanel(shinyDirButton("dir", "Chose directory", "Upload")),
            wellPanel(textInput(inputId = "ind","Individual"),
                      actionButton("random","Random individual",class = "btn-primary"),
                      actionButton("prev","Previous individual",class = "btn-primary"),
                      actionButton("nexta","Next individual",class = "btn-primary")),
            radioButtons("PAF","Type of PAF",choices = c("Sequence" = "seq","Individual" = "indi"),
                         selected = "seq"),
            sliderInput(inputId = "xlim","x_lim",0.2,2.5,value = 0.5,0.1),
            numericInput(inputId = "PAF_threshold","PAF threshold [%]",value = 2),
            sliderInput(inputId = "hei","Plot height [px]",500,2000,value = 1000,100),
            textOutput("druk"),textOutput("num"),textOutput("pres"),textOutput("pafy")),
        mainPanel(plotOutput("tree"))
    )
)

# Define server logic ####
server <- function(input, output,session) {
    
    #shinyDirChoose(input, 'dir', roots = c(Shared = '/home/tomek/Dropbox/MHC_projekt/Shared_MHC_project'))
    shinyDirChoose(input, 'dir', roots = c(root = path_var))
    dir <- reactive(paste(path_var,paste(unlist(input$dir)[-c(1,length(unlist(input$dir)))],collapse = "/"),sep = "/"))
    #output$dir <- renderPrint(dir())
    
    output$druk <- renderPrint({
        print(dir())
    })
    
    hei <- reactive(input$hei)
    
    observeEvent(input$random,{
        b <- list.files(path = dir(),pattern = "modified.csv")
        whole_p <- paste(dir(),b,sep = '/')
        mod <- read.csv(whole_p,stringsAsFactors = F,check.names = F)
        nms <- colnames(mod)[-c(1:9)]
        
        updateTextInput(session,"ind",value = as.numeric(sample(nms,1,replace = F)))
    })
    
    observeEvent(input$nexta,{
        b <- list.files(path = dir(),pattern = "modified.csv")
        whole_p <- paste(dir(),b,sep = '/')
        mod <- read.csv(whole_p,stringsAsFactors = F,check.names = F)
        nms <- colnames(mod)[-c(1:9)]
        
        if((as.character(input$ind) %in% nms) == F){
            indx <- 0
        } else{
            indx <- which(nms == as.character(input$ind)) 
        }
        
        if(indx == length(nms)){
            indx <- 1
        }
        
        updateTextInput(session,"ind",value = nms[indx+1])
    })
    
    observeEvent(input$prev,{
        b <- list.files(path = dir(),pattern = "modified.csv")
        whole_p <- paste(dir(),b,sep = '/')
        mod <- read.csv(whole_p,stringsAsFactors = F,check.names = F)
        nms <- colnames(mod)[-c(1:9)]

        if((as.character(input$ind) %in% nms) == F){
            indx <- 2
        } else{
            indx <- which(nms == as.character(input$ind))
        }

        if(indx == 1){
            indx <- length(nms)
        }

        updateTextInput(session,"ind",value = nms[indx-1])
    })

    output$tree <- renderPlot({
        b <- list.files(path = dir(),pattern = "tree.new")
        whole_p <- paste(dir(),b,sep = '/')
        plik <- ape::read.tree(whole_p)

        b <- list.files(path = dir(),pattern = "modified.csv")
        whole_p <- paste(dir(),b,sep = '/')
        mod <- read.csv(whole_p,stringsAsFactors = F,check.names = F)

        b <- list.files(path = dir(),pattern = "appendix.csv")
        whole_p <- paste(dir(),b,sep = '/')
        appxa <- read.csv(whole_p,stringsAsFactors = F,check.names = F)

        colnm <- as.character(input$ind)
        freqs <- adding_freqs_info(dt = mod[,-c(1:3)],appx = appxa)[[2]]
        pres <- freqs$NAME[!is.na(freqs[,colnm])]

        output$num <- renderText({
            print(paste("Number of alleles:\n",length(pres)))
        })

        output$pres <- renderText({
            print(paste("Present alleles:\n",paste(pres,collapse = " ")))
        })

        if(as.character(input$PAF) == "seq"){
            pafy <- (freqs$MAX_FREQ[!is.na(freqs[,colnm])])*100
        } else {
            pafy <- freqs[!is.na(freqs[,colnm]),colnm]*100
        }

        output$pafy <- renderText({
            print(paste("PAFs:\n",paste(round(pafy,2),collapse = " ")))
        })

        create_gg_tree(plik,PAF_thr = input$PAF_threshold,present = pres,PAFs = pafy,x_lim = input$xlim,height_im = input$hei)
    },width = 600,height = function(){input$hei}
)
    #height = function(){reactive(input$hei)}
}

# Run the application 
shinyApp(ui = ui, server = server)
