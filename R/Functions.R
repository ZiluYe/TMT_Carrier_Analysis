
library(data.table)
library(magrittr)
library(tidyverse)
library(seqinr)
library(ggpubr)
library(plyr)
library(corrplot)
library(ggprism)

theme_set(theme_bw() +theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()))

Bind <- function(Files){
  dataset <- data.table()
  
  for (i in 1:length(Files)){
    temp_data <- fread(Files[i], stringsAsFactors = F)
    temp_data[,Sample := Files[i]]
    dataset <- rbindlist(list(dataset, temp_data), use.names = FALSE) 
  }
  
  return(dataset)
}

PepNumber <- function (sequence, enzyme = "trypsin", missed = 0) {
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  if (enzyme == "trypsin") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_string <- gsub("KP", "!P", seq_string)
    seq_string <- gsub("RP", "!P", seq_string)
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "trypsin.strict") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (length(stop) == 0) 
    warning("sequence does not contain cleavage sites")
  cleave <- function(sequence, start, stop, misses) {
    peptide <- substring(sequence, start, stop)
    mc <- rep(misses, times = length(peptide))
    result <- data.frame(peptide, start, stop, mc, stringsAsFactors = FALSE)
    return(result)
  }
  start <- c(1, start)
  stop <- c(stop, end_position)
  results <- cleave(sequence, start, stop, 0)
  if (missed > 0) {
    for (i in 1:missed) {
      start_tmp <- start[1:(length(start) - i)]
      stop_tmp <- stop[(1 + i):length(stop)]
      peptide <- cleave(sequence, start_tmp, stop_tmp, 
                        i)
      results <- rbind(results, peptide)
    }
  }
  results$length <- str_length(results$peptide)
  
  results <- results[results$length >= 6 & results$length <= 30,]
  
  return(dim(results)[1])
}



FASTAtoDF <- function(Fasta){
  DB <- read.fasta(Fasta, as.string = TRUE, seqtype = "AA")
  DB <- as.data.frame(do.call(rbind, DB))
  DB$Accessions <- row.names.data.frame(DB)
  row.names(DB) <- 1:length(DB[,1])
  colnames(DB)[1] <- "Sequence_full"
  DB$Sequence_full <- as.character(DB$Sequence_full)
  DB$Accessions <- str_split_fixed(DB$Accessions, "\\|", 3)[,2]
  return(DB)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
