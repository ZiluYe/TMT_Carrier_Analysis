
library(OrgMassSpecR)
library(data.table)
library(magrittr)
library(tidyverse)
library(seqinr)
library(ggpubr)

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

Proteins <- fread("data/proteinGroups.txt")

Proteins <- Proteins[`Potential contaminant` == "" & Reverse == ""]%>%
  .[,Accessions := str_split_fixed(`Majority protein IDs`,"\\|",3)[,2]]


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

Fasta <- FASTAtoDF("data/SProt_Human_reviewed_March_2019.fasta")

library(plyr)

Fasta$count <- aaply(Fasta$Sequence_full,1,PepNumber)

Proteins2 <- merge(Proteins,Fasta, by = "Accessions")%>%data.table()

ggplot(Proteins2,aes(x = log10(iBAQ), y = log10(Intensity / count))) + geom_point()+ 
  stat_cor(method = "pearson",r.digits = 3, label.sep = "\n")

