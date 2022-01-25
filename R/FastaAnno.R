
FASTAtoDF <- function(Fasta){
  DB <- seqinr::read.fasta(Fasta, as.string = TRUE, seqtype = "AA")
  DB <- as.data.frame(do.call(rbind, DB))
  DB$Accessions <- row.names.data.frame(DB)
  row.names(DB) <- 1:length(DB[,1])
  colnames(DB)[1] <- "Sequence_full"
  DB$Accessions <- str_split_fixed(DB$Accessions, "\\|", 3)[,2]
  return(data.table(DB))
}


Ecoli <- FASTAtoDF("Ecoli_UP000000625_83333.fasta")%>%.[,Organism := "Ecoli"]
Human <- FASTAtoDF("Human_UP000005640_9606.fasta")%>%.[,Organism := "Human"]
Yeast <- FASTAtoDF("Yeast_UP000002311_559292.fasta")%>%.[,Organism := "Yeast"]

HEY <- rbind(Human[,.(Accessions,Organism)],Ecoli[,.(Accessions,Organism)],Yeast[,.(Accessions,Organism)])

#fwrite(HEY,"tables/FastaAnno.txt")

PSMs <- fread("tables/PSMs.txt")

PSMs[,ID := 1:dim(PSMs)[1]]

Uniprot <- splitstackshape::cSplit(PSMs[,.(`Master Protein Accessions`)], 1, sep = ";", type.convert = FALSE)

Uniprot[,ID := 1:dim(PSMs)[1]]

Uniprot <- melt(Uniprot,id.vars = "ID")

Uniprot <- merge(Uniprot,HEY,by.x = "value",by.y = "Accessions")

Uniprot <- Uniprot[,.(ID,Organism)]%>%unique()

Uniprot <- dcast(Uniprot,ID ~ Organism)

Uniprot[is.na(Uniprot) == TRUE] <- "No"

Uniprot[,Organism := "Nonunique"]%>%
  .[Human != "No" & Yeast == "No" & Ecoli == "No",Organism := "Human"]%>%
  .[Human == "No" & Yeast != "No" & Ecoli == "No",Organism := "Yeast"]%>%
  .[Human == "No" & Yeast == "No" & Ecoli != "No",Organism := "Ecoli"]


PSMs <- merge(PSMs,Uniprot[,.(ID,Organism)],by = "ID", all.x = TRUE)

fwrite(PSMs,"tables/PSMs_Annotated.txt")
