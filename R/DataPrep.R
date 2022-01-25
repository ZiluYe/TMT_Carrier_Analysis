
library(tidyverse)
library(data.table)
library(magrittr)

Files <- dir("data/Results/")

Files <- str_replace_all(Files,"2020","data/Results/2020")

# Index <- str_detect(Files,"PSMs|QuanSpectra|Proteins")
# 
# file.remove(Files[Index == FALSE])

#create a list of the files from your target directory
# F1 <- Files[str_detect(Files,"[12]_Proteins")]




Proteins_Int_Raw <- Bind(Files[str_detect(Files,"[12s%]_Proteins")])%>%.[,Method := "Int_Raw"]

Proteins_Int_Corrected <- Bind(Files[str_detect(Files,"1\\)_Proteins")])%>%.[,Method := "Int_Corrected"]

Proteins_SN_Corrected <- Bind(Files[str_detect(Files,"2\\)_Proteins")])%>%.[,Method := "SN_Corrected"]

Proteins_SN_Raw <- Bind(Files[str_detect(Files,"3\\)_Proteins")])%>%.[,Method := "SN_Raw"]

Proteins <- rbind(Proteins_Int_Raw,Proteins_Int_Corrected,Proteins_SN_Raw,Proteins_SN_Corrected)

colnames(Proteins) <- str_remove_all(colnames(Proteins),"F7 ")

colnames(Proteins)[32:47] <- str_c(str_split_fixed(colnames(Proteins)[32:47], " ", 7)[,1],
                                   str_split_fixed(colnames(Proteins)[32:47], " ", 7)[,2],
                                   str_split_fixed(colnames(Proteins)[32:47], " ", 7)[,6],sep = " ")

Proteins[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_Proteins.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]


fwrite(Proteins,"tables/Proteins.txt")



PSMs_Int_Raw <- Bind(Files[str_detect(Files,"[12s%]_PSMs")])%>%.[,Method := "Int_Raw"]

PSMs_Int_Corrected <- Bind(Files[str_detect(Files,"1\\)_PSMs")])%>%.[,Method := "Int_Corrected"]

PSMs_SN_Corrected <- Bind(Files[str_detect(Files,"2\\)_PSMs")])%>%.[,Method := "SN_Corrected"]

PSMs_SN_Raw <- Bind(Files[str_detect(Files,"3\\)_PSMs")])%>%.[,Method := "SN_Raw"]

PSMs <- rbind(PSMs_Int_Raw,PSMs_Int_Corrected,PSMs_SN_Raw,PSMs_SN_Corrected)

PSMs[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_Proteins.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]

fwrite(PSMs,"tables/PSMs.txt")


QuanSpectra_Int_Raw <- Bind(Files[str_detect(Files,"[12s%]_QuanSpectra")])%>%.[,Method := "Int_Raw"]

QuanSpectra_Int_Corrected <- Bind(Files[str_detect(Files,"1\\)_QuanSpectra")])%>%.[,Method := "Int_Corrected"]

QuanSpectra_SN_Corrected <- Bind(Files[str_detect(Files,"2\\)_QuanSpectra")])%>%.[,Method := "SN_Corrected"]

QuanSpectra_SN_Raw <- Bind(Files[str_detect(Files,"3\\)_QuanSpectra")])%>%.[,Method := "SN_Raw"]

QuanSpectra <- rbind(QuanSpectra_Int_Raw,QuanSpectra_Int_Corrected,QuanSpectra_SN_Raw,QuanSpectra_SN_Corrected)

QuanSpectra[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_Proteins.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]

fwrite(QuanSpectra,"tables/QuanSpectra.txt")







Files <- dir("data/HEY_Minora/")

Files <- str_replace_all(Files,"2020","data/HEY_Minora/2020")


Proteins <- Bind(Files[str_detect(Files,"[12s%]_Proteins")])

colnames(Proteins)[16] <- "TMT_Minora"

Proteins[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_Proteins.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]


fwrite(Proteins,"tables/Proteins_Minora.txt")


