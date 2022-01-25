
source("R/Functions.R")


Files <- dir("data/H126/")

Files <- str_replace_all(Files,"2020","data/H126/2020")

LFQ <- Files[str_detect(Files,"_LFQ")]
LFQ <- LFQ[str_detect(LFQ,"_Proteins")]

Area <- Bind(LFQ[str_detect(LFQ,"[12]_Proteins")])%>%.[,Method := "LFQ_Area"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F5 na Sample`,Sample,Method)]

Int <- Bind(LFQ[str_detect(LFQ,"1\\)_Proteins")])%>%.[,Method := "LFQ_Int"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F5 na Sample`,Sample,Method)]


TMT <- Files[str_detect(Files,"_TMT")]
TMT <- TMT[str_detect(TMT,"_Proteins")]

TMT_Area <- Bind(TMT[str_detect(TMT,"[12]_Proteins")])%>%.[,Method := "TMT_Area"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 na Sample`,Sample,Method)]

TMT_Int <- Bind(TMT[str_detect(TMT,"5\\)_Proteins")])%>%.[,Method := "TMT_Int"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 na Sample`,Sample,Method)]

SN_Corrected <- Bind(TMT[str_detect(TMT,"1\\)_Proteins")])%>%.[,Method := "SN_Corrected"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 126 Sample`,Sample,Method)]

SN_Raw <- Bind(TMT[str_detect(TMT,"2\\)_Proteins")])%>%.[,Method := "SN_Raw"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 126 Sample`,Sample,Method)]

Int_Raw <- Bind(TMT[str_detect(TMT,"3\\)_Proteins")])%>%.[,Method := "Int_Raw"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 126 Sample`,Sample,Method)]

Int_Corrected <- Bind(TMT[str_detect(TMT,"4\\)_Proteins")])%>%.[,Method := "Int_Corrected"]%>%
  .[Master == "IsMasterProtein",.(Accession,`Abundance F11 126 Sample`,Sample,Method)]

H126 <- rbind(Area,Int,TMT_Area,TMT_Int,SN_Raw,SN_Corrected,Int_Raw,Int_Corrected,use.names = FALSE)

colnames(H126)[2] <- "Abundance"

H126[,Amount := str_split_fixed(Sample, "Hela_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",2)[,1]]%>%
  .[,Rep := "1"]%>%
  .[Sample %like% "_2",Rep := "2"]

fwrite(H126,"tables/H126.txt")
