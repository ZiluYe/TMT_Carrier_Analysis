
source("R/Functions.R")

library(ggprism)


Files <- dir("data/TurboTMT/")

Files <- Files[str_detect(Files, "SA")]

Files <- str_replace_all(Files,"2021","data/TurboTMT/2021")


Proteins_SN <- Bind(Files[str_detect(Files,"01_Proteins")])%>%.[,Method := "SN"]

Proteins_Int <- Bind(Files[str_detect(Files,"1\\)_Proteins")])%>%.[,Method := "Int"]

Proteins <- rbind(Proteins_Int,Proteins_SN)

colnames(Proteins) <- str_remove_all(colnames(Proteins),"F3 ")


Proteins[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_Proteins.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]%>%
  .[,Sample := str_replace_all(Sample,"_15K_turboTMT-OFF_NCE29", "_15K_turboTMT-ON_NCE29")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE32", "_15K_turboTMT-ON_NCE32")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE35", "_15K_turboTMT-ON_NCE35")]

Proteins[,Condition := str_split_fixed(Sample,"500ng_",2)[,2]]



PSMs_SN <- Bind(Files[str_detect(Files,"01_PSMs")])%>%.[,Method := "SN"]

PSMs_Int <- Bind(Files[str_detect(Files,"1\\)_PSMs")])%>%.[,Method := "Int"]

PSMs <- rbind(PSMs_Int,PSMs_SN)

colnames(PSMs) <- str_remove_all(colnames(PSMs),"F3 ")


PSMs[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_PSMs.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]%>%
  .[,Sample := str_replace_all(Sample,"_15K_turboTMT-OFF_NCE29", "_15K_turboTMT-ON_NCE29")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE32", "_15K_turboTMT-ON_NCE32")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE35", "_15K_turboTMT-ON_NCE35")]

PSMs[,Condition := str_split_fixed(Sample,"500ng_",2)[,2]]




MSMS_SN <- Bind(Files[str_detect(Files,"01_MSMS")])%>%.[,Method := "SN"]

MSMS_Int <- Bind(Files[str_detect(Files,"1\\)_MSMS")])%>%.[,Method := "Int"]

MSMS <- rbind(MSMS_Int,MSMS_SN)

MSMS[,Sample := str_split_fixed(Sample,"/",3)[,3]]%>%
  .[,Sample := str_remove_all(Sample,"_MSMSSpectrumInfo.txt")]%>%
  .[,Sample := str_remove_all(Sample,"-\\([123]\\)")]%>%
  .[,Sample := str_replace_all(Sample,"_15K_turboTMT-OFF_NCE29", "_15K_turboTMT-ON_NCE29")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE32", "_15K_turboTMT-ON_NCE32")]%>%
  .[,Sample := str_replace_all(Sample,"_60K_turboTMT-OFF_NCE35", "_15K_turboTMT-ON_NCE35")]


MSMS[,Condition := str_split_fixed(Sample,"500ng_",2)[,2]]

P1 <- Proteins[Master == "IsMasterProtein" & Method == "SN",.(Sample,Condition,Accession)]%>%
  .[,Category := "Proteins"]%>%unique()

P2 <- PSMs[Confidence == "High" & Method == "SN",.(Sample,Condition,`Master Protein Accessions`)]%>%
  .[,Category := "PSMs"]

P22 <- MSMS[Method == "SN",.(Sample,Condition,`First Scan`)]%>%
  .[,Category := "MSMS"]

Numbers <- rbind(P1,P2,P22,use.names=FALSE)

Numbers[,Condition := str_remove_all(Condition,"_01")]%>%
  .[,Condition := str_remove_all(Condition,"TMT")]%>%
  .[,Category := factor(Category, levels = c("Proteins", "PSMs", "MSMS"))]


# ggplot(Numbers, aes(x = Condition, fill = Category)) +
#   geom_bar(position = position_dodge(),alpha = 0.9) +
#   scale_fill_brewer(palette = "Set2")+
#   geom_text(stat='count', aes(label=..count..), position = position_dodge(width = 0.9),
#             vjust=-0.1, size = 4) + theme_prism()


## Quant

P3 <- PSMs[Confidence == "High"]%>%
  .[,Condition := str_remove_all(Condition,"_01")]%>%
  .[,Condition := str_remove_all(Condition,"TMT")]

Index <- apply(P3[,32:47], 1, function(x) sum(is.na(x)))

P3[,MV := Index]

p <- P3[Method == "SN"]

# ggplot(p, aes(x = Condition, fill = MV >= 8)) +
#   geom_bar(position = position_dodge(),alpha = 0.9) +
#   scale_fill_brewer(palette = "Set2")+
#   geom_text(stat='count', aes(label=..count..), position = position_dodge(width = 0.9),
#             vjust=-0.1, size = 4)


P3[,Precursor := str_c(`Annotated Sequence`,Charge,sep = "_")]

P3 <- P3[,rank := rank(-`Abundance 126`),by=.(Sample,Method,Precursor)]

P3 <- P3[rank == 1]

# P3%>%
#   ggplot(aes(x=Condition, y=log10(`Abundance 126`),group = Precursor,color = Condition)) +
#   geom_line(color = "grey") +
#   geom_point() + facet_grid(Method~.,scales = "free")



P4 <- P3[,c(54,55,57,32:47)]%>%unique()

P4[ ,`126/134N` := `Abundance 126`/`Abundance 134N`]%>%
  .[ ,`127N/134N` := `Abundance 127N`/`Abundance 134N`]%>%
  .[ ,`127C/134N` := `Abundance 127C`/`Abundance 134N`]%>%
  .[ ,`128N/134N` := `Abundance 128N`/`Abundance 134N`]%>%
  .[ ,`128C/134N` := `Abundance 128C`/`Abundance 134N`]%>%
  .[ ,`129N/134N` := `Abundance 129N`/`Abundance 134N`]%>%
  .[ ,`129C/134N` := `Abundance 129C`/`Abundance 134N`]%>%
  .[ ,`130N/134N` := `Abundance 130N`/`Abundance 134N`]%>%
  .[ ,`130C/134N` := `Abundance 130C`/`Abundance 134N`]%>%
  .[ ,`131N/134N` := `Abundance 131N`/`Abundance 134N`]%>%
  .[ ,`131C/134N` := `Abundance 131C`/`Abundance 134N`]%>%
  .[ ,`132N/134N` := `Abundance 132N`/`Abundance 134N`]%>%
  .[ ,`132C/134N` := `Abundance 132C`/`Abundance 134N`]%>%
  .[ ,`133N/134N` := `Abundance 133N`/`Abundance 134N`]%>%
  .[ ,`133C/134N` := `Abundance 133C`/`Abundance 134N`]

P5 <- melt(P4[,c(1:3,20:34)],id.vars = c("Method","Condition","Precursor"))

P6 <- dcast(P5,Method + Precursor + variable ~ Condition, value.var = "value")[Method == "SN"]

cols <- colnames(P6)[4:dim(P6)[2]]

P6[,(cols) := lapply(.SD, log2),.SDcols = cols]

# cor(P6[,4:dim(P6)[2]], method = "pearson", use = "na.or.complete")%>%
#   corrplot.mixed(upper = "ellipse", lower = "number", tl.pos = "lt")


# cor(P6[variable %like% "127C",4:dim(P6)[2]], method = "pearson", use = "complete.obs")%>%
#   corrplot.mixed(upper = "ellipse", lower = "number", tl.pos = "lt")
