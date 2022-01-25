
source("R/Functions.R")

Fasta <- fread("tables/Fasta.txt")

H126 <- fread("tables/H126.txt")

LFQ <- H126[Sample %like% "LFQ"]


H126 <- H126[!(Sample %like% "15K") & Rep == "1"]

ggplot(H126, aes(x = Method, fill = factor(Amount))) + geom_bar(position = position_dodge())

H126 <- merge(H126,Fasta,by = "Accession")

H126 <- H126[,Abundance := Abundance / count]

H2 <- dcast(H126,Accession + Amount ~ Method + Rep,value.var = "Abundance")

H2 <- H2[Amount == "100ng",
         .(Accession,Amount,Int_Raw_1,SN_Raw_1,LFQ_Int_1,TMT_Int_1)]

colnames(H2) <- str_remove_all(colnames(H2),"_1")

res <- cor(H2[,3:dim(H2)[2]], method = "spearman", use = "complete.obs")

corrplot.mixed(res, upper = "ellipse", lower = "number", tl.pos = "lt")



LFQ[,Res := "60K"]%>%
  .[Sample %like% "15K", Res := "15K"]

LFQ <- merge(LFQ,Fasta,by = "Accession")

LFQ <- LFQ[,Abundance := Abundance / count]

H2 <- dcast(LFQ,Accession + Amount ~ Method + Res+ Rep,value.var = "Abundance")


res <- cor(H2[Amount == "100ng",3:dim(H2)[2]], method = "spearman", use = "complete.obs")

corrplot.mixed(res, upper = "ellipse", lower = "number", tl.pos = "lt")
