
source("R/Functions.R")


Fasta <- FASTAtoDF("Human_UP000005640_9606.fasta")

Fasta$count <- aaply(Fasta$Sequence_full,1,PepNumber)

Fasta <- Fasta[,c("Accessions","count")]%>%data.table()

colnames(Fasta)[1] <- "Accession"

fwrite(Fasta, "tables/Fasta.txt")

All <- fread("tables/Proteins.txt")

All <- All[Master == "IsMasterProtein"]%>%
  .[,Human := str_detect(`Description`, "Homo sapiens")]%>%
  .[,Yeast := str_detect(`Description`, "Saccharomyces cerevisiae")]%>%
  .[,Ecoli := str_detect(`Description`, "Escherichia coli")]%>%
  .[,Organism := "Nonunique"]%>%
  .[Human == TRUE & Yeast == FALSE & Ecoli == FALSE,Organism := "Human"]%>%
  .[Human == FALSE & Yeast == TRUE & Ecoli == FALSE,Organism := "Yeast"]%>%
  .[Human == FALSE & Yeast == FALSE & Ecoli == TRUE,Organism := "Ecoli"]

All <- All[!(Sample %like% "K_|0ms")]

All <- All[Organism != "" & Organism != "Nonunique"]

All[,Raw := Sample]%>%
  .[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Amount := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",3)[,2]]%>%
  .[,Amount := factor(Amount, levels = c("50pg","100pg","200pg"))]%>%
  .[,AGC := str_split_fixed(`Raw`, "pg_",2)[,2]]%>%
  .[,AGC := str_split_fixed(AGC, "_",2)[,1]]%>%
  .[,AGC := factor(AGC, levels = c("AGC50%","AGC300%"))]%>%
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]%>%
  .[,Organism := factor(Organism, levels = c("Human","Ecoli","Yeast"))]

colnames(All) <- str_replace_all(colnames(All), "Abundance ", "SN")
colnames(All) <- str_remove_all(colnames(All), " Sample")

All[,ID := 1:dim(All)[1]]

All[,Av14SN := sum(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,
                     SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                   na.rm = TRUE)/14, by = ID]




Human <- All[Organism == "Human"]


Human <- Human[,.(Accession,SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N,Method,Raw)]%>%unique()

colnames(Human) <- str_remove_all(colnames(Human),"SN")

H1 <- melt(Human,id.vars = c("Accession","Method","Raw"))

H1 <- merge(H1,Fasta,by = "Accession")

H1 <- H1[,value := value / count]

H1 <- H1[Method %like% "Raw"]%>%.[,Method := str_remove_all(Method,"_Raw")]

H1[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Amount := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",3)[,2]]%>%
  .[,Amount := factor(Amount, levels = c("50pg","100pg","200pg"))]%>%
  .[,AGC := str_split_fixed(`Raw`, "pg_",2)[,2]]%>%
  .[,AGC := str_split_fixed(AGC, "_",2)[,1]]%>%
  .[,AGC := factor(AGC, levels = c("AGC50%","AGC300%"))]%>%
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]

H2 <- dcast(H1[Ratio == "98"],Accession + Raw ~ Method + variable,value.var = "value")


res <- cor(H2[,3:dim(H2)[2]], method = "spearman", use = "complete.obs")

corrplot.mixed(res, upper = "ellipse", lower = "number", tl.pos = "lt", cl.lim = c(0.8, 1), is.corr = FALSE)



## To calculate all

H1 <- melt(Human,id.vars = c("Accession","Method","Raw"))

H1 <- merge(H1,Fasta,by = "Accession")

H1 <- H1[,value := value / count]

H2 <- dcast(H1,Accession + Raw ~ Method + variable,value.var = "value")

Files <- H2$Raw%>%unique()

HALL <- data.table()

for (i in 1:length(Files)) {
  H0 <- H2[Raw == Files[i]]
  H0 <- cor(H0[,3:dim(H0)[2]], method = "spearman", use = "complete.obs")%>%data.table(keep.rownames = TRUE)
  H0 <- melt(H0,id.vars = "rn")
  H0[,Raw := Files[i]]
  HALL <- rbind(H0,HALL)
}

#xx <- melt(HALL,id.vars = c("rn","Raw"))

HALL[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Amount := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",3)[,2]]%>%
  .[,Amount := factor(Amount, levels = c("50pg","100pg","200pg"))]%>%
  .[,AGC := str_split_fixed(`Raw`, "pg_",2)[,2]]%>%
  .[,AGC := str_split_fixed(AGC, "_",2)[,1]]%>%
  .[,AGC := factor(AGC, levels = c("AGC50%","AGC300%"))]%>%
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]%>%
  .[,x1 := str_split_fixed(rn,"_",3)[,1]]%>%
  .[,x2 := str_split_fixed(rn,"_",3)[,2]]%>%
  .[,x3 := str_split_fixed(rn,"_",3)[,3]]%>%
  .[,y1 := str_split_fixed(variable,"_",3)[,1]]%>%
  .[,y2 := str_split_fixed(variable,"_",3)[,2]]%>%
  .[,y3 := str_split_fixed(variable,"_",3)[,3]]%>%
  .[,M12 := str_c(x1,x2,y1,y2,sep = "_")]%>%
  .[,variable := as.character(variable)]


old <- HALL$rn%>%unique()

new <- 1:56

HA <- HALL[,I1 := factor(rn, levels = old, labels = new)%>%as.numeric()]%>%
  .[,I2 := factor(variable, levels = old, labels = new)%>%as.numeric()]

HA <- HA[I1 > I2]


HA <- HALL[x1 == "Int" & y1 == "SN"]

# ggplot(HA[Amount == "100pg" & x2 == "Raw" & y2 == "Raw"], aes(x = M12, y = value,fill= M12)) + 
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(size=0.1,width = 0.2,aes(color = (x3 == "128C") | (y3 == "128C") ))+
#   facet_grid(AGC ~Ratio, scales = "free_y") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=5))


ggplot(HA[Amount == "100pg" & x2 == y2], aes(x = Ratio, y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.1,width = 0.2,aes(color = (x3 == "128C") | (y3 == "128C")),alpha = 0.7)+
  facet_grid(AGC ~ M12)  + ylab("Spearman R")+
  scale_color_brewer(palette="Set2",name="With128C")
