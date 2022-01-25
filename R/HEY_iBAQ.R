
source("R/Functions.R")

Fasta <- fread("tables/Fasta.txt")

P1 <- fread("tables/Proteins.txt")%>%.[Master == "IsMasterProtein"]

P1 <- P1[,c(4:5,16:31,50:51)]

P1 <- melt(P1,id.vars = c("Accession","Description","Sample","Method"))

P1[,variable := str_split_fixed(variable," ",3)[,2]]

P1 <- dcast(P1,Accession + Description + Sample ~ variable + Method, value.var = "value")

P2 <- fread("tables/Proteins_Minora.txt")%>%.[Master == "IsMasterProtein"]

Proteins <- merge(P1,P2[,.(Accession,Description,Sample,TMT_Minora)],
                  by = c("Accession","Description","Sample"))

Human <- melt(Proteins,id.vars = c("Accession","Description","Sample"))

Human <- merge(Human,Fasta,by = "Accession")

Human <- Human[,value := value / count]

Human <- dcast(Human,Accession + Description + Sample ~ variable, value.var = "value")

# res <- cor(Human[,4:23], method = "spearman", use = "complete.obs")
# 
# corrplot.mixed(res, upper = "ellipse", lower = "number", tl.pos = "lt")



## To calculate all


Files <- Human$Sample%>%unique()

HALL <- data.table()

for (i in 1:length(Files)) {
  H0 <- Human[Sample == Files[i]]
  H0 <- cor(H0[,4:dim(H0)[2]], method = "spearman", use = "complete.obs")%>%data.table(keep.rownames = TRUE)
  H0 <- melt(H0,id.vars = "rn")
  H0[,Sample := Files[i]]
  HALL <- rbind(H0,HALL)
}


HALL <- HALL[!(Sample %like% "K_|0ms")]

HM <- HALL[(rn %like% "Minora") & !(variable %like% "Minora")]


HM[,Booster := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Amount := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",3)[,2]]%>%
  .[,Amount := factor(Amount, levels = c("50pg","100pg","200pg"))]%>%
  .[,AGC := str_split_fixed(Sample, "pg_",2)[,2]]%>%
  .[,AGC := str_split_fixed(AGC, "_",2)[,1]]%>%
  .[,AGC := factor(AGC, levels = c("AGC50%","AGC300%"))]%>%
  .[,Rep := str_split_fixed(Sample, "%_",2)[,2]]%>%
  .[,y1 := str_split_fixed(variable,"_",3)[,1]]%>%
  .[,y2 := str_split_fixed(variable,"_",3)[,2]]%>%
  .[,y3 := str_split_fixed(variable,"_",3)[,3]]

fwrite(HM[Booster == "HEY" & (variable %like% "126|127|128") & y1 != "127C"],"tables/Fig4.txt")

pdf("figures/MS1vsReporterIons_v2.pdf",width = 9,height = 5)

# HM[Booster == "HEY" & (variable %like% "12") & y3 == "Raw"]%>%
# ggplot(aes(x = y1, y = value, color = Ratio)) + 
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(AGC ~ y2)  + ylab("Spearman R")+
#   scale_color_brewer(palette="Set2",name="Ratio") + theme_prism()
# 
# HM[Booster == "HEY" & (variable %like% "12")]%>%
#   ggplot(aes(x = y1, y = value, color = Ratio)) + 
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(AGC ~ y3 + y2)  + ylab("Spearman R")+
#   scale_color_brewer(palette="Set2",name="Ratio")
# 
# HM[Booster == "HEY" & (variable %like% "12") & y1 != "127C"]%>%
#   ggplot(aes(x = y1, y = value, color = Ratio)) + 
#   geom_boxplot(outlier.shape = NA) +
#   facet_grid(AGC ~ y3 + y2)  + ylab("Spearman R")+ 
#   theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=6))+
#   scale_color_brewer(palette="Set2",name="Ratio")

HM[Booster == "HEY" & (variable %like% "126|127|128") & y1 != "127C"]%>%
  ggplot(aes(x = y1, y = value, color = Ratio)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_grid(AGC ~ y3 + y2)  + ylab("Spearman R")+ 
  theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=6))+
  scale_color_brewer(palette="Set2",name="Ratio") + theme_prism(base_size = 10,axis_text_angle = 45)+
  ylim(0.75,1)

HM[Booster == "HEY" & (variable %like% "126|127|128") & y1 != "127C"]%>%
  ggplot(aes(x = Ratio, y = value, color = Ratio)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(size=0.3,width = 0.2) +
  facet_grid(AGC ~ y3 + y2 + y1,scales = "free")  + ylab("Spearman R")+ 
  theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=6))+
  scale_color_brewer(palette="Set2",name="Ratio") + theme_prism(base_size = 10,axis_text_angle = 45)+
  ylim(0.75,1)

dev.off()

# HM[(variable %like% "126") & Booster %like% "H"]%>%
#   ggplot(aes(x = Booster, y = value, color = Ratio)) + 
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(size=0.1,alpha = 0.7,aes(x = Booster, y = value, color = Ratio))+
#   facet_grid(AGC + Amount ~ variable)  + ylab("Spearman R")+
#   scale_color_brewer(palette="Set2",name="With128C")

#xx <- melt(HALL,id.vars = c("rn","Raw"))

HALL <- HALL[rn != "TMT_Minora" & variable != "TMT_Minora"]

HALL[,Booster := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Amount := str_split_fixed(Sample, "TMT_",2)[,2]]%>%
  .[,Amount := str_split_fixed(Amount, "_",3)[,2]]%>%
  .[,Amount := factor(Amount, levels = c("50pg","100pg","200pg"))]%>%
  .[,AGC := str_split_fixed(Sample, "pg_",2)[,2]]%>%
  .[,AGC := str_split_fixed(AGC, "_",2)[,1]]%>%
  .[,AGC := factor(AGC, levels = c("AGC50%","AGC300%"))]%>%
  .[,Rep := str_split_fixed(Sample, "%_",2)[,2]]%>%
  .[,x1 := str_split_fixed(rn,"_",3)[,1]]%>%
  .[,x2 := str_split_fixed(rn,"_",3)[,2]]%>%
  .[,x3 := str_split_fixed(rn,"_",3)[,3]]%>%
  .[,y1 := str_split_fixed(variable,"_",3)[,1]]%>%
  .[,y2 := str_split_fixed(variable,"_",3)[,2]]%>%
  .[,y3 := str_split_fixed(variable,"_",3)[,3]]%>%
  .[,M23 := str_c(x2,x3,y2,y3,sep = "_")]%>%
  .[x2 == "SN", M23 := str_c(y2,y3,x2,x3,sep = "_")]%>%
  .[,variable := as.character(variable)]


old <- HALL$rn%>%unique()

new <- 1:64

HA <- HALL[,I1 := factor(rn, levels = old, labels = new)%>%as.numeric()]%>%
  .[,I2 := factor(variable, levels = old, labels = new)%>%as.numeric()]

HA <- HA[I1 > I2]

fwrite(HA[!(x1 %like% "126|127C") & !(y1 %like% "126|127C") & x2 != y2 & Amount == "100pg" & x3 == y3 & Booster == "HEY"],"tables/Fig.S6.txt")

pdf("figures/IntvsSN.pdf",width = 10,height = 5)

HA[!(x1 %like% "126|127C") & !(y1 %like% "126|127C") & x2 != y2 & Amount == "100pg" & x3 == y3 & Booster == "HEY"]%>%
  ggplot(aes(x = Ratio, y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.1,width = 0.2,aes(color = (x1 == "128C") | (y1 == "128C")),alpha = 0.7)+
  facet_grid(AGC ~ M23)  + ylab("Spearman R")+
  scale_color_brewer(palette="Set2",name="With128C") + theme_prism()

HA[!(x1 %like% "126|127C") & !(y1 %like% "126|127C") & x2 != y2 & Amount == "100pg" & x3 == y3 & Booster == "HEY" & AGC == "AGC300%"]%>%
  ggplot(aes(x = Ratio, y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.1,width = 0.2,aes(color = (x1 == "128C") | (y1 == "128C")),alpha = 0.7)+
  facet_grid(. ~ M23)  + ylab("Spearman R")+
  scale_color_brewer(palette="Set2",name="With128C") + theme_prism(base_size = 10)


dev.off()

HA[,Channel := "Others"]%>%
  .[(x1 == "126") | (y1 == "126"),Channel := "126"]%>%
  .[(x1 == "127C") | (y1 == "127C"),Channel := "127C"]%>%
  .[(x1 == "128C") | (y1 == "128C"),Channel := "128C"]


HA[x2 != y2 & Amount == "100pg" & x3 == y3 & Booster == "HEY"]%>%
  ggplot(aes(x = Ratio, y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.1,width = 0.2,aes(color = Channel),alpha = 0.7)+
  facet_grid(AGC ~ M23)  + ylab("Spearman R")+
  scale_color_brewer(palette="Set2",name="With128C")
