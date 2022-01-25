
source("R/Functions.R")

# SC: SN_Corrected

ScanInfo <- fread("data/ScanInfo_MQ.txt")

All <- fread("tables/PSMs_Annotated.txt")[!(Sample %like% "K_|0ms")]

All <- All[Organism != "" & Organism != "Nonunique"]

All[,Identifier := str_remove(`Spectrum File`,"\\.raw")]%>%
  .[,Identifier := str_c(`Identifier`,`First Scan`,sep = "__")]%>%
  .[,ID := 1:dim(All)[1]]

All <- All[Confidence == "High",]%>%unique()

All[,Raw := str_split_fixed(Identifier,"__",2)[,1]]%>%
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

All[,Av12SN := sum(c(SN128N,SN129N,SN129C,SN130N,SN130C,
                    SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                  na.rm = TRUE)/12, by = ID]

All <- merge(All, ScanInfo, by = "Identifier")

# SCH: SN_Corrected Human

SCH <- All[Organism == "Human"]

# This part is for normalization. However, normalization doesn't increase number of PSMs <= 20% in CV too much. 
# So no normalization was conducted.

# Files <- SCH$Raw%>%unique()
# 
# SCH0 <- data.table()
# 
# for (i in 1:length(Files)) {
#   SCH1 <- SCH[Raw == Files[i]]
# 
#   x <- colSums(SCH1[,c(35,37:49)],na.rm = TRUE)%>%median()
# 
#   m <- SCH1[,c(35,37:49)]
# 
#   m <- sweep(m, 2, colSums(m,na.rm = TRUE), FUN="/")
#   m <- scale(m, center=FALSE, scale=colSums(m,na.rm = TRUE)/x)%>%data.table()
# 
#   SCH1[,c(35,37:49)] <- m
# 
#   SCH1[,mean := mean(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
#                       SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
#     .[,sd := sd(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
#                   SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
#     .[,CV := sd/mean * 100]%>%
#     .[,ID := 1:dim(SCH1)[1]]
# 
#   SCH0 <- rbind(SCH0,SCH1)
# }
# 
# x <- colSums(SCH[,c(35,37:49)],na.rm = TRUE)%>%median()
# 
# m <- SCH[,c(35,37:49)]
# 
# m <- sweep(m, 2, colSums(m,na.rm = TRUE), FUN="/")
# m <- scale(m, center=FALSE, scale=colSums(m,na.rm = TRUE)/x)%>%data.table()
# 
# SCH[,c(35,37:49)] <- m

SCH[,mean := mean(c(SN128N,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                    SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN128N,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]%>%
  .[,ID := 1:dim(SCH)[1]]

table(SCH$CV >= 20)

# SCH[Booster == "HEY" & Method %like% "SN"]%>%
# ggplot(aes(x = CV,fill = Method)) + geom_density(alpha = 0.2) + 
#   facet_grid(Amount ~ Ratio) + xlim(0,100)

SCH[Booster == "HEY" & Method %like% "SN_Corrected" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = AGC,fill =(CV <= 20))) + 
  geom_bar(stat = "count",position = position_stack())+ 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.2, size = 4)+ 
  facet_grid(Amount ~ Ratio)

Index <- apply(SCH, 1, function(x) sum(is.na(x)))

SCH2 <- SCH[Index <= 4]


fwrite(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Corrected" & Amount == "50pg"],"tables/Fig1de_Fig.S3a.txt")

pdf("figures/CV_PSMs_Av12.pdf",width = 10,height = 5)

ggplot(SCH2[Av12SN >= 1 & Booster == "HEY" & Amount == "50pg" & is.na(CV) == FALSE], 
       aes(x = Method, fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~Ratio,scales = "free") + theme_prism() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))


ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & AGC == "AGC300%" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(Method~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")


ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Corrected" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n") + ggtitle("SN_Corrected")

ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Raw" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")+ ggtitle("SN_Raw")

ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "Int_Corrected" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")+ ggtitle("Int_Corrected")

ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "Int_Raw" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")+ ggtitle("Int_Raw")

ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Corrected" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV), color = log10(RawOvFtT))) + 
  geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")+ ggtitle("RawOvFtT")
  

dev.off()

pdf("figures/CV_PSMs_Av12_v2.pdf",width = 10,height = 10)

p1 <- ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Corrected" & Amount == "50pg"], 
       aes(x = log10(RawOvFtT), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(Ratio~AGC) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n") + ggtitle("SN_Corrected")

p2 <- ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "SN_Corrected" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(Ratio~AGC) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n") + ggtitle("SN_Corrected")


p3 <- ggplot(SCH2[Av12SN >= 1 & Booster %like% "No|HEY" & Method == "Int_Corrected" & Amount == "50pg"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue", alpha = 0.3)+
  facet_grid(Ratio~AGC) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism() + 
  theme(legend.position = "none")+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n")+ ggtitle("Int_Corrected")

multiplot(p1, p2, p3, cols=3)

dev.off()
