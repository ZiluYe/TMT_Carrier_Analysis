
source("R/Functions.R")

# All: SN_Corrected

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
#   x <- colSums(SCH1[,c(17,19:31)],na.rm = TRUE)%>%median()
# 
#   m <- SCH1[,c(17,19:31)]
# 
#   m <- sweep(m, 2, colSums(m,na.rm = TRUE), FUN="/")
#   m <- scale(m, center=FALSE, scale=colSums(m,na.rm = TRUE)/x)%>%data.table()
# 
#   SCH1[,c(17,19:31)] <- m
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
# x <- colSums(SCH[,c(17,19:31)],na.rm = TRUE)%>%median()
# 
# m <- SCH[,c(17,19:31)]
# 
# m <- sweep(m, 2, colSums(m,na.rm = TRUE), FUN="/")
# m <- scale(m, center=FALSE, scale=colSums(m,na.rm = TRUE)/x)%>%data.table()
# 
# SCH[,c(17,19:31)] <- m

SCH[,mean := mean(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                    SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]%>%
  .[,ID := 1:dim(SCH)[1]]

table(SCH$CV >= 20)

# SCH[Booster == "HEY" & Method %like% "SN"]%>%
#   ggplot(aes(x = CV,fill = Method)) + geom_density(alpha = 0.2) + 
#   facet_grid(Amount ~ Ratio) + xlim(0,100)
# 
# SCH[Booster == "HEY" & Method %like% "SN_Corrected" & is.na(CV) == FALSE]%>%
#   ggplot(aes(x = AGC,fill =(CV <= 20))) + 
#   geom_bar(stat = "count",position = position_stack())+ 
#   geom_text(stat='count', aes(label=..count..), position = position_stack(),
#             vjust=0.2, size = 4)+ 
#   facet_grid(Amount ~ Ratio)
# 
# ggplot(SCH[Amount == "100pg"], aes(x = Ratio, y = CV, color = Ratio)) + 
#   geom_boxplot(outlier.shape = NA) +
#   # geom_jitter(size=1,width = 0.2)+
#   facet_grid(AGC ~Method, scales = "free_y") + ylim(0,100)+ 
#   geom_hline(yintercept = 20, size = 1, linetype = "dashed")

# ggplot(SCH0[Av14SN >= 1 & Booster == "HEY" & ID %% 10 == 1], 
#        aes(x = log2(Av14SN), y = log2(CV), color = log10(RawOvFtT))) + 
#   geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
#   facet_grid(Amount~Ratio,scales = "free") + 
#   geom_hline(yintercept = log2(20),size = 1, linetype = "dashed")
# 
# ggplot(SCH0[Av14SN >= 1 & Booster == "HEY" & ID %% 10 == 1], 
#        aes(color = log2(Av14SN), y = log2(CV), x = log10(RawOvFtT))) + 
#   geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
#   facet_grid(Amount~Ratio,scales = "free") + 
#   geom_hline(yintercept = log2(20),size = 1, linetype = "dashed")



pdf("figures/CV_Proteins.pdf",width = 10,height = 5)

ggplot(SCH[Booster == "HEY" & Amount == "200pg" & is.na(CV) == FALSE], 
       aes(x = Method, fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 3)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~Ratio,scales = "free") + theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=6))

ggplot(SCH[Av14SN >= 1 & Booster == "HEY" & Method == "SN_Corrected" & Amount == "200pg"], 
       aes(x = log2(Av14SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue")+scale_color_gradient(low = "white",high = "red")+
  facet_grid(AGC~Ratio,scales = "free") + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed")

dev.off()
