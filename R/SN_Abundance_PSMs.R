
source("R/Functions.R")

# SC: SN_Corrected

ScanInfo <- fread("data/ScanInfo_MQ.txt")

SC <- fread("tables/PSMs_Annotated.txt")[Method == "SN_Corrected" & Sample %like% "K_"]

SC <- SC[Confidence == "High"]%>%unique()

SC[,Identifier := str_remove(`Spectrum File`,"\\.raw")]%>%
  .[,Identifier := str_c(`Identifier`,`First Scan`,sep = "__")]

SC[,Raw := Sample]%>%
  .[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Resolution := str_split_fixed(`Raw`, "100pg_",2)[,2]]%>%
  .[,Resolution := str_split_fixed(Resolution, "_",2)[,1]]%>%
  .[,Resolution := factor(Resolution, levels = c("30K","45K","60K","120K"))]%>%
  .[,NCE := str_split_fixed(`Raw`, "K_",2)[,2]]%>%
  .[,NCE := str_split_fixed(NCE, "_",2)[,1]]

colnames(SC) <- str_replace_all(colnames(SC), "Abundance ", "SN")

SC[,Av14SN := sum(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,
                    SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                  na.rm = TRUE)/14, by = ID]

SC <- merge(SC, ScanInfo, by = "Identifier")


# ggplot(ScanInfo[1:1000,],aes(x = log10(`Total ion current` * `Ion injection time` / 1000), y = log10(RawOvFtT))) +
#   geom_point()

# ggplot(SC,aes(x = log10(Av14SN), y = log10(RawOvFtT), color = NCE)) + geom_point()

ggplot(SC[Resolution == "60K"],aes(x = NCE, fill = NCE)) + geom_bar() + 
  facet_grid(. ~ Ratio)

ggplot(SC[NCE == "32%"],aes(x = Resolution, fill = Resolution)) + geom_bar() + 
  facet_grid(. ~ Ratio)



SC[,Precursor := str_c(`Annotated Sequence`,Charge,sep = "_")]

SC2 <- SC[,.(Raw,`First Scan`,Precursor,Organism,Av14SN,SN126,RawOvFtT,XCorr)]%>%unique()

SC2[,rank := rank(-Av14SN),by=.(Raw,Precursor)]

SC2 <- SC2[rank == 1]

SC2[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Booster := str_split_fixed(Booster, "[0-9]",2)[,1]]%>%
  .[,Ratio := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
  .[,Ratio := str_split_fixed(Ratio, "_",2)[,1]]%>%
  .[,Ratio := str_remove_all(Ratio, "[HEY]")]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]%>%
  .[,Resolution := str_split_fixed(`Raw`, "100pg_",2)[,2]]%>%
  .[,Resolution := str_split_fixed(Resolution, "_",2)[,1]]%>%
  .[,Resolution := factor(Resolution, levels = c("30K","45K","60K","120K"))]%>%
  .[,NCE := str_split_fixed(`Raw`, "K_",2)[,2]]%>%
  .[,NCE := str_split_fixed(NCE, "_",2)[,1]]

P1 <- SC2[,.(Precursor)]%>%unique()

P1[,ID := 1:dim(P1)[1]]%>%
  .[,Sequence := str_split_fixed(Precursor,"_",2)[,1]%>%toupper()]%>%
  .[,Length := str_length(Sequence)]

SC2 <- merge(SC2, P1, by = "Precursor")

SC2[Length <= 10,LengthRange := "6-10"]%>%
  .[Length > 10, LengthRange := "11-15"]%>%
  .[Length > 15, LengthRange := "16-20"]%>%
  .[Length > 20, LengthRange := ">20"]%>%
  .[,LengthRange := factor(LengthRange, levels = c("6-10","11-15","16-20",">20"))]

#fwrite(SC2, "tables/NCE_SC2.txt")

SC3 <- SC2[Resolution == "60K"]

SC3 <- dcast(SC3,Precursor + ID + Ratio ~ NCE, value.var = "Av14SN")

SC3[,`29%` := `29%`/`26%`]%>%
  .[,`32%` := `32%`/`26%`]%>%
  .[,`35%` := `35%`/`26%`]%>%
  .[,`38%` := `38%`/`26%`]%>%
  .[,`26%` := `26%`/`26%`]

SC3 <- melt(SC3,id.vars = c("Precursor","ID","Ratio"))



SC4 <- SC2[NCE == "32%"]

SC4 <- dcast(SC4,Precursor + ID + Ratio ~ Resolution, value.var = "Av14SN")

SC4[,`45K` := `45K`/`30K`]%>%
  .[,`60K` := `60K`/`30K`]%>%
  .[,`120K` := `120K`/`30K`]%>%
  .[,`30K` := `30K`/`30K`]

SC4 <- melt(SC4,id.vars = c("Precursor","ID","Ratio"))




SC2.2 <- copy(SC2)

T1 <- table(SC2.2[,.(ID)])%>%data.table()

colnames(T1)[1] <- "ID"

SC2.2 <- merge(SC2.2, T1[,ID := as.integer(ID)], by = "ID")

SC2.2[N == 48 & Resolution == "60K" & Av14SN >= 1 & ID %% 20 == 4]%>%
  ggplot(aes(x=NCE, y=log10(Av14SN),group = Precursor,color = NCE)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 

fwrite(SC2.2[N == 48 & Resolution == "60K" & Av14SN >= 1 & ID %% 20 == 4],"tables/Fig3b.txt")

pdf("figures/NCE_Lines.pdf",width = 10,height = 5)

SC2[Resolution == "60K" & Av14SN >= 1 & ID %% 10 == 1]%>%
  ggplot(aes(x=NCE, y=log10(Av14SN),group = Precursor,color = NCE)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 

SC2.2[N == 48 & Resolution == "60K" & Av14SN >= 1 & ID %% 20 == 4]%>%
  ggplot(aes(x=NCE, y=log10(Av14SN),group = Precursor,color = NCE)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 

SC2[Resolution == "60K" & Av14SN >= 1 & ID %% 10 == 1]%>%
  ggplot(aes(x=NCE, y=XCorr,group = Precursor,color = NCE)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 


dev.off()


pdf("figures/NCE_Box.pdf",width = 10,height = 2)

SC3%>%
  ggplot(aes(x=variable, y=log2(value),color = variable)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10) + ylim(-2,5)

dev.off()

fwrite(SC2[Resolution == "60K" & Av14SN >= 1],"tables/Fig.S5a.txt")

pdf("figures/XCorr_Box.pdf",width = 10,height = 5)

SC2[Resolution == "60K" & Av14SN >= 1]%>%
  ggplot(aes(x=NCE, y=XCorr,color = NCE)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10) + ylim(0,5)

dev.off()

fwrite(SC2.2[N == 48 & NCE == "32%" & Av14SN >= 1 & ID %% 10 == 1],"tables/Fig3d.txt")

pdf("figures/Resolution.pdf",width = 10,height = 5)

SC2[NCE == "32%" & Av14SN >= 1 & ID %% 10 == 1]%>%
  ggplot(aes(x=Resolution, y=log10(Av14SN),group = Precursor,color = Resolution)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 

SC2.2[N == 48 & NCE == "32%" & Av14SN >= 1 & ID %% 10 == 1]%>%
  ggplot(aes(x=Resolution, y=log10(Av14SN),group = Precursor,color = Resolution)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio) + theme_prism(base_size = 10) 

SC2[NCE == "32%" & Av14SN >= 1]%>%
  ggplot(aes(x=Resolution, y=XCorr,color = Resolution)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10) + ylim(0,5)

dev.off()

pdf("figures/Resolution_Box.pdf",width = 10,height = 2)

SC4%>%
  ggplot(aes(x=variable, y=log2(value),color = variable)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10) + ylim(-2,5)


dev.off()









SC2[Resolution == "60K" & Ratio == "42" & Av14SN >= 1 & ID %% 5 == 1]%>%
  ggplot(aes(x=NCE, y=log10(Av14SN),group = Precursor,color = NCE)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(.~LengthRange)

SC2[NCE == "32%" & Av14SN >= 1 & ID %% 10 == 1]%>%
  ggplot(aes(x=Resolution, y=log10(Av14SN),group = Precursor,color = Resolution)) +
  geom_line(color = "grey") +
  geom_point() + facet_grid(. ~ Ratio)




SC[,mean := mean(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                    SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]%>%
  .[,ID := 1:dim(SC)[1]]

table(SC$CV >= 20)

#fwrite(SC, "tables/NCE_SC.txt")

fwrite(SC[Resolution == "60K" & is.na(CV) == FALSE], "tables/Fig3a.txt")
fwrite(SC[NCE == "32%" & is.na(CV) == FALSE], "tables/Fig3c.txt")

pdf("figures/NCE_Res_Numbers.pdf",width = 10,height = 5)

SC[Resolution == "60K" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = NCE,fill =(CV <= 20))) + scale_fill_brewer(palette = "Set2")+ 
  geom_bar(stat = "count",position = position_stack())+ 
  facet_grid(. ~ Ratio) + theme_prism(base_size = 10)

SC[NCE == "32%" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = Resolution,fill =(CV <= 20))) + scale_fill_brewer(palette = "Set2")+ 
  geom_bar(stat = "count",position = position_stack())+ 
  facet_grid(. ~ Ratio)+ theme_prism(base_size = 10)


# SC[Resolution == "60K" & is.na(CV) == FALSE]%>%
#   ggplot(aes(x = NCE,fill =(CV <= 20))) + scale_fill_brewer(palette = "Set2")+ 
#   geom_bar(stat = "count",position = position_stack())+ 
#   geom_text(stat='count', aes(label=..count..), position = position_stack(),
#             vjust=0.2, size = 3)+ 
#   facet_grid(. ~ Ratio)

dev.off()

SC[NCE == "32%" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = Resolution,fill =(CV <= 20))) + scale_fill_brewer(palette = "Set2")+ 
  geom_bar(stat = "count",position = position_stack())+ 
  facet_grid(. ~ Ratio)+ theme_prism()

SC[NCE == "32%"]%>%
  ggplot(aes(x = Resolution,fill =(CV <= 20))) + 
  geom_bar(stat = "count",position = position_stack())+ 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.2, size = 3)+ 
  facet_grid(. ~ Ratio)


ggplot(SC[Av14SN >= 1 & Resolution == "60K"], 
       aes(x = log2(Av14SN), y = log2(CV), color = log10(RawOvFtT))) + 
  geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
  facet_grid(NCE~Ratio,scales = "free") + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed")





