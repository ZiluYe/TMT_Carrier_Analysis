
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

All[,Av14 := sum(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,
                     SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                   na.rm = TRUE)/14, by = ID]

All <- merge(All, ScanInfo, by = "Identifier")


A2 <- All[`Parent intensity fraction` >= 0.98 & Av14 > 1 & Method %like% "SN_Raw"]

A2[Organism != "Ecoli",Y1 := log2(SN131C/SN132N)]%>%
  .[Organism != "Ecoli",Y2 := log2(SN129C/SN131N)]%>%
  .[Organism != "Ecoli",Y3 := log2(SN128C/SN130C)]%>%
  .[Organism != "Ecoli",Y4 := log2(SN133N/SN130N)]%>%
  .[Organism == "Ecoli",Y1 := log2(SN128N/SN128C)]%>%
  .[Organism == "Ecoli",Y2 := log2(SN129N/SN130C)]%>%
  .[Organism == "Ecoli",Y3 := log2(SN129C/SN131C)]%>%
  .[Organism == "Ecoli",Y4 := log2(SN130N/SN132C)]

A2 <- A2[is.na(Y1+ Y2+ Y3+Y4) == FALSE & is.infinite(Y1+ Y2+ Y3+Y4) == FALSE]

A2[,pvalue := t.test(c(Y1,Y2,Y3,Y4))[["p.value"]], by = ID]

A2[,FC := median(c(Y1,Y2,Y3,Y4),na.rm = TRUE), by = ID]

tA2 <- A2[Booster != "H" & AGC == "AGC300%"]

tA2[Organism == "Human",Index := (abs(FC) >= 0.5 & abs(FC) <= 1.5 & pvalue <= 0.05)]%>%
  .[Organism == "Yeast", Index := (FC >= 0.5 & FC <= 1.5 & pvalue <= 0.05)]%>%
  .[Organism == "Ecoli", Index := (FC <= -0.5 & FC >= -1.5 & pvalue <= 0.05)]

# ggplot(tA2,aes(x = Organism, fill = Index)) + geom_bar(stat = "count",position = "fill")+ 
#   facet_grid(Ratio ~ .)+ 
#   geom_text(stat='count', aes(label=..count..), position = position_stack(),
#             vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set2")


tA2 <- table(tA2[,.(Ratio,Index,Organism)])%>%data.table()%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]

fwrite(tA2[Organism != "Human"], "tables/Fig2a.txt")

pdf("figures/HEY_Volcano_Numbers_SN.pdf",width = 6,height = 6)

ggplot(tA2,aes(x = Organism, y = N, fill = Index)) + geom_bar(stat = "identity")+
  facet_wrap(.~Ratio)+
  geom_text(stat='identity', aes(label=N), position = position_stack(),
            vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set2") + theme_prism()+ 
  theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=8))


ggplot(tA2[Organism != "Human"],aes(x = Ratio, y = N, fill = Index)) + geom_bar(stat = "identity")+
  facet_grid(Organism~.)+
  geom_text(stat='identity', aes(label=N), position = position_stack(),
            vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set3") + theme_prism()


dev.off()

tA3 <- dcast(tA2, Ratio + Organism ~ Index, value.var = "N")

tA3[,Freq := `TRUE` / (`TRUE` + `FALSE`)]%>%
  .[,Ratio := factor(Ratio, levels = c("No126","14","42","98","210","434"))]
# 
# ggplot(tA3,aes(x = Organism, y = Freq)) +
#   facet_grid(Ratio ~ .)+
#   geom_text(stat='identity', aes(label=Freq, y = 0.75), size = 4)+ scale_fill_brewer(palette = "Set2")

tA3[,Freq := round(Freq * 100, digits = 1)]%>%
  .[,Freq := str_c(Freq,"%")]%>%
  .[,N := str_c("N = ", `TRUE`)]

tA3[Organism == "Ecoli",x := -1]%>%
  .[Organism == "Human", x := 0]%>%
  .[Organism == "Yeast", x := 1]

fwrite(A2[Booster != "H" & AGC == "AGC300%"],"tables/Fig2b.txt")

pdf("figures/HEY_Volcano_SN.pdf",width = 7,height = 12)

A2[Booster != "H" & AGC == "AGC300%"]%>%
  ggplot(aes(x = FC, y = -log10(pvalue), color = Organism)) + 
  geom_point(alpha = 0.3,aes(size = log10(Av14))) +
  geom_vline(xintercept = c(-1.5,-0.5, 0.5,1.5), col = "darkgrey", size = 1, linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "darkgrey", size = 1, linetype = "dashed")+ 
  facet_grid(Ratio ~ .) + xlim(-2,2) + 
  scale_size(range = c(0.01,2)) + theme_prism(base_size = 20)+ 
  geom_text(data = tA3,aes(x = x, y = 6,label = Freq),color = "black",size = 5)+ 
  geom_text(data = tA3,aes(x = x, y = 4.5,label = N),color = "black",size = 5)

dev.off()


A2[Organism == "Yeast", Index := (FC >= 0.5 & FC <= 1.5 & pvalue <= 0.05)]

fwrite(A2[Booster != "H" & AGC == "AGC300%" & Organism == "Yeast"], "tables/Fig2c.txt")

pdf("figures/Yeast_Regulated.pdf",width = 12,height = 5)

A2[Booster != "H" & AGC == "AGC300%" & Organism == "Yeast"]%>%
  ggplot(aes(x = FC, y = log10(SN133N), color = Index)) + 
  geom_point() +
  geom_vline(xintercept = c(1), col = "darkgrey", size = 1, linetype = "dashed") +  
  facet_grid(. ~ Ratio)  + theme_prism()

dev.off()

pdf("figures/Yeast_Regulated_Boxplot.pdf",width = 12,height = 2.5)

A2[Booster != "H" & AGC == "AGC300%" & Organism == "Yeast"]%>%
  ggplot(aes(x = Index, y = log10(SN133N), color = Index)) + 
  geom_boxplot()  +  
  facet_grid(. ~ Ratio)  + theme_prism()

dev.off()


grob <- grobTree(textGrob("Scatter plot", x=0.1,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=13, fontface="italic")))

# SCH: SN_Corrected Human

Yeast <- All[Organism == "Yeast" & `Parent intensity fraction` >= 0.98 & Av14 > 0 & Method %like% "SN"]%>%
  .[,SN_Index := `SN133C`/`SN134N`]%>%
  .[order(Av14,decreasing = TRUE)]


#Yeast <- Yeast[SN_Index <= 0.2]

Yeast[,rank := rank(Av14, ties.method="first")]

Index <- str_detect(colnames(Yeast), "1")

Y0 <- cbind(Yeast[,.(Identifier,Method)], Yeast[,..Index])

Y1 <- melt(Y0,id.vars = c("Identifier","Method"))

Y1 <- merge(Y1,Yeast[,.(Identifier,rank,Method,Av14,SN_Index)],by=c("Identifier","Method"))

Y1[,Channel := str_remove(variable,"SN")]%>%
  .[,Channel := factor(Channel, 
                       levels = c("126","127N","127C","128N","128C","129N","129C","130N",
                                  "130C","131N","131C","132N","132C","133N","133C","134N"))]%>%
  .[,Type := str_split_fixed(Method,"_",2)[,2]]


Y1 <- Y1[Channel != "126" & Channel != "127C"]

Y1 <- Y1[,Max := max(value,na.rm = TRUE), by = .(Type,rank)]%>%
  .[,value := value/Max * 100]

Y1[,Raw := str_split_fixed(Identifier,"__",2)[,1]]%>%
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
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]

# ggplot(Y1[Type == "Corrected" & Booster == "EY"], aes(x = Channel, y = value, color = Channel)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(size=0.1,width = 0.2,alpha = 0.3) + theme(legend.position = "none")+
#   labs(x = 'Channels', y = 'Relative intensity')+facet_grid(. ~ Ratio)

ggplot(Y1[Type == "Corrected" & Booster == "EY"], aes(x = Channel, y = value, color = Channel)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = 'Channels', y = 'Relative intensity')+facet_grid(. ~ Ratio)

ggplot(Y1[Type == "Corrected" & Booster == "EY" & Channel == "127N"], 
       aes(x = log10(Av14) >= 2, y = value,color = log10(Av14) >= 2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.1,width = 0.2,alpha = 0.3) +
  labs(x = 'Channels', y = 'Relative intensity')+facet_grid(. ~ Ratio)


ggplot(Y1[Type == "Corrected" & Booster == "EY" & Channel == "130N"], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1,alpha = 0.3)+ scale_color_gradient(low = "white",high = "red") +
  labs(x = 'log10(Av14)', y = 'Relative intensity of 130N')+facet_grid(. ~ Ratio)

Y2 <- data.table(value = c(50,90,80,70,60,50,40,30,20,10,0,100,1,99),
                 Channel = c("127N","128N","128C","129N","129C","130N",
                             "130C","131N","131C","132N","132C","133N","133C","134N"))

Y2[,Channel := factor(Channel,levels = c("127N","128N","128C","129N","129C","130N",
                                         "130C","131N","131C","132N","132C","133N","133C","134N"))]

fwrite(Y1[Type == "Corrected" & (Booster == "EY" | Booster == "No") & Av14 >= 1], "tables/FigS4.txt")

pdf("figures/Yeast_Quant.pdf",width = 20,height = 12)


ggplot(Y1[Type == "Corrected" & (Booster == "EY" | Booster == "No") & Av14 >= 1], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1)+ scale_color_gradient(low = "white",high = "red") +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed") + theme_prism(base_size = 20)

ggplot(Y1[Type == "Corrected" & (Booster == "EY" | Booster == "No") & Av14 >= 1], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1,alpha = 0.3)+ scale_color_gradient(low = "white",high = "red") +
  geom_violin(fill = NA) +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed") + theme_prism(base_size = 20)


ggplot(Y1[Type == "Raw" & (Booster == "EY" | Booster == "No") & Av14 >= 1], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1)+ scale_color_gradient(low = "white",high = "red") +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed") + theme_prism(base_size = 20)

ggplot(Y1[Type == "Raw" & (Booster == "EY" | Booster == "No") & Av14 >= 1], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1,alpha = 0.3)+ scale_color_gradient(low = "white",high = "red") +
  geom_violin(fill = NA) +
  labs(x = 'log10(Av14)', y = 'Relative intensity_Raw')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed") + theme_prism(base_size = 20)

dev.off()

ggplot(Y1[rank == 13400],
       aes(x = Channel,y=value))+geom_bar(stat = "Identity",fill = "skyblue")+
  facet_grid(Type ~.,scales = "free_y") 


Y3 <- Y1[Method == "SN_Corrected"]

Y3[is.na(value) == TRUE, value := 0]

Y3 <- merge(Y3,Y2,by = "Channel")

Y3[,value := value.x + 100 - value.y]

Y3 <- dcast(Y3,Identifier + Av14 ~ variable, value.var = "value")

Y3[,ID := 1:dim(Y3)[1]]

Y3[,mean := mean(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                   SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]

Y3[,Raw := str_split_fixed(Identifier,"__",2)[,1]]%>%
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
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]


ggplot(Y3[(Booster == "EY"| Booster == "No") & is.na(CV) == FALSE], 
       aes(x = Ratio, fill = CV <= 10)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~.,scales = "free")







gplot(Yeast[(Booster == "EY"| Booster == "No")], 
       aes(x = Ratio, y = log2(`SN132N`/`SN133C`),fill = Ratio)) + geom_boxplot()+
  scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~.,"free_y")















Ecoli <- All[Organism == "Ecoli" & `Parent intensity fraction` >= 0.9 & Av14 > 1 & Method %like% "SN"]%>%
  .[,SN_Index := `SN133C`/`SN134N`]%>%
  .[order(Av14,decreasing = TRUE)]

Ecoli <- Ecoli[SN_Index >= 5]

Ecoli[,rank := rank(Av14, ties.method="first")]

Index <- str_detect(colnames(Ecoli), "1")

Y0 <- cbind(Ecoli[,.(Identifier,Method)], Ecoli[,..Index])

Y1 <- melt(Y0,id.vars = c("Identifier","Method"))

Y1 <- merge(Y1,Ecoli[,.(Identifier,rank,Av14,SN_Index)],by="Identifier", allow.cartesian=TRUE)

Y1[,Channel := str_remove(variable,"SN")]%>%
  .[,Channel := factor(Channel, 
                       levels = c("126","127N","127C","128N","128C","129N","129C","130N",
                                  "130C","131N","131C","132N","132C","133N","133C","134N"))]%>%
  .[,Type := str_split_fixed(Method,"_",2)[,2]]


Y1 <- Y1[Channel != "126" & Channel != "127C"]

Y1 <- Y1[,Max := max(value,na.rm = TRUE), by = .(Type,rank)]%>%
  .[,value := value/Max * 100]

Y1[,Raw := str_split_fixed(Identifier,"__",2)[,1]]%>%
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
  .[,Rep := str_split_fixed(`Raw`, "%_",2)[,2]]

Y2 <- data.table(value = c(50,10,20,30,40,50,60,70,80,90,100,0,99,1),
                 Channel = c("127N","128N","128C","129N","129C","130N",
                             "130C","131N","131C","132N","132C","133N","133C","134N"))

Y2[,Channel := factor(Channel,levels = c("127N","128N","128C","129N","129C","130N",
                                         "130C","131N","131C","132N","132C","133N","133C","134N"))]

pdf("figures/Ecoli_Quant.pdf",width = 10,height = 8)


ggplot(Y1[Type == "Corrected" & (Booster == "EY" | Booster == "No")], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1)+ scale_color_gradient(low = "white",high = "red") +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed")

ggplot(Y1[Type == "Corrected" & (Booster == "EY" | Booster == "No")], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1,alpha = 0.3)+ scale_color_gradient(low = "white",high = "red") +
  geom_violin(fill = NA) +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed")


ggplot(Y1[Type == "Raw" & (Booster == "EY" | Booster == "No")], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1)+ scale_color_gradient(low = "white",high = "red") +
  labs(x = 'log10(Av14)', y = 'Relative intensity')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed")

ggplot(Y1[Type == "Raw" & (Booster == "EY" | Booster == "No")], 
       aes(x = log10(Av14), y = value,color = log10(Av14))) +
  geom_point(size=0.1,alpha = 0.3)+ scale_color_gradient(low = "white",high = "red") +
  geom_violin(fill = NA) +
  labs(x = 'log10(Av14)', y = 'Relative intensity_Raw')+facet_grid(Ratio ~ Channel) +
  geom_hline(data = Y2,aes(yintercept = value), size = 0.5, linetype = "dashed")

dev.off()

ggplot(Y1[rank == 13400],
       aes(x = Channel,y=value))+geom_bar(stat = "Identity",fill = "skyblue")+
  facet_grid(Type ~.,scales = "free_y") 