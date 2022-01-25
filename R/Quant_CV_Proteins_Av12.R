
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

All[,Av12SN := sum(c(SN128N,SN129N,SN129C,SN130N,SN130C,
                    SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                  na.rm = TRUE)/12, by = ID]

# SCH: SN_Corrected Human

SCH <- All[Organism == "Human"]

SCH[,mean := mean(c(SN128N,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                    SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN128N,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]%>%
  .[,ID := 1:dim(SCH)[1]]

table(SCH$CV >= 20)


pdf("figures/Proteins_ID&CV_Av12.pdf",width = 5,height = 5)

SCH[Booster %like% "No|HEY" & AGC == "AGC300%" & Amount == "50pg" & Method == "SN_Corrected" & Rep == "2" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = Ratio, fill = CV <= 20)) + 
  geom_bar(position = position_stack(),alpha = 0.9) + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.2, size = 5)+ scale_fill_brewer(palette = "Set2") + theme_prism()


dev.off()

pdf("figures/CV_Proteins_Av12.pdf",width = 10,height = 5)

ggplot(SCH[Booster == "HEY" & Amount == "50pg" & is.na(CV) == FALSE  & Rep == "2"], 
       aes(x = Method, fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 3)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~Ratio,scales = "free") + theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=6))+ 
  theme_prism() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))

ggplot(SCH[Booster %like% "HEY|No" & Amount == "50pg" & is.na(CV) == FALSE  & Rep == "2" & Method %like% "SN"], 
       aes(x = str_split_fixed(Method,"_",2)[,2], fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 5)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~Ratio,scales = "free") + theme_prism()

ggplot(SCH[Booster %like% "HEY|No" & Amount == "50pg" & is.na(CV) == FALSE  & Rep == "2" & Method %like% "SN"], 
       aes(x = str_split_fixed(Method,"_",2)[,2], fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 5)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC~Ratio,scales = "free") + theme_prism() + theme(legend.position = "none")

ggplot(SCH[Av12SN >= 1 & Booster %like% "HEY|No" & Method == "SN_Corrected" & Amount == "50pg"  & Rep == "2"], 
       aes(x = log2(Av12SN), y = log2(CV))) + 
  geom_point(size = 0.3, color = "skyblue")+scale_color_gradient(low = "white",high = "red")+
  facet_grid(AGC~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed") + theme_prism()+ 
  stat_cor(method = "spearman",r.digits = 3, label.sep = "\n") + ggtitle("Proteins: SN_Corrected")

dev.off()

fwrite(SCH[Av12SN >= 1 & Booster %like% "HEY|No" & Method == "SN_Corrected" & Amount == "50pg"  & Rep == "2"],"tables/Fig.S3bc")
