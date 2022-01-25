
source("R/Functions.R")

# SC: SN_Corrected

ScanInfo <- fread("data/ScanInfo_MQ.txt")[!(Identifier %like% "K_|0ms")]

All <- fread("tables/PSMs_Annotated.txt")[!(Sample %like% "K_|0ms")]

All <- All[Organism != "" & Organism != "Nonunique"]

All[,Identifier := str_remove(`Spectrum File`,"\\.raw")]%>%
  .[,Identifier := str_c(`Identifier`,`First Scan`,sep = "__")]%>%
  .[,ID := 1:dim(All)[1]]

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

All[,Av14SN := sum(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,
                    SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                  na.rm = TRUE)/14, by = ID]%>%
  .[, Sum16SN := sum(c(SN126,SN127N,SN127C,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,
                    SN131N,SN131C,SN132N,SN132C,SN133N,SN133C,SN134N),
                  na.rm = TRUE), by = ID]

All <- merge(All, ScanInfo, by = "Identifier", all.y = TRUE)

# SCH: SN_Corrected Human

SCH <- All[Organism == "Human"]


SCH[,mean := mean(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                    SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,sd := sd(c(SN127N,SN128N,SN128C,SN129N,SN129C,SN130N,SN130C,SN131N,SN131C,
                SN132N,SN132C,SN133N,SN133C,SN134N),na.rm = TRUE), by = ID]%>%
  .[,CV := sd/mean * 100]%>%
  .[,ID := 1:dim(SCH)[1]]

SCH[Booster == "HEY" & Method %like% "SN_Corrected" & is.na(CV) == FALSE]%>%
  ggplot(aes(x = AGC,fill =(CV <= 20))) + 
  geom_bar(stat = "count",position = position_stack())+ 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.2, size = 4)+ 
  facet_grid(Amount ~ Ratio)

 

#pdf("figures/CV_PSMs.pdf",width = 10,height = 5)

ggplot(SCH[Av14SN >= 1 & Booster == "HEY" & Amount == "200pg" & is.na(CV) == FALSE], 
       aes(x = Method, fill = CV <= 20)) + geom_bar() + 
  geom_text(stat='count', aes(label=..count..), position = position_stack(),
            vjust=0.1, size = 4)+ scale_fill_brewer(palette = "Set2")+
  facet_grid(AGC + Confidence~Ratio,scales = "free")

ggplot(SCH[Av14SN >= 1 & Booster == "HEY" & Method == "SN_Corrected" & Amount == "200pg"], 
       aes(x = log2(Av14SN), y = log2(CV), color = log10(RawOvFtT))) + 
  geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
  facet_grid(Confidence~Ratio) + 
  geom_hline(yintercept = log2(20), size = 1, linetype = "dashed")

dev.off()



ggplot(SCH[Av14SN >= 1 & Booster == "HEY" & Method == "SN_Corrected" & Amount == "200pg"], 
       aes(x = XCorr, y = log2(RawOvFtT - Sum16SN ), color = log10(Sum16SN))) + 
  geom_point(size = 0.3)+scale_color_gradient(low = "white",high = "red")+
  facet_grid(.~Confidence) 
