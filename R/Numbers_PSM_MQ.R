
source("R/Functions.R")

PSM <- fread("../../../CPR/TMT_BoosterEffect/TMTdata/data/combined/txt/msmsScans.txt")

PSM <- PSM[(`Raw file` %like% "K_")]

PSM[,Identifier := str_c(`Raw file`,`Scan number`,sep = "__")]

PSM[,Human := str_detect(`Proteins`, "HUMAN")]%>%
  .[,Yeast := str_detect(`Proteins`, "YEAST")]%>%
  .[,Ecoli := str_detect(`Proteins`, "ECOLI")]%>%
  .[,Organism := "Nonunique"]%>%
  .[Human == TRUE & Yeast == FALSE & Ecoli == FALSE,Organism := "Human"]%>%
  .[Human == FALSE & Yeast == TRUE & Ecoli == FALSE,Organism := "Yeast"]%>%
  .[Human == FALSE & Yeast == FALSE & Ecoli == TRUE,Organism := "Ecoli"]%>%
  .[,Organism := factor(Organism, levels = c("Human", "Ecoli", "Yeast", "Nonunique"))]

SC <- PSM[Identified == "+" & Organism == "Human"]%>%unique()

SC[,Raw := str_split_fixed(Identifier, "__",2)[,1]]%>%
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

fwrite(SC[Resolution == "60K"],"tables/Fig.S5b.txt")

pdf("figures/Andromeda_Box.pdf",width = 10,height = 5)

SC[Resolution == "60K"]%>%
  ggplot(aes(x=NCE, y= log10(Score),color = NCE)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10)


SC[NCE == "32%"]%>%
  ggplot(aes(x=Resolution, y= log10(Score),color = Resolution)) +
  geom_boxplot(outlier.size = 1) + facet_grid(. ~ Ratio)+
  stat_summary(fun=median, geom="point", shape=20, size=3, color="red", fill="red") + 
  theme_prism(base_size = 10)

dev.off()

