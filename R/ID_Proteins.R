
source("R/Functions.R")

# SC: SN_Corrected

SC <- fread("tables/Proteins.txt")[Method == "SN_Corrected" & !(Sample %like% "K_|0ms")]

SC <- SC[Master == "IsMasterProtein"]


SC[,Human := str_detect(`Description`, "Homo sapiens")]%>%
  .[,Yeast := str_detect(`Description`, "Saccharomyces cerevisiae")]%>%
  .[,Ecoli := str_detect(`Description`, "Escherichia coli")]%>%
  .[,Organism := "Nonunique"]%>%
  .[Human == TRUE & Yeast == FALSE & Ecoli == FALSE,Organism := "Human"]%>%
  .[Human == FALSE & Yeast == TRUE & Ecoli == FALSE,Organism := "Yeast"]%>%
  .[Human == FALSE & Yeast == FALSE & Ecoli == TRUE,Organism := "Ecoli"]


SC <- SC[Organism != "" & Organism != "Nonunique"]

SC[,Raw := Sample]%>%
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


# ggplot(SC[Rep == 1 & Amount == "200pg" & AGC == "AGC300%"], aes(x = Ratio, fill = Organism)) + 
#   geom_bar(position = position_dodge(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
#   facet_grid(Organism ~ Booster, scales = "free")

p1 <- SC[Rep == 2 & Amount == "50pg" & AGC == "AGC300%" & Booster != "No"]

p2 <- SC[Rep == 2 & Amount == "50pg" & AGC == "AGC300%" & Booster == "No"]

p2 <- table(p2$Organism)%>%data.table()

colnames(p2)[1] <- "Organism"


p3 <- copy(p2)

p3 <- p3[Organism != "Human", Organism := "EY"]%>%
  .[,N := sum(N), by = Organism]%>%unique()


fwrite(p1,"tables/Fig1b.txt")

pdf("figures/Identified_Proteins&Organism.pdf",width = 10,height = 5)

ggplot(p1, aes(x = Booster, fill = Ratio)) + 
  geom_bar(position = position_dodge(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(.~Organism, scales = "free")+
  geom_hline(data = p2,aes(yintercept = N), size = 0.5, linetype = "dashed") + theme_prism()

ggplot(p1, aes(x = Booster, fill = Ratio)) + 
  geom_bar(position = position_dodge(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(.~ Organism == "Human", scales = "free")+
  geom_hline(data = p3,aes(yintercept = N), size = 0.5, linetype = "dashed") + theme_prism()


dev.off()

fwrite(SC[Booster == "H" | Booster %like% "No" & Organism == "Human"], "tables/Fig.S1a.txt")

pdf("figures/Identified_Proteins.pdf",width = 10,height = 5)

ggplot(SC[Rep == 1 & Amount == "200pg"], aes(x = Ratio, fill = AGC)) + 
  geom_bar(position = position_dodge(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(Organism ~ Booster, scales = "free")+ theme_prism()

ggplot(SC[Booster == "HEY" | Booster %like% "No" & Organism == "Human"], aes(x = Ratio, fill = Rep)) + 
  geom_bar(position = position_stack(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(AGC~Amount, scales = "free") + theme_prism()

ggplot(SC[Booster == "H" | Booster %like% "No" & Organism == "Human"], aes(x = Ratio, fill = Rep)) + 
  geom_bar(position = position_stack(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(AGC~Amount, scales = "free") + theme_prism()

dev.off()


SC2 <- SC[,.(Accession,Organism,Raw)]%>%unique()

SC2[,Booster := str_split_fixed(`Raw`, "TMT_",2)[,2]]%>%
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


SC2 <- SC2[(Booster == "H" | Booster %like% "No") & Organism == "Human" & AGC == "AGC300%" & Amount == "100pg"]

fwrite(SC2, "tables/Fig.S1b.txt")

lt <- list(NC_1 = SC2[Ratio == "No126" & Rep == "1"]$Accession%>%unique(),
           NC_2 = SC2[Ratio == "No126" & Rep == "2"]$Accession%>%unique(),
           `14_1` = SC2[Ratio == "14" & Rep == "1"]$Accession%>%unique(),
           `14_2` = SC2[Ratio == "14" & Rep == "2"]$Accession%>%unique(),
           `42_1` = SC2[Ratio == "42" & Rep == "1"]$Accession%>%unique(),
           `42_2` = SC2[Ratio == "42" & Rep == "2"]$Accession%>%unique(),
           `98_1` = SC2[Ratio == "98" & Rep == "1"]$Accession%>%unique(),
           `98_2` = SC2[Ratio == "98" & Rep == "2"]$Accession%>%unique(),
           `210_1` = SC2[Ratio == "210" & Rep == "1"]$Accession%>%unique(),
           `210_2` = SC2[Ratio == "210" & Rep == "2"]$Accession%>%unique(),
           `434_1` = SC2[Ratio == "434" & Rep == "1"]$Accession%>%unique(),
           `434_2` = SC2[Ratio == "434" & Rep == "2"]$Accession%>%unique())


library(ComplexHeatmap)
library(RColorBrewer)

m = make_comb_mat(lt)

m = m[comb_size(m) >= 10]

pdf("figures/UpSet_Proteins.pdf",width = 12,height = 6)

UpSet(m, comb_col = brewer.pal(12, "Paired")[comb_degree(m)])

dev.off()
