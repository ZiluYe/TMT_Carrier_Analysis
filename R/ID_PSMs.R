
source("R/Functions.R")

# SC: SN_Corrected

SC <- fread("tables/PSMs_Annotated.txt")[Method == "SN_Corrected" & !(Sample %like% "K_|0ms")]

SC <- SC[Organism != "" & Organism != "Nonunique"]

SC[,Identifier := str_remove(`Spectrum File`,"\\.raw")]%>%
  .[,Identifier := str_c(`Identifier`,`First Scan`,sep = "__")]

SC <- SC[Confidence == "High",.(Identifier,Organism)]%>%unique()

SC[,Raw := str_split_fixed(Identifier,"__",2)[,1]]%>%
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


ggplot(SC[Rep == 1 & Amount == "200pg"], aes(x = Ratio, fill = AGC)) + 
  geom_bar(position = position_dodge(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(Organism ~ Booster, scales = "free")

ggplot(SC[Booster == "HEY"| Booster %like% "No"], aes(x = Ratio, fill = Rep)) + 
  geom_bar(position = position_stack(),alpha = 0.9) + scale_fill_brewer(palette = "Set2") + 
  facet_grid(.~Amount, scales = "free")
