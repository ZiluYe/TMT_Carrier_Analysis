
source("R/Functions.R")

Files <- fread("data/20200820_EXPL8_EVO1_ZY_SA_44min_TMT_270xR_SN_QuanSpectra.txt")

Files <- Files[!(`Spectrum File` %like% "150ms|500ms"),.(`Spectrum File`)]%>%unique()

Files[,Raw := str_remove(`Spectrum File`,".raw")]%>%
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
  .[,Resolution := str_split_fixed(`Raw`, "100pg_",2)[,2]]%>%
  .[,Resolution := str_split_fixed(Resolution, "_",2)[,1]]%>%
  .[,Resolution := factor(Resolution, levels = c("30K","45K","60K","120K"))]%>%
  .[,NCE := str_split_fixed(`Raw`, "K_",2)[,2]]%>%
  .[,NCE := str_split_fixed(NCE, "_",2)[,1]]%>%
  .[is.na(AGC), AGC := "AGC300%"]%>%
  .[is.na(Resolution), Resolution := "60K"]%>%
  .[NCE == "", NCE := "32%"]


writexl::write_xlsx(Files,"tables/RawFiles.xlsx")


PSMs <- fread("tables/PSMs_Annotated.txt")[!(Sample %like% "150ms|500ms") & Confidence == "High"]

PSMs[,`:=` (ID = NULL, `PSMs Workflow ID` = NULL, `PSMs Peptide ID` = NULL, Checked = NULL, Confidence = NULL)]

Proteins <- fread("tables/Proteins.txt")

Proteins <- Proteins[Master == "IsMasterProtein"]%>%
  .[,Human := str_detect(`Description`, "Homo sapiens")]%>%
  .[,Yeast := str_detect(`Description`, "Saccharomyces cerevisiae")]%>%
  .[,Ecoli := str_detect(`Description`, "Escherichia coli")]%>%
  .[,Organism := "Nonunique"]%>%
  .[Human == TRUE & Yeast == FALSE & Ecoli == FALSE,Organism := "Human"]%>%
  .[Human == FALSE & Yeast == TRUE & Ecoli == FALSE,Organism := "Yeast"]%>%
  .[Human == FALSE & Yeast == FALSE & Ecoli == TRUE,Organism := "Ecoli"]

Proteins <- Proteins[!(Sample %like% "150ms|500ms")]

writexl::write_xlsx(list(RawFiles = Files,
                         Proteins = Proteins,
                         PSMs_SNR_Corrected = PSMs[Method == "SN_Corrected"],
                         PSMs_SNR_Raw = PSMs[Method == "SN_Raw"],
                         PSMs_Intensity_Corrected = PSMs[Method == "Int_Corrected"],
                         PSMs_Intensity_Raw = PSMs[Method == "Int_Raw"]),"tables/SupplementaryTable1.xlsx")

