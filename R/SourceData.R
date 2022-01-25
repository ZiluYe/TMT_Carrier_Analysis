
library(data.table)
library(magrittr)
library(writexl)

Fig1b <- fread("tables/Fig1b.txt")
Fig1b <- Fig1b[,.(Accession,Booster,Ratio,Organism)]%>%unique()
colnames(Fig1b)[2] <- "Carrier"

Fig1c <- fread("tables/Fig1c.txt")
Fig1c <- Fig1c[,.(Accession,Ratio,CV)]%>%unique()

Fig1d1 <- fread("tables/Fig1de_Fig.S2a.txt")
Fig1d1 <- Fig1d1[,.(Identifier,Av14SN,CV)]
Fig1d1[,Channels := "Av14"]
colnames(Fig1d1)[2] <- "Average SN"


Fig1d2 <- fread("tables/Fig1de_Fig.S3a.txt")
Fig1d2 <- Fig1d2[,.(Identifier,Av12SN,CV)]
Fig1d2[,Channels := "Av12"]
colnames(Fig1d2)[2] <- "Average SN"

Fig1d_1e_S2a_S2b <- rbind(Fig1d1,Fig1d2)

Fig1f <- fread("tables/Fig1f.txt")
Fig1f <- Fig1f[,.(Identifier,Comment,variable,value)]%>%unique()
Fig1f <- Fig1f[Identifier %like% "13"]%>%unique()

Fig2a <- fread("tables/Fig2a.txt")

Fig2b <- fread("tables/Fig2b.txt")
Fig2b <- Fig2b[,.(Identifier,pvalue,FC,Organism)]%>%unique()

Fig2c <- fread("tables/Fig2c.txt")
Fig2c <- Fig2c[,.(Identifier,pvalue,FC,Organism,Index,SN133N)]%>%unique()

Fig3a <- fread("tables/Fig3a.txt")
Fig3a <- Fig3a[,.(Identifier,mean,sd,CV,NCE)]%>%unique()

Fig3b <- fread("tables/Fig3b.txt")

Fig3c <- fread("tables/Fig3c.txt")
Fig3c <- Fig3c[,.(Identifier,mean,sd,CV,Resolution)]%>%unique()

Fig3d <- fread("tables/Fig3d.txt")

Fig4 <- fread("tables/Fig4.txt")
Fig4 <- Fig4[,.(Sample,rn,variable,value)]%>%unique()
colnames(Fig4) <- c("Sample","MS1","MS2","Correlation")

FigS1 <- fread("tables/Fig.S1a.txt")
FigS1 <- FigS1[,.(Accession,Raw,Organism)]%>%unique()

FigS2bc <- fread("tables/Fig.S2bc")
FigS2bc <- FigS2bc[,.(Raw,Accession,Av14SN,Organism,CV)]%>%unique()

FigS3bc <- fread("tables/Fig.S3bc")
FigS3bc <- FigS3bc[,.(Raw,Accession,Av12SN,Organism,CV)]%>%unique()

FigS4 <- fread("tables/FigS4.txt")
FigS4 <- FigS4[,.(Identifier,variable,value,Av14)]%>%unique()

FigS5a <- fread("tables/Fig.S5a.txt")
FigS5a <- FigS5a[,.(Raw,Precursor,NCE,XCorr)]%>%unique()

FigS5b <- fread("tables/Fig.S5b.txt")
FigS5b <- FigS5b[,.(Raw,`Scan number`,NCE,Score)]%>%unique()

FigS6 <- fread("tables/Fig.S6.txt")
FigS6 <- FigS6[,.(Sample,rn,variable,value)]%>%unique()
colnames(FigS6) <- c("Sample","RI1","RI2","Correlation")

write_xlsx(list(Fig1b = Fig1b,
                Fig1c = Fig1c,
                Fig1d_1e_S2a_S2b  = Fig1d_1e_S2a_S2b,
                Fig1f  = Fig1f,
                Fig2a  = Fig2a,
                Fig2b  = Fig2b,
                Fig2c  = Fig2c,
                Fig3a  = Fig3a,
                Fig3b  = Fig3b,
                Fig3c  = Fig3c,
                Fig3d  = Fig3d,
                Fig4  = Fig4,
                FigS1  = FigS1,
                FigS2bc  = FigS2bc,
                FigS3bc  = FigS3bc,
                FigS4  = FigS4,
                FigS5a  = FigS5a,
                FigS5b  = FigS5b,
                FigS6  = FigS6),"tables/SourceData.xlsx")

