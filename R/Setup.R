
source("R/Functions.R")


HEY <- data.table(Channel = as.factor(c("127N","128N","128C","129N","129C","130N","130C",
                                        "131N","131C","132N","132C","133N","133C","134N")),
                  Human = c(100,100,100,100,100,100,100,100,100,100,100,100,100,100),
                  Yeast = c(50,90,80,70,60,50,40,30,20,10,0,100,1,99),
                  Ecoli = c(50,10,20,30,40,50,60,70,80,90,100,0,99,1))

HEY <- melt(HEY,id.vars = "Channel")

HEY[,Channel := factor(Channel,levels = c("127N","128N","128C","129N","129C","130N","130C",
                                          "131N","131C","132N","132C","133N","133C","134N"))]%>%
  .[,Channel := fct_rev(Channel)]%>%
  .[,variable := factor(variable, levels = c("Yeast","Human","Ecoli"))]


pdf("figures/HEY_Setup_v2.pdf",width = 4,height = 8)

ggplot(HEY,aes(y = Channel, x = value,fill = variable)) + 
  geom_bar(stat = "identity",position = position_stack()) + 
  geom_text(stat='identity', aes(label=value), position = position_stack(),hjust = 1, size = 4) + 
  theme_prism()

 dev.off()
