library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
DeleInfo=read.table("Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_SolMsaNonMajorDele_Genotype_GERP2.75.txt.filter.rmBothMajor.txt",header=T)

both_non_dele<-DeleInfo_rmBothMajor[DeleInfo_rmBothMajor$A626==0 & DeleInfo_rmBothMajor$E463==0 & !is.na(DeleInfo_rmBothMajor$E463) & !is.na(DeleInfo_rmBothMajor$A626) & DeleInfo_rmBothMajor$Lan_DeleFreq >0,]
both_dele<-DeleInfo_rmBothMajor[DeleInfo_rmBothMajor$A626==2 & DeleInfo_rmBothMajor$E463==2 & !is.na(DeleInfo_rmBothMajor$A626) & !is.na(DeleInfo_rmBothMajor$E463),]
A626_dele<-DeleInfo_rmBothMajor[DeleInfo_rmBothMajor$A626==2 & !is.na(DeleInfo_rmBothMajor$A626),]
E463_dele<-DeleInfo_rmBothMajor[DeleInfo_rmBothMajor$E463==2 & !is.na(DeleInfo_rmBothMajor$E463),]

nrow(both_dele)
nrow(both_non_dele)
nrow(A626_dele)
nrow(E463_dele)

summary(both_dele$Lan_DeleFreq)
summary(both_dele$Can_DeleFreq)
summary(both_non_dele$Lan_DeleFreq)
summary(both_non_dele$Can_DeleFreq)
summary(A626_dele$Lan_DeleFreq)
summary(A626_dele$Can_DeleFreq)
summary(E463_dele$Lan_DeleFreq)
summary(E463_dele$Can_DeleFreq)

both_dele_land <- both_dele
both_dele_land$type <- "F1 homo dele"
both_dele_land$population <- "Lanrace"
both_dele_land$DeleFreq <- both_dele_land$Lan_DeleFreq

both_dele_scand <- both_dele
both_dele_scand$type <- "F1 homo dele"
both_dele_scand$population <- "Scand"
both_dele_scand$DeleFreq <- both_dele_scand$Can_DeleFreq

both_non_dele_land <- both_non_dele
both_non_dele_land$type <- "Both neutral"
both_non_dele_land$population <- "Lanrace"
both_non_dele_land$DeleFreq <- both_non_dele_land$Lan_DeleFreq

both_non_dele_scand <- both_non_dele
both_non_dele_scand$type <- "Both neutral"
both_non_dele_scand$population <- "Scand"
both_non_dele_scand$DeleFreq <- both_non_dele_scand$Can_DeleFreq

A626_dele_land <- A626_dele
A626_dele_land$type <- "A6-26 homo dele"
A626_dele_land$population <- "Lanrace"
A626_dele_land$DeleFreq <- A626_dele_land$Lan_DeleFreq

A626_dele_scand <- A626_dele
A626_dele_scand$type <- "A6-26 homo dele"
A626_dele_scand$population <- "Scand"
A626_dele_scand$DeleFreq <- A626_dele_scand$Can_DeleFreq

E463_dele_land <- E463_dele
E463_dele_land$type <- "E4-63 homo dele"
E463_dele_land$population <- "Lanrace"
E463_dele_land$DeleFreq <- E463_dele_land$Lan_DeleFreq

E463_dele_scand <- E463_dele
E463_dele_scand$type <- "E4-63 homo dele"
E463_dele_scand$population <- "Scand"
E463_dele_scand$DeleFreq <- E463_dele_scand$Can_DeleFreq

input_list <- list(
  both_dele_land[, c("type", "population", "DeleFreq")],
  both_dele_scand[, c("type", "population", "DeleFreq")]
)

input_list[[3]] <- both_non_dele_land[, c("type", "population", "DeleFreq")]
if(exists("both_non_dele_scand")) {
  input_list[[4]] <- both_non_dele_scand[, c("type", "population", "DeleFreq")]
}

input_list[[length(input_list)+1]] <- A626_dele_land[, c("type", "population", "DeleFreq")]
if(exists("A626_dele_scand")) {
  input_list[[length(input_list)+1]] <- A626_dele_scand[, c("type", "population", "DeleFreq")]
}

input_list[[length(input_list)+1]] <- E463_dele_land[, c("type", "population", "DeleFreq")]
if(exists("E463_dele_scand")) {
  input_list[[length(input_list)+1]] <- E463_dele_scand[, c("type", "population", "DeleFreq")]
}

input <- do.call(rbind, input_list)

input$type <- factor(input$type, 
                     levels = c("F1 homo dele", "A6-26 homo dele", "E4-63 homo dele", "Both neutral"))


pdf(file="Fig5.pdf", width=8, height=6)
ggplot(input, aes(x = type, y = DeleFreq, fill = population)) +
  geom_boxplot(outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Lanrace" = "palegreen4", 
                               "Scand" = "orange2"),
                    name = "Group") +
  theme_minimal() +
  labs(x = "", y = "Deleterious allele frequency") +
  theme(panel.border = element_rect(color = 'black', linewidth = 1.5),
    panel.grid.major = element_line(linewidth = 1.2),
    panel.grid.minor = element_line(linewidth = 0.6),
    axis.title = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 13, color = "black",angle=30,hjust=1,vjust=1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = c(0.65,0.85),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.background = element_blank())
dev.off()




























































































































































































































