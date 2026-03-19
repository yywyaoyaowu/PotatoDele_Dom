library(ggplot2)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
setwd("/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/08_PCA/eigenvec")

options(digits = 3)
eigenvel_pca338 <- read.table("pca338.eigenval")
eigenvel_pca233 <- read.table("Pre_pca233.eigenval")

PC5_explain_pca338 <- round((eigenvel_pca338$V1[c(1:5)] / sum(eigenvel_pca338$V1)) * 100, digits = 2)
PC5_explain_pca233 <- round((eigenvel_pca233$V1[c(1:5)] / sum(eigenvel_pca233$V1)) * 100, digits = 2)

PC5_explainCum <- data.frame("expanded" = cumsum(PC5_explain_pca338), "Previous" = cumsum(PC5_explain_pca233))
PC5_explainCum$Index <- 1:nrow(PC5_explainCum)
df_longPC5 <- melt(PC5_explainCum, id.vars = "Index", variable.name = "Variable", value.name = "Value")
print(head(df_longPC5))


pdf("explained_variancesPC5.pdf", width = 7, height = 7)
ggplot(df_longPC5, aes(x = factor(Index), y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  scale_fill_manual(values = c("expanded" = "#F0A419", "Previous" = "#F3C97F")) +
  labs(title = "Comparison of expanded vs previous PCA", 
       x = "Principal Components", 
       y = "Percentage of Cumulative explained variances(%)", 
       fill = "Variable") +
  theme_bw() +
  theme(
    panel.border = element_rect(color='black',linewidth=1.5),
    panel.grid.major = element_line(linewidth = 1.2),
    panel.grid.minor = element_line(linewidth = 0.6),
    axis.title = element_text(size=15,color="black"),
    axis.text = element_text(size=13,color="black"),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    legend.direction = "horizontal",
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text=element_text(size=13),
    plot.title = element_text(size=15,color="black", hjust = 0.5, vjust = 1),
    legend.background = element_rect(color = "black", linewidth = 0.5, linetype = "solid"))
dev.off()



PC3_explain_pca338 <- round((eigenvel_pca338$V1[c(1:3)] / sum(eigenvel_pca338$V1)) * 100, digits = 2)
PC3_explain_pca233 <- round((eigenvel_pca233$V1[c(1:3)] / sum(eigenvel_pca233$V1)) * 100, digits = 2)

PC3_explainCum <- data.frame("expanded" = cumsum(PC3_explain_pca338), "Previous" = cumsum(PC3_explain_pca233))
PC3_explainCum$Index <- 1:nrow(PC3_explainCum)
df_longPC3 <- melt(PC3_explainCum, id.vars = "Index", variable.name = "Variable", value.name = "Value")
print(head(df_longPC3))

pdf("explained_variancesPC3.pdf", width = 5, height = 7)
ggplot(df_longPC3, aes(x = factor(Index), y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  scale_fill_manual(values = c("expanded" = "#F0A419", "Previous" = "#F3C97F")) +
  labs(title = "Comparison of expanded vs previous PCA", 
       x = "Principal Components", 
       y = "Percentage of Cumulative explained variances(%)", 
       fill = "Variable") +
  theme_bw() +
  theme(
    panel.border = element_rect(color='black',linewidth=1.5),
    panel.grid.major = element_line(linewidth = 1.2),
    panel.grid.minor = element_line(linewidth = 0.6),
    axis.title = element_text(size=15,color="black"),
    axis.text = element_text(size=13,color="black"),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    legend.direction = "horizontal",
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    plot.title = element_text(size=15,color="black", hjust = 0.5, vjust = 1),
    legend.background = element_rect(color = "black", linewidth = 0.5, linetype = "solid"))
dev.off()

