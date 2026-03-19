library(ggplot2)
library(scatterplot3d)
setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/08_PCA")

options(digits = 3)
eigenvel=read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile.pca338.eigenval")
PC1_explain <- round((eigenvel$V1[1] / sum(eigenvel$V1)) * 100, digits = 2)
PC2_explain <- round((eigenvel$V1[2] / sum(eigenvel$V1)) * 100, digits = 2)

df=read.table("eigenvec_pca338", header=T)
df <- df[df$Group != "inbred" & df$Group != "founder", ]

df[df$Group == "Previous_candolleanum", ]$Group <- 'Previous_candolleanum'
df[df$Group == "New_candolleanum", ]$Group <- 'New_candolleanum'
df[df$Group == "stenotomum", ]$Group <- 'Landrace'
df[df$Group == "goniocalyx", ]$Group <- 'Landrace'
df[df$Group == "phureja", ]$Group <- 'Landrace'
df[df$Group == "ajanhuiri", ]$Group <- 'Landrace'

df$Group<-factor(df$Group,levels=c("Previous_candolleanum","New_candolleanum","Landrace"))

cbbPalette <- c("#000000", "#F3C97F", "#F0A419", "palegreen4", "#8CA0CB", "#56B4E9", "#8CC63E", "#009E73", 
                "#0072B2",  "#CC79A7", "#8CC63E", "#94CDEE", "#ECA2CF", "#D55E00")
cbbPalette4=c(cbbPalette[c(2,3,4)])

pdf("PCA_2.pdf", width = 3.45, height = 1.9)
ggplot(data=df, aes(PC1, PC2)) +
  geom_point(aes(color=Group), alpha=0.9, size=0.3) +
  theme_bw(base_family = "serif") + scale_color_manual(values = cbbPalette4) +
  theme(panel.grid = element_blank(),
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6)) +
  xlab(paste0("PC1 (", PC1_explain, "%)")) + ylab(paste0("PC2 (", PC2_explain, "%)")) +
  guides(color = guide_legend(override.aes = list(size = 1)))
dev.off()
