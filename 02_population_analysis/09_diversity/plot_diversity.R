library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/09_diversity")

genome.pi.CanPre <- read.table("CanPre.windowed.pi", header=T,sep="\t")
genome.pi.Can <- read.table("Can.windowed.pi", header=T,sep="\t")
genome.pi.Lan <- read.table("Lan.windowed.pi", header=T,sep="\t")

chr.length <- read.table("Solanum_tuberosumDM.length.txt", header=T, sep="\t")
colnames(chr.length) <- c("CHROM", "length", "chr")

colnames(genome.pi.CanPre) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.CanPre")
colnames(genome.pi.Can) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Can")
colnames(genome.pi.Lan) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Lan")

genome.pi.CanMerge <- inner_join(genome.pi.CanPre, genome.pi.Can, by="CHROM_BIN_START_BIN_END")
genome.pi.CanMerge <- genome.pi.CanMerge[, -c(5,7,8,9,10)]
colnames(genome.pi.CanMerge) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "pi.CanPre", "pi.Can")

genome.pi <- inner_join(genome.pi.CanMerge, genome.pi.Lan, by="CHROM_BIN_START_BIN_END")
genome.pi <- genome.pi[, -c(1,7,8,9,10)]
colnames(genome.pi) <- c("CHROM", "BIN_START", "BIN_END", "pi.CanPre", "pi.Can", "pi.Lan")

plot_data.CanPre <- data.frame(c(genome.pi$pi.CanPre), rep("Previous candolleanum", length(c(genome.pi$pi.CanPre))))
colnames(plot_data.CanPre) <- c("Pi", "Group")
plot_data.Can <- data.frame(c(genome.pi$pi.Can), rep("Candolleanum", length(c(genome.pi$pi.Can))))
colnames(plot_data.Can) <- c("Pi", "Group")
plot_data.Lan <- data.frame(c(genome.pi$pi.Lan), rep("Landrace", length(c(genome.pi$pi.Lan))))
colnames(plot_data.Lan) <- c("Pi", "Group")
plot_data <- rbind(plot_data.Lan, plot_data.Can, plot_data.CanPre)

cbbPalette <- c("#000000", "orange2", "#F3C97F", "palegreen4", "#56B4E9", "#8CC63E", "#009E73", 
                "#0072B2",  "#CC79A7", "#8CC63E", "#94CDEE", "#ECA2CF", "#D55E00")
col_key <- as.data.frame(t(cbbPalette[c(2,3,4)]))
colnames(col_key) <- c("Candolleanum", "Previous candolleanum", "Landrace")

pdf("diversity_compare.pdf", width = 2.1, height = 1.9)
ggplot(plot_data, aes(x=Group, y=Pi, fill=Group)) + geom_boxplot(outliers = F, lwd = 0.2) + 
  theme_bw(base_family = "serif") + scale_x_discrete(limits = c("Candolleanum", "Previous candolleanum", "Landrace")) + 
  scale_fill_manual(values = c("Previous candolleanum" = "#F3C97F",
                               "Candolleanum" = "orange2",
                               "Landrace" = "palegreen4")) + 
  theme(panel.grid.major=element_line(size=0),
        panel.grid.minor=element_line(size=0),
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(size=6),
        axis.text.x=element_text(size=4),
        axis.text.y=element_text(size=6),
        legend.position = "none") +
  ylim(0, 0.018) + xlab("") + ylab("Pi")
dev.off()
