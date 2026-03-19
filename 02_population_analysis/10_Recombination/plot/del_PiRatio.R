library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)
options(scipen = 0)

Args <- commandArgs(T)
del_input <- Args[1]
dom_input <- Args[2]
snp_input <- Args[3]

del <- read.table(del_input, header=T,sep="\t")
Lan_del <- del[del$Land_DeleFreq > 0 & del$GERP >= 2.75, ]
snp <- read.table(snp_input, header=F, sep = "\t")
colnames(snp) <- c("CHROM", "BIN_START", "BIN_END", "snp_num")
pi <- read.table(dom_input, header=T, sep = "\t")
pi <- pi[pi$Can_PI > 0.0001, ]

df_Lan_del <- Lan_del[, c(1,2,3)]
df_Lan_del$del <- 1
df_windows <- pi[, c(1,2,3)]
colnames(df_windows) <- c("CHROM", "BIN_START", "BIN_END")

setDT(df_windows)
setDT(df_Lan_del)

setkey(df_windows, CHROM, BIN_START, BIN_END)
setkey(df_Lan_del, Chrom, PosStart, position)

overlap_del_num <- foverlaps(
  df_windows, 
  df_Lan_del,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("Chrom", "PosStart", "position"),
  type = "any", 
  nomatch = NA)

overlap_del <- overlap_del_num %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    del_num = sum(del, na.rm = TRUE),
    .groups = "drop")

pi$CHROM_BIN_START_BIN_END <- paste(pi$CHROM, pi$BIN_START, pi$BIN_END, sep="_")
overlap_del$CHROM_BIN_START_BIN_END <- paste(overlap_del$CHROM, overlap_del$BIN_START, overlap_del$BIN_END, sep="_")

plot_data <- inner_join(pi, overlap_del, by="CHROM_BIN_START_BIN_END")
plot_data <- plot_data[, c(1,2,3,6,11)]
colnames(plot_data) <- c("CHROM", "BIN_START", "BIN_END", "PiRation", "del_num")
plot_data$del_num_perKb <- plot_data$del_num / 100
#plot_data <- plot_data[plot_data$del_num > 0 & plot_data$PiRation > 0, ]
plot_data <- na.omit(plot_data)

### dSNPs/Kb
pdf("Lan_PiRatio_del_PerKb.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = PiRation, y = del_num_perKb)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x, 
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$PiRation, na.rm = TRUE) / 3,
    label.y = max(plot_data$del_num_perKb, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Pi Ratio",
    y = "Deleterious density (dSNPs/Kb)",
    title = "Landrace Pi Ratio vs. Deleterious Density"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.text.y.right = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
  )
dev.off()

### dSNPs/SNP
plot_data$CHROM_BIN_START_BIN_END <- paste(plot_data$CHROM, plot_data$BIN_START, plot_data$BIN_END, sep="_")
snp$CHROM_BIN_START_BIN_END <- paste(snp$CHROM, snp$BIN_START, snp$BIN_END, sep="_")

snp <- snp[snp$CHROM_BIN_START_BIN_END %in% plot_data$CHROM_BIN_START_BIN_END, ]

plot_data$snp_num <- snp$snp_num
plot_data <- as.data.frame(plot_data)
plot_data$del_snp <- plot_data$del_num / plot_data$snp_num
plot_data <- na.omit(plot_data)

pdf("Lan_PiRatio_del_PerSNP.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = PiRation, y = del_snp)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$PiRation, na.rm = TRUE) / 3,
    label.y = max(plot_data$del_snp, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Pi Ratio",
    y = "Deleterious density (dSNPs/SNP)",
    title = "Landrace Pi Ratio vs. Deleterious Density"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.text.y.right = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
  )
dev.off()


