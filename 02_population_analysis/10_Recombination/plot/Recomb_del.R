library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

Args <- commandArgs(T)
del_input <- Args[1]
snp_input <- Args[2]
recomb_input <- Args[3]

### input
del <- read.table(del_input, header=T,sep="\t")
Lan_del <- del[del$Land_DeleFreq > 0 & del$GERP >= 2.75, ]
snp <- read.table(snp_input, header=F, sep = "\t")
colnames(snp) <- c("CHROM", "BIN_START", "BIN_END", "snp_num")
Recomb <- read.table(recomb_input, header=T,sep="\t")
Recomb$recombRate <- Recomb$recombRate * 100000000
#Recomb <- Recomb[Recomb$recombRate != 0, ]

df_recomb <- Recomb
df_Lan_del <- Lan_del[, c(1,2,3)]
df_Lan_del$del <- 1
df_windows <- Recomb[, c(1,2,3)]
colnames(df_windows) <- c("CHROM", "BIN_START", "BIN_END")

setDT(df_windows)
setDT(df_recomb)
setDT(df_Lan_del)

setkey(df_recomb, chrom, start, end)
setkey(df_windows, CHROM, BIN_START, BIN_END)
setkey(df_Lan_del, Chrom, PosStart, position)

overlap_Recomb <- foverlaps(
  df_windows,
  df_recomb,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("chrom", "start", "end"),
  type = "any", 
  nomatch = NA)
overlap_del_num <- foverlaps(
  df_windows, 
  df_Lan_del,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("Chrom", "PosStart", "position"),
  type = "any", 
  nomatch = NA)

overlap_mean_Recomb <- overlap_Recomb %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    mean_recomb = mean(recombRate, na.rm = TRUE),
    .groups = "drop")
overlap_del <- overlap_del_num %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    del_num = sum(del, na.rm = TRUE),
    .groups = "drop")

overlap_mean_Recomb$CHROM_BIN_START_BIN_END <- paste(overlap_mean_Recomb$CHROM, overlap_mean_Recomb$BIN_START, overlap_mean_Recomb$BIN_END, sep="_")
overlap_del$CHROM_BIN_START_BIN_END <- paste(overlap_del$CHROM, overlap_del$BIN_START, overlap_del$BIN_END, sep="_")

plot_data <- inner_join(overlap_mean_Recomb, overlap_del, by="CHROM_BIN_START_BIN_END")
plot_data <- as.data.frame(plot_data)
plot_data <- plot_data[, c(1,2,3,4,9)]
colnames(plot_data) <- c("CHROM", "BIN_START", "BIN_END", "recombRate", "del_num")
plot_data$del_num_perKb <- plot_data$del_num / ((plot_data$BIN_END - plot_data$BIN_START) / 1000)

### dSNPs/Kb
pdf("Lan_Recomb_del_prelength.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = recombRate, y = del_num_perKb)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x, 
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$recombRate, na.rm = TRUE) / 3,
    label.y = max(plot_data$del_num_perKb, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Recombination Rate (cM/Mb)",
    y = "Deleterious number (dSNPs/Kb)",
    title = "Landrace Recombination Rate vs. Deleterious number"
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

plot_data$CHROM_BIN_START_BIN_END <- paste(plot_data$CHROM, plot_data$BIN_START, plot_data$BIN_END, sep="_")
snp$CHROM_BIN_START_BIN_END <- paste(snp$CHROM, snp$BIN_START, snp$BIN_END, sep="_")

snp <- snp[snp$CHROM_BIN_START_BIN_END %in% plot_data$CHROM_BIN_START_BIN_END, ]

plot_data$snp_num <- snp$snp_num
plot_data <- as.data.frame(plot_data)
plot_data$del_snp <- plot_data$del_num / plot_data$snp_num
plot_data <- na.omit(plot_data)

### dSNPs/SNP
pdf("Lan_Recomb_del_per_snp.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = recombRate, y = del_snp)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x, 
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$recombRate, na.rm = TRUE) / 3,
    label.y = max(plot_data$del_snp, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Recombination Rate (cM/Mb)",
    y = "Deleterious density (dSNPs/SNP)",
    title = "Landrace Recombination Rate vs. Deleterious density"
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


plot_data$snp_num_perKb <- plot_data$snp_num / (plot_data$BIN_END - plot_data$BIN_START)
### Recomb and snpDensity
pdf("Lan_Recomb_snpDensity.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = recombRate, y = snp_num_perKb)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$recombRate, na.rm = TRUE) / 3,
    label.y = max(plot_data$snp_num_perKb, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Recombination Rate (cM/Mb)",
    y = "SNP density (SNPs/Kb)",
    title = "Landrace Recombination Rate vs. SNP density"
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


### Deleterious density and snpDensity
pdf("Lan_del_snpDensity.pdf", width = 8, height = 8)
ggplot(plot_data, aes(x = del_num_perKb, y = snp_num_perKb)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$del_num_perKb, na.rm = TRUE) / 3,
    label.y = max(plot_data$snp_num_perKb, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Deleterious density (dSNPs/Kb)",
    y = "SNP density (SNPs/Kb)",
    title = "Landrace Deleterious density vs. SNP density"
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




