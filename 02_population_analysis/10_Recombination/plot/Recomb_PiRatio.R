library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

Args <- commandArgs(T)
pi_input <- Args[1]
Recomb_input <- Args[2]
prefix <- Args[3]

genome.pi <- read.table(pi_input, header=T, sep="\t")
genome.pi <- genome.pi[genome.pi$Can_PI > 0.0001, ]
Recomb <- read.table(Recomb_input, header=T,sep="\t")
Recomb$recombRate <- Recomb$recombRate * 100000000
#Recomb <- Recomb[Recomb$recombRate != 0, ]

df_recomb <- Recomb
df_windows <- genome.pi[, c(1,2,3)]

setDT(df_recomb)
setDT(df_windows)

setkey(df_recomb, chrom, start, end)
setkey(df_windows, CHROM, BIN_START, BIN_END)

overlap <- foverlaps(
  df_windows, 
  df_recomb,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("chrom", "start", "end"),
  type = "any", 
  nomatch = NA
)

overlap_Recomb <- overlap %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    mean_recomb = mean(recombRate, na.rm = TRUE),
    .groups = "drop")

genome.pi$CHROM_BIN_START_BIN_END <- paste(genome.pi$CHROM, genome.pi$BIN_START, genome.pi$BIN_END, sep="_")
overlap_Recomb$CHROM_BIN_START_BIN_END <- paste(overlap_Recomb$CHROM, overlap_Recomb$BIN_START, overlap_Recomb$BIN_END, sep="_")

plot_data <- inner_join(genome.pi, overlap_Recomb, by="CHROM_BIN_START_BIN_END")
plot_data <- plot_data[, c(1,2,3,6,11)]
colnames(plot_data) <- c("CHROM", "BIN_START", "BIN_END", "PIRatio", "recombRate")
plot_data <- na.omit(plot_data)
#plot_data <- plot_data[plot_data$recombRate != 0 & plot_data$PIRatio != 0, ]

pdf(paste(prefix, "Recomb_PiRatio.pdf", sep = '_'), width = 8, height = 8)
ggplot(plot_data, aes(x = recombRate, y = PIRatio)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x, 
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = max(plot_data$recombRate, na.rm = TRUE) / 3,
    label.y = max(plot_data$PIRatio, na.rm = TRUE),
    size = 10,
    color = "black"
  ) +
  labs(
    x = "Recombination Rate (cM/Mb)",
    y = "PI Ratio",
    title = paste(prefix, "Recombination Rate vs. PI Ratio", sep = ' ')
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





