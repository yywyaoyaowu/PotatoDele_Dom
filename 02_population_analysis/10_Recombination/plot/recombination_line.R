library(ggplot2)
library(dplyr)

Args <- commandArgs(T)
Recombination <- Args[1]
chr_length <- Args[2]
prefix <- Args[3]

data <- read.table(Recombination, header = T)
data$pos <- data$start + (data$end - data$start) / 2
data$recombRate <- data$recombRate * 100000000

threshold_top5 <- quantile(data$recombRate, 0.95, na.rm = TRUE)
data$is_top5 <- data$recombRate >= threshold_top5

chr.length = read.table(chr_length, header = T, sep = "\t")
chr.length$chr_asscum_len = c(0, cumsum(chr.length$length)[-nrow(chr.length)])
data$asscum_pos <- data$pos +
  chr.length$chr_asscum_len[match(data$chr, chr.length$DM_V6_chr)]
chr_boundaries <- data %>%
  group_by(chrom) %>%
  summarise(
    start = min(asscum_pos),
    end = max(asscum_pos),
    mid = mean(asscum_pos))

chr_boundaries$chr_num <- as.numeric(gsub("chr", "", chr_boundaries$chrom))
data$chr_num <- as.numeric(gsub("chr", "", data$chrom))
odd_color <- "#1C60AC"
even_color <- "#5C7FAC"
data$chr_color <- ifelse(data$chr_num %% 2 == 1, odd_color, even_color)

pdf(paste0(prefix, "_recombination.pdf"), width = 6.8, height = 2)
ggplot() + geom_col(data = data, aes(x = asscum_pos, y = recombRate, fill = chr_color), 
                    color = data$chr_color, alpha = 0.8, linewidth = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10), name = "Recombination Rate (cM/Mb)") +
  scale_x_continuous(expand = c(0, 0), breaks = c(chr_boundaries$mid), labels = c(chr_boundaries$chrom)) +
  geom_hline(yintercept = 2.15, linetype = "dashed") +
  theme_bw() + xlab("") + labs(title = prefix) +
  theme(
	plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8))
dev.off() 



