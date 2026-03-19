library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

Args <- commandArgs(T)
pi_input <- Args[1]
Recomb_input <- Args[2]
bed_input <- Args[3]
prefix <- Args[4]

genome.pi <- read.table(pi_input, header=T,sep="\t")
gene_bed <- read.table(bed_input, header=F,sep="\t")
colnames(gene_bed) <- c("chrom", "start", "end", "geneID", "a", "strand")
gene_bed$base_num <- gene_bed$end - gene_bed$start + 1
Recomb <- read.table(Recomb_input, header=T,sep="\t")
#Recomb <- Recomb[Recomb$recombRate != 0, ]
Recomb$recombRate <- Recomb$recombRate * 100000000

df_windows <- genome.pi[, c(2,3,4)]
df_recomb <- Recomb
df_gene <- gene_bed

# 转换为data.table进行高效区间计算
setDT(df_windows)
setDT(df_recomb)
setDT(df_gene)

# 设置区间键
setkey(df_windows, CHROM, BIN_START, BIN_END)
setkey(df_recomb, chrom, start, end)
setkey(df_gene, chrom, start, end)

# 执行区间重叠计算
overlap_gene_num <- foverlaps(
  df_windows, 
  df_gene,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("chrom", "start", "end"),
  type = "any", 
  nomatch = NA)
overlap_recomb <- foverlaps(
  df_windows, 
  df_recomb,
  by.x = c("CHROM", "BIN_START", "BIN_END"),
  by.y = c("chrom", "start", "end"),
  type = "any", 
  nomatch = NA)

overlap_num <- overlap_gene_num %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    gene_base_num = sum(base_num, na.rm = TRUE),
    .groups = "drop")
overlap_mean_Recomb <- overlap_recomb %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    mean_recomb = mean(recombRate, na.rm = TRUE),
    .groups = "drop")

overlap_mean_Recomb$CHROM_BIN_START_BIN_END <- paste(overlap_mean_Recomb$CHROM, overlap_mean_Recomb$BIN_START, overlap_mean_Recomb$BIN_END, sep="_")
overlap_num$CHROM_BIN_START_BIN_END <- paste(overlap_num$CHROM, overlap_num$BIN_START, overlap_num$BIN_END, sep="_")

plot_data <- inner_join(overlap_mean_Recomb, overlap_num, by="CHROM_BIN_START_BIN_END")
plot_data <- plot_data[, c(1,2,3,4,9)]
colnames(plot_data) <- c("CHROM", "BIN_START", "BIN_END", "recombRate", "gene_base_num")
plot_data <- na.omit(plot_data)
plot_data$gene_base_num <- plot_data$gene_base_num / 100000
#plot_data <- plot_data[plot_data$recombRate != 0, ]

pdf(paste(prefix, "Recomb_geneDensity.pdf", sep = '_'), width = 8, height = 8)
ggplot(plot_data, aes(x = recombRate, y = gene_base_num)) +
  geom_point(alpha = 0.6, size = 2, color = "#4E79A7") +
  geom_smooth(
    method = "lm", formula = y ~ x, 
    se = TRUE, color = "#E15759", linewidth = 1.5
  ) +
  # 添加R²和p值标签（位置可调）
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",           # 方法可选pearson/spearman
    label.x = max(plot_data$recombRate, na.rm = TRUE) / 3,  # 左下角
    label.y = max(plot_data$gene_base_num, na.rm = TRUE),       # 顶部对齐
    size = 10,                    # 字体大小
    color = "black"              # 字体颜色
  ) +
  labs(
    x = "Recombination Rate (cM/Mb)",
    y = "Gene Density",
    title = paste(prefix, "Recombination Rate vs. Gene Density", sep = ' ')
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




