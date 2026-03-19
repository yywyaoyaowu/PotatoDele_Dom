library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

Args <- commandArgs(T)
input <- Args[1]
output <- Args[2]

setwd("/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/dSNP_density")
input <- "DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.plot.bed"
output <- "DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.pdf"

your_data <- read.table(input, header = T)
colnames(your_data) <- c("Chromosome", "Start", "End", "A626", "E463", "Hotspot_AllDele_CNDMin")
your_data$Chromosome <- as.character(your_data$Chromosome)

your_data <- your_data %>%
  mutate(Mid = (Start + End) / 2,
         Width = End - Start,
         Chromosome = factor(Chromosome))

df_long <- your_data %>%
  pivot_longer(
    cols = c("A626", "E463", "Hotspot_AllDele_CNDMin"),
    names_to = "Variable",
    values_to = "Value"
  )
df_long$Chromosome <- paste(df_long$Chromosome, df_long$Variable, sep = "_")

df_long$Mid_Mb <- df_long$Mid / 1e6
df_long$Width_Mb <- df_long$Width / 1e6

custom_order <- c("chr12_A626", "chr12_E463","chr12_Hotspot_AllDele_CNDMin",
                  "chr11_A626", "chr11_E463","chr11_Hotspot_AllDele_CNDMin",
                  "chr10_A626", "chr10_E463","chr10_Hotspot_AllDele_CNDMin", 
                  "chr09_A626", "chr09_E463","chr09_Hotspot_AllDele_CNDMin", 
                  "chr08_A626", "chr08_E463","chr08_Hotspot_AllDele_CNDMin", 
                  "chr07_A626", "chr07_E463","chr07_Hotspot_AllDele_CNDMin", 
                  "chr06_A626", "chr06_E463","chr06_Hotspot_AllDele_CNDMin", 
                  "chr05_A626", "chr05_E463","chr05_Hotspot_AllDele_CNDMin", 
                  "chr04_A626", "chr04_E463","chr04_Hotspot_AllDele_CNDMin", 
                  "chr03_A626", "chr03_E463","chr03_Hotspot_AllDele_CNDMin", 
                  "chr02_A626", "chr02_E463","chr02_Hotspot_AllDele_CNDMin", 
                  "chr01_A626", "chr01_E463","chr01_Hotspot_AllDele_CNDMin")

existing_chr <- intersect(custom_order, unique(df_long$Chromosome))
df_long$Chromosome <- factor(df_long$Chromosome, levels = existing_chr)

chromosome_bounds <- read.table("Solanum_tuberosumDM.bed.txt", header = F)
colnames(chromosome_bounds) <- c("Chromosome", "Start_Mb", "End_Mb")
chromosome_bounds_plot <- chromosome_bounds %>%
  # 为每个组创建副本
  tidyr::expand_grid(group = c("A626", "E463", "Hotspot_AllDele_CNDMin")) %>%
  # 计算中点和宽度
  mutate(
    Mid_Mb = (Start_Mb + End_Mb) / 2,
    Width_Mb = End_Mb - Start_Mb,
    # 创建组合的染色体名称
    Chromosome = paste(Chromosome, group, sep = "_")
  ) %>%
  arrange(group, Chromosome)

pdf(output, width = 10, height = 6)
ggplot(df_long, aes(x = Mid_Mb, y = Chromosome, fill = Value)) +
  geom_tile(aes(width = Width_Mb), height = 0.8) +
  geom_tile(
    data = chromosome_bounds_plot,
    aes(x = Mid_Mb, y = Chromosome, width = Width_Mb),
    height = 0.8,
    fill = NA,
    color = "gray",
    linewidth = 0.8) +
  scale_x_continuous(
    labels = function(x) paste0(x, "Mb"),
    breaks = seq(0, max(df_long$Mid_Mb)* 1.05, by = 10),
    limits = c(0, max(df_long$Mid_Mb) * 1.05),  # 增加5%的边距
    expand = expansion(mult = c(0, 0.05))       # 控制轴扩展
  ) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Genomic Position (Mb)", y = "Chromosome",
       fill = "Deleterious Burden", 
       title = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 10))
  )
dev.off()









