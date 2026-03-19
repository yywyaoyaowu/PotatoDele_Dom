library(ggplot2)
library(dplyr)
library(stringr)
myColor <- c("#AD1414", "#FFA12C", "#FFFE93", "#D3632C", "#200000",
             "#243285", "#4DA4BA", "#85BBDA", "#1C60AC", "#5C7FAD",
             "#9590CB", "#8361B3", "#432374", "#8A58C8", "#412850",
             "orange2", "#548B54")
Args <- commandArgs(T)
Lan_del_Conservation2 <- Args[1]
Lan_del_Conservation2.75 <- Args[2]
Lan_del_Conservation3.5 <- Args[3]
Can_del_Conservation2 <- Args[4]
Can_del_Conservation2.75 <- Args[5]
Can_del_Conservation3.5 <- Args[6]
output <- Args[7]

Lan_Conservation2 <- read.table(Lan_del_Conservation2,header=F)
Lan_Conservation2.75 <- read.table(Lan_del_Conservation2.75,header=F)
Lan_Conservation3.5 <- read.table(Lan_del_Conservation3.5,header=F)

colnames(Lan_Conservation2)=c("chr","pos","het")
colnames(Lan_Conservation2.75)=c("chr","pos","het")
colnames(Lan_Conservation3.5)=c("chr","pos","het")

Lan_Conservation2$Conservation <- "Lan_Conservation2"
Lan_Conservation2.75$Conservation <- "Lan_Conservation2.75"
Lan_Conservation3.5$Conservation <- "Lan_Conservation3.5"

Can_Conservation2 <- read.table(Can_del_Conservation2,header=F)
Can_Conservation2.75 <- read.table(Can_del_Conservation2.75,header=F)
Can_Conservation3.5 <- read.table(Can_del_Conservation3.5,header=F)

colnames(Can_Conservation2)=c("chr","pos","het")
colnames(Can_Conservation2.75)=c("chr","pos","het")
colnames(Can_Conservation3.5)=c("chr","pos","het")

Can_Conservation2$Conservation <- "Can_Conservation2"
Can_Conservation2.75$Conservation <- "Can_Conservation2.75"
Can_Conservation3.5$Conservation <- "Can_Conservation3.5"

dele_het <- rbind(Can_Conservation2, Can_Conservation2.75, Can_Conservation3.5,
                  Lan_Conservation2, Lan_Conservation2.75, Lan_Conservation3.5)
Lan_dele_het <- dele_het[dele_het$Conservation %in% c("Lan_Conservation2", "Lan_Conservation2.75", "Lan_Conservation3.5"), ]
Can_dele_het <- dele_het[dele_het$Conservation %in% c("Can_Conservation2", "Can_Conservation2.75", "Can_Conservation3.5"), ]

Lan_dele_het <- Lan_dele_het %>%
  mutate(het_interval = cut(het, 
                            breaks = seq(0, 1, by = 0.05),
                            include.lowest = TRUE))

Can_dele_het <- Can_dele_het %>%
  mutate(het_interval = cut(het, 
                            breaks = seq(0, 1, by = 0.05),
                            include.lowest = TRUE))

Lan_dele_het_Conservation_proportion <- Lan_dele_het %>%
  group_by(Conservation) %>%
  mutate(het_interval = cut(het, 
                            breaks = seq(0, 1, by = 0.2),
                            include.lowest = TRUE)) %>%
  group_by(Conservation, het_interval) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Conservation) %>%
  mutate(proportion = count / sum(count)) %>% 
  ungroup()

Can_dele_het_Conservation_proportion <- Can_dele_het %>%
  group_by(Conservation) %>%
  mutate(het_interval = cut(het, 
                            breaks = seq(0, 1, by = 0.2),
                            include.lowest = TRUE)) %>%
  group_by(Conservation, het_interval) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Conservation) %>%
  mutate(proportion = count / sum(count)) %>% 
  ungroup()

Lan_dele_het_Conservation_proportion$Group_het_interval <- paste("Lan", Lan_dele_het_Conservation_proportion$het_interval, sep = "_")
Can_dele_het_Conservation_proportion$Group_het_interval<- paste("Can", Can_dele_het_Conservation_proportion$het_interval, sep = "_")
dele_het_Conservation_proportion <- rbind(Lan_dele_het_Conservation_proportion, Can_dele_het_Conservation_proportion)

dele_het_Conservation_proportion$Conservation <- factor(
  dele_het_Conservation_proportion$Conservation,
  levels = c("Lan_Conservation2", "Lan_Conservation2.75", "Lan_Conservation3.5",
             "Can_Conservation2", "Can_Conservation2.75", "Can_Conservation3.5"))
dele_het_Conservation_proportion$Group_het_interval <- factor(
  dele_het_Conservation_proportion$Group_het_interval,
  levels = c("Lan_[0,0.2]", "Lan_(0.2,0.4]", "Lan_(0.4,0.6]", "Lan_(0.6,0.8]", "Lan_(0.8,1]",
             "Can_[0,0.2]", "Can_(0.2,0.4]", "Can_(0.4,0.6]", "Can_(0.6,0.8]", "Can_(0.8,1]"))

dele_het_Conservation_proportion$Group <- str_split(dele_het_Conservation_proportion$Conservation, pattern = "_", simplify = TRUE)[,1]
dele_het_Conservation_proportion$Conservation <- str_split(dele_het_Conservation_proportion$Conservation, pattern = "_", simplify = TRUE)[,2]

dele_het_Conservation_proportion$Conservation <- factor(
  dele_het_Conservation_proportion$Conservation,
  levels = c("Conservation2", "Conservation2.75", "Conservation3.5"))
dele_het_Conservation_proportion$Group <- factor(
  dele_het_Conservation_proportion$Group,
  levels = c("Can", "Lan"))


pdf(output, width=6, height=6)
ggplot(dele_het_Conservation_proportion, aes(x = Group, y = proportion, fill = het_interval)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Conservation, nrow = 1) +
  labs(x = "",
       y = "Proportion (%)",
       fill = "Heterozygous deleterious ratio") +
  scale_fill_manual(values = c("[0,0.2]" = "#85BBDA",
                               "(0.2,0.4]" = "#4DA4BA", 
                               "(0.4,0.6]" = "#5C7FAD",
                               "(0.6,0.8]" = "#1C60AC",
                               "(0.8,1]" = "#243285")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right")
dev.off()

write.table(dele_het_Conservation_proportion, "dele_het_Conservation_proportion.txt", row.names = F, quote = F)
