library(ggplot2)
library(reshape2)
library(phytools)
library(ape)
library(magrittr)
library(dendextend)
library(viridis)
library(dplyr)
library(phylogram)

setwd("~/文件/马铃薯重测序-有害突变/dom_del_final/07_structure")
me_tree <- read.tree("DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.increaseNode.newick")
label <- rev(me_tree$tip.label)
write.table(label, "sortName.txt", quote = F, row.names = F)

sample_group <- read.table("sample_info_final")
colnames(sample_group) <- c("samples", "Group")
sample_group <- sample_group[!(sample_group$Group == "inbred" | sample_group$Group == "founder"), ]
sample_group$col <- c("NA")
sample_group[sample_group$Group=="stenotomum_subsp_stenotomum", ]$col <- rep("#FF0000",dim(sample_group[sample_group$Group=="stenotomum_subsp_stenotomum", ])[1])
sample_group[sample_group$Group=="stenotomum_subsp.goniocalyx", ]$col <- rep("#F8B55A",dim(sample_group[sample_group$Group=="stenotomum_subsp.goniocalyx", ])[1])
sample_group[sample_group$Group=="phureja", ]$col <- rep("#B97EC0",dim(sample_group[sample_group$Group=="phureja", ])[1])
sample_group[sample_group$Group=="ajanhuiri", ]$col <- rep("#F701FD",dim(sample_group[sample_group$Group=="ajanhuiri", ])[1])

sample_group[sample_group$Group=="New_candolleanum", ]$col <- rep("#B6E05D",dim(sample_group[sample_group$Group=="New_candolleanum", ])[1])
sample_group[sample_group$Group=="Previous_candolleanum", ]$col <- rep("#87A7FF",dim(sample_group[sample_group$Group=="Previous_candolleanum", ])[1])
sample_group[sample_group$col=="NA", ]$col <- "#000000"

label <- label[!(label=="A626" | label=="E463" | label=="E8669" | 
                   label=="C10-20" | label=="C151" | label=="RH")]


#### K = 2 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.2.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample, df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]
s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#FFFF00", "#FF0000")  ###设置调色板

df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "2", sep = " ")

aql <- aql[!(aql$samples=="A626" | aql$samples=="E463" | aql$samples=="E8669" | 
               aql$samples=="C10-20" | aql$samples=="C151" | aql$samples=="RH"), ]

pdf("admixture_k2_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab) #+ xlab("")
dev.off()


#### K = 3 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.3.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample, df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]
s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#00FFFF","#FFFF00","#FF0000")  ###设置调色板

df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2", "V3")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "3", sep = " ")

aql <- aql[!(aql$samples=="A626" | aql$samples=="E463" | aql$samples=="E8669" | 
               aql$samples=="C10-20" | aql$samples=="C151" | aql$samples=="RH"), ]

pdf("admixture_k3_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab)
dev.off()



#### K = 4 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.4.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample, df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]
s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#FFFF00","#FF0000","#00FF00","#00FFFF")  ###设置调色板

df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2", "V3", "V4")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "4", sep = " ")

aql <- aql[!(aql$samples=="A626" | aql$samples=="E463" | aql$samples=="E8669" | 
               aql$samples=="C10-20" | aql$samples=="C151" | aql$samples=="RH"), ]

pdf("admixture_k4_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab)
dev.off()



#### K = 5 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.5.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample,df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]

s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#00FF00","#00FFFF","#FFFF00","#FF0000","#FF00FF")  ###设置调色板


df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2", "V3", "V4", "V5")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "5", sep = " ")

pdf("admixture_k5_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab)
dev.off()




#### K = 6 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.6.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample,df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]

s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#00FFFF","#0000FF","#FFFF00","#00FF00","#FF00FF","#FF0000")  ###设置调色板

df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2", "V3", "V4", "V5", "V6")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "6", sep = " ")

pdf("admixture_k6_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab)
dev.off()


#### K = 7 #######
df_order <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.7.Q", header = F, sep = " ")
df_sample <- read.table("DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.nosex", header = F, sep = "\t") ##存放个体名称的文件
df_sample  <- df_sample[ ,1]
df_order <- data.frame(df_sample,df_order)
df_order <- df_order[!(df_order$df_sample=="A626" | df_order$df_sample=="E463" | df_order$df_sample=="E8669" | 
                         df_order$df_sample=="C10-20" | df_order$df_sample=="C151" | df_order$df_sample=="RH"), ]

s <- data.frame(df_sample=c(label), order=c(1:336))
s <- merge(df_order, s)

df_sample <- label
s <- data.frame(samples=c(df_sample), order=c(1:336))
s <- merge(sample_group, s)
col_text_x <- s$col[order(s$order)]
col_panel <- c("#FF0000","#0000FF","#F4BA1B","#00FF00","#FF00FF","#00FFFF", "#FFFF00")  ###设置调色板

df_new <- df_order[order(s$order),]
colnames(df_new) <- c("samples", "V1", "V2", "V3", "V4", "V5", "V6")
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K = ", "6", sep = " ")

pdf("admixture_k7_2.pdf", width = 48, height = 8)
ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
  scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.title.y = element_text(size = 13),
        panel.grid = element_blank()) + ylab(y_lab)
dev.off()




axis.text.x = element_text(angle = 90, size = 10, colour = col_text_x),
