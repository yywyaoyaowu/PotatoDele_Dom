library(ggplot2)
library(dplyr)
library(data.table)
setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/09_diversity")

genome.pi.CanPre <- read.table("CanPre.windowed.pi", header=T,sep="\t")
genome.pi.Lan <- read.table("Lan.windowed.pi", header=T,sep="\t")
genome.pi.Can <- read.table("Can.windowed.pi", header=T,sep="\t")

chr.length <- read.table("Solanum_tuberosumDM.length.txt", header=T, sep="\t")
colnames(chr.length) <- c("CHROM", "length", "chr")

colnames(genome.pi.Lan) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Lan")
colnames(genome.pi.Can) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Can")

genome.pi <- inner_join(genome.pi.Lan, genome.pi.Can, by="CHROM_BIN_START_BIN_END")
genome.pi <- genome.pi[, -c(1,7,8,9,10)]
colnames(genome.pi) <- c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Lan", "pi.Can")
genome.pi$pi.ratio <- genome.pi$pi.Can / genome.pi$pi.Lan
genome.pi$mid.pos <- (genome.pi$BIN_START + genome.pi$BIN_END)/2

all.chrs=as.character(unique(genome.pi$CHROM))
chr.length$chr_asscum_len=c(0,cumsum(chr.length$length)[-nrow(chr.length)])
options(scipen = 200)
genome.pi$asscum_pos <- genome.pi$mid.pos+
  chr.length$chr_asscum_len[match(genome.pi$CHROM,chr.length$CHROM)]
chr.mid.ascum=0.5*chr.length$length+chr.length$chr_asscum_len

all.win.infor=NULL
win.len = 2000000

plot_data <- data.table()
for (i in 1:nrow(chr.length)){
  genome.pi_chr <- genome.pi[genome.pi$CHROM == chr.length$CHROM[i] & 
                               genome.pi$BIN_END >= 1 & 
                               genome.pi$BIN_START <= chr.length$length[i], ]
  pi.max=max(c(genome.pi_chr$pi.Lan, genome.pi_chr$pi.Can),na.rm=T)
  
  chr=chr.length$CHROM[i]
  win.num=floor(chr.length[i,2]/win.len)
  win.start=seq(0,chr.length[i,2],win.len)
  win.end=win.start+win.len
  win.infor=data.frame(chr,win.start,win.end)
  win.infor[,4:6]=NA
  
  for (n in 1:length(win.end)){
    one.window=subset(genome.pi_chr,mid.pos>=win.start[n] & mid.pos<win.end[n])
    win.infor[n,4:6]=apply(one.window[,c(5:7)],2,mean)
  }
  names(win.infor)[4:6]=names(one.window[,c(5:7)])
  win.infor$mid.pos=apply(win.infor[,c(2,3)],1,mean)
  
  all.win.infor=rbind(all.win.infor,win.infor)
  genome.pi_chr=win.infor
  plot_data <- rbind(plot_data, genome.pi_chr)
}

plot_data$asscum_pos <- plot_data$mid.pos+
  chr.length$chr_asscum_len[match(plot_data$chr,chr.length$CHROM)]
chr_centers <- plot_data %>%
  group_by(chr) %>%
  summarize(center = mean(range(asscum_pos)))  # 计算染色体中心位置

pdf(paste0("Pi_WholeGenomeAddRatio.pdf"), height = 5, width = 25)
p <- ggplot(plot_data) +
  geom_line(aes(x=asscum_pos/1000000, y=pi.Lan, color="Landrace")) +
  geom_line(aes(x=asscum_pos/1000000, y=pi.Can, color="Candolleanum")) +
  geom_line(aes(x=asscum_pos/1000000, y=pi.ratio/473.3639, color="Pi Ratio"), linewidth=1) +
  geom_line(aes(x=asscum_pos/1000000, y=2.23/473.3639, color="Pi Ratio = 2.23"), linetype=2, linewidth=1) +
  theme_bw() + scale_y_continuous(name = "Pi", sec.axis = sec_axis(~.*473.3639, name="Pi Ratio")) +
  scale_x_continuous(breaks = chr_centers$center/1000000,labels = chr_centers$chr) +
  scale_color_manual(values = c("Landrace"="palegreen4", "Candolleanum"="orange2", 
                                "Pi Ratio"="black", "Pi Ratio = 2.23" = "gray50"), name = "") +
  labs(x = "chromosome (Mb)", y = "Pi Ratio") + 
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
print(p)
dev.off()

all.win.infor=NULL
win.len = 500000


for (i in 1:nrow(chr.length)){
  genome.pi_chr <- genome.pi[genome.pi$CHROM == chr.length$CHROM[i] & 
                               genome.pi$BIN_END >= 1 & 
                               genome.pi$BIN_START <= chr.length$length[i], ]
  pi.max=max(c(genome.pi_chr$pi.Lan, genome.pi_chr$pi.Can),na.rm=T)
  
  chr=chr.length$CHROM[i]
  win.num=floor(chr.length[i,2]/win.len)
  win.start=seq(0,chr.length[i,2],win.len)
  win.end=win.start+win.len
  win.infor=data.frame(chr,win.start,win.end)
  win.infor[,4:6]=NA
  for (n in 1:length(win.end)){
    one.window=subset(genome.pi_chr,mid.pos>=win.start[n] & mid.pos<win.end[n])
    win.infor[n,4:6]=apply(one.window[,c(5:7)],2,mean)
  }
  names(win.infor)[4:6]=names(one.window[,c(5:7)])
  win.infor$mid.pos=apply(win.infor[,c(2,3)],1,mean)
  
  all.win.infor=rbind(all.win.infor,win.infor)
  genome.pi_chr=win.infor
  
  pdf(paste0("Pi_",chr.length$CHROM[i],".pdf"), height = 5, width = 15)
  p <- ggplot(genome.pi_chr) +
    geom_line(aes(x=mid.pos/1000000, y=pi.Lan, color="Landrace")) +
    geom_line(aes(x=mid.pos/1000000, y=pi.Can, color="Candolleanum")) +
    geom_line(aes(x=mid.pos/1000000, y=pi.ratio/473.3639, color="Pi Ratio"), linewidth=1) +
    geom_line(aes(x=mid.pos/1000000, y=2.23/473.3639, color="Pi Ratio = 2.23"), linetype=2, linewidth=1) +
    theme_bw() + scale_y_continuous(name = "Pi", sec.axis = sec_axis(~.*473.3639, name="Pi Ratio")) +
    scale_color_manual(values = c("Landrace"="palegreen4", "Candolleanum"="orange2", 
                                  "Pi Ratio"="black", "Pi Ratio = 2.23" = "gray50"), name = "") +
    labs(x = paste0(chr.length$CHROM[i], "(Mb)"), y = "Pi Ratio") + 
    #ggtitle(paste0("Pi_",chr.length$CHROM[i])) +
    theme(plot.title = element_text(size = 30, hjust = 0.5),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  print(p)
  dev.off()
}


