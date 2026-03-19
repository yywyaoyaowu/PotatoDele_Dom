library(ggplot2)
library(dplyr)
setwd("/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/het")

het.Lan <- read.table("Lan.het.bed", header=F, sep="\t")
het.Can <- read.table("Can.het.bed", header=F, sep="\t")
colnames(het.Lan)=c("chr", "pos-1", "pos", "het")
colnames(het.Can)=c("chr", "pos-1", "pos", "het")

window_size <- 2000000
chr_order <- unique(het.Lan$chr)
chr_lengths <- het.Lan %>%
  group_by(chr) %>%
  summarise(max_pos = max(pos, na.rm = TRUE)) %>%
  arrange(factor(chr, levels = chr_order)) %>%
  mutate(cum_start = cumsum(lag(max_pos, default = 0)))

het.Lan.win <- het.Lan %>%
  mutate(window = floor(pos / window_size) * window_size) %>%
  group_by(chr, window) %>%
  summarise(mean_het = mean(het, na.rm = TRUE)) %>%
  left_join(chr_lengths %>% select(chr, cum_start), by = "chr") %>%
  mutate(cum_pos = cum_start + window) %>%
  select(-cum_start)

het.Can.win <- het.Can %>%
  mutate(window = floor(pos / window_size) * window_size) %>%
  group_by(chr, window) %>%
  summarise(mean_het = mean(het, na.rm = TRUE)) %>%
  left_join(chr_lengths %>% select(chr, cum_start), by = "chr") %>%
  mutate(cum_pos = cum_start + window) %>%
  select(-cum_start)

het.Lan.win$Group <- "Lan"
het.Can.win$Group <- "Can"

het.win <- rbind(het.Lan.win, het.Can.win)
het.max=max(c(het.win$mean_het),na.rm=T)

pdf(paste("het_WholeGenome_2Mb.pdf",sep=""),width = 6, height = 2)
par(mar=c(2,3,1,1))
plot(het.Lan.win$cum_pos/10000000,het.Lan.win$mean_het,main="",pch=10,bg="grey",xlab="", ylab="",cex=1,font=2,col="palegreen4",
     type="l",cex.lab=1,font.lab=1,cex.axis=0.6,lwd=0.6,ylim=c(0,het.max*1.2), xaxt="n")
lines(het.Can.win$cum_pos/10000000,het.Can.win$mean_het,pch=10,cex=0.2,font=2,type="l",col="orange2",lwd=0.6)
text(x = c(4.43, 11.1, 16.735, 23.15675, 29, 35, 41, 47 ,53, 59.5, 65, 70), y=-0.02,xpd=T,
     labels = c('Chr01','Chr02','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09','Chr10','Chr11','Chr12'), cex=0.6, font=2)
axis(side=1,at=c(0.000000, 8.859169, 13.469460, 19.540217, 26.463850, 32.023820, 37.932978, 43.696909, 49.619509, 56.379539, 62.483955, 67.161693, 73.1166933),
     labels = c("", "", "", "", "", "", "", "", "", "", "", "", ""), lwd=1, line=0)
legend("topright",c("Candolleanum", "Landrace"),bty="n",lwd=2,cex=0.6,text.font =1.5,ncol=1,xpd=T,
       lty = c(1, 1), col=c("orange2", "palegreen4"))
box()
dev.off()

