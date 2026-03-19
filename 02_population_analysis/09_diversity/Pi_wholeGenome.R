library(ggplot2)
library(dplyr)
setwd("/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/09_diversity")

genome.pi.Lan <- read.table("Lan.windowed.pi", header=T, sep="\t")
genome.pi.Can <- read.table("Can.windowed.pi", header=T, sep="\t")
genome.pi.CanPre <- read.table("CanPre.windowed.pi", header=T, sep="\t")

colnames(genome.pi.Lan) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Lan")
colnames(genome.pi.Can) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.Can")
colnames(genome.pi.CanPre) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi.CanPre")

genome.pi <- inner_join(genome.pi.Lan, genome.pi.Can, by="CHROM_BIN_START_BIN_END")
genome.pi <- genome.pi[, -c(5,7,8,9,10)]
colnames(genome.pi) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "pi.Lan", "pi.Can")

genome.pi <- inner_join(genome.pi, genome.pi.CanPre, by="CHROM_BIN_START_BIN_END")
genome.pi <- genome.pi[, -c(7,8,9,10)]
colnames(genome.pi) <- c("CHROM_BIN_START_BIN_END", "CHROM", "BIN_START", "BIN_END", "pi.Lan", "pi.Can", "pi.CanPre")

t.test(genome.pi$pi.CanPre, genome.pi$pi.Can)$p.value

all.chrs <- as.character(unique(genome.pi$CHROM))
chr.length=read.table("Solanum_tuberosumDM.length.txt", header=T, sep="\t")
chr.length$chr_asscum_len=c(0,cumsum(chr.length$length)[-nrow(chr.length)])
options(scipen = 200)
genome.pi$pos <- (genome.pi$BIN_START + genome.pi$BIN_END - 1) / 2
genome.pi$asscum_pos <- genome.pi$pos +
  chr.length$chr_asscum_len[match(genome.pi$CHROM,chr.length$DM_V6_chr)]
chr.mid.ascum=0.5*chr.length$length+chr.length$chr_asscum_len

win.len = 2000000
asscum_len = genome.pi$asscum_pos[length(genome.pi$asscum_pos)]
win.num = floor(asscum_len/win.len)
win.start=seq(0, asscum_len, win.len)
win.end = win.start + win.len
win.infor = data.frame(win.start, win.end)
win.infor[, 3:5] = NA
for (n in 1:length(win.end)){
  one.window = subset(genome.pi, asscum_pos>=win.start[n] & asscum_pos<win.end[n])
  win.infor[n, 3:5] = apply(one.window[, c(5,6,7)], 2, mean)
}
names(win.infor)[3:5] = names(one.window[,c(5,6,7)])
win.infor$mid.pos = (win.infor$win.start + win.infor$win.end) / 2

win.infor$pi.Lan <- win.infor$pi.Lan * 1000
win.infor$pi.Can <- win.infor$pi.Can * 1000
win.infor$pi.CanPre <- win.infor$pi.CanPre * 1000
pi.max=max(c(win.infor$pi.Lan, win.infor$pi.Can),na.rm=T)
pdf(paste("Pi_WholeGenome_2Mb.pdf",sep=""),width = 6, height = 2)
par(mar=c(2,3,1,1))
plot(win.infor$mid.pos/10000000,win.infor$pi.Lan,main="",pch=10,bg="grey",xlab="", ylab="",cex=1,font=2,col="palegreen4",
     type="l",cex.lab=1,font.lab=1,cex.axis=0.6,lwd=0.6,ylim=c(0,pi.max), xaxt="n")
lines(win.infor$mid.pos/10000000,win.infor$pi.Can,pch=10,cex=0.2,font=2,type="l",col="orange2",lwd=0.6)
lines(win.infor$mid.pos/10000000,win.infor$pi.CanPre,pch=10,cex=0.2,font=2,type="l",col="#F3C97F",lwd=0.6)
mtext(expression(pi~(x~10^-3)), side=2, at=0.5*pi.max, line=2, cex=0.8, col="black", font=2)
text(x = c(4.43, 11.1, 16.735, 23.15675, 29, 35, 41, 47 ,53, 59.5, 65, 70), y=-1.2,xpd=T,
     labels = c('Chr01','Chr02','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09','Chr10','Chr11','Chr12'), cex=0.6, font=2)
axis(side=1,at=c(0.000000, 8.859169, 13.469460, 19.540217, 26.463850, 32.023820, 37.932978, 43.696909, 49.619509, 56.379539, 62.483955, 67.161693, 73.1166933),
     labels = c("", "", "", "", "", "", "", "", "", "", "", "", ""), lwd=1, line=0)
#legend("bottomright",c("All Candolleanum", "Previous Candolleanum", "Landrace"),bty="n",lwd=2,cex=2,text.font =1.5,ncol=1,xpd=T,
#       lty = c(1, 1, 1), col=c("orange2", "#F3C97F", "palegreen4"))
box()
dev.off()
 
