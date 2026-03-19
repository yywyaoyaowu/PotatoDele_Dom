
setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/06_LD")

pdf("LDdecayResult.pdf", width = 5, height = 5)

read.table("LDdecayResult.CandolleanumPre")->ECandolleanumPre;
plot(ECandolleanumPre[,1]/1000, ECandolleanumPre[,2], type="l", col="#F3C97F", cex=1, main="", xlab="Distance(Kb)",
     xlim=c(0,300), ylim=c(0,0.4), ylab=expression(r^{2}), bty="n", lwd=2, cex.lab=1.2, cex.axis=1.2)

read.table("LDdecayResult.Candolleanum")->ECandolleanum;
lines(ECandolleanum[,1]/1000, ECandolleanum[,2], col="#EE9A00", lwd=2)

read.table("LDdecayResult.Landrace")->ELandrace;
lines(ELandrace[,1]/1000, ELandrace[,2],col="#548B54", lwd=2)

legend("topright",c("Candolleanum","Previous Candolleanum","Landrace"),col=c("#EE9A00","#F3C97F","#548B54"),
       cex=1.2,lty=c(1,1,1),bty="n",lwd=2);

dev.off()



setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/06_LD")

pdf("LDdecayResult_selective_region.pdf", width = 5, height = 5)

read.table("LDdecayResult.Landrace_Selective_region")->ESelective_region;
plot(ESelective_region[,1]/1000, ESelective_region[,2], type="l", col="#1C60AC", cex=1, main="", xlab="Distance(Kb)",
     xlim=c(0,300), ylim=c(0,0.5), ylab=expression(r^{2}), bty="n", lwd=2, cex.lab=1.2, cex.axis=1.2)

read.table("LDdecayResult.Landrace_Whole_genome")->EWhole_genome;
lines(EWhole_genome[,1]/1000, EWhole_genome[,2], col="#D3632C", lwd=2)

legend("topright",c("Selective region","Whole genome"),col=c("#1C60AC","#D3632C"),
       cex=1.2,lty=c(1,1,1),bty="n",lwd=2);

dev.off()


