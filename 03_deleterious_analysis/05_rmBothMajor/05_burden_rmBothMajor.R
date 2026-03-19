input<-read.table("Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_SolMsaNonMajorDele_Genotype_GERP2.75.txt.filter.txt",header=T)
#freq change
input$DeleIncrease<-input$Lan_DeleFreq-input$Can_DeleFreq
#select region
select_0.05<-read.table("Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.merge.bed",header=F)
colnames(select_0.05)<-c("chrom","start","end")
input$SelectiveRegion0.05 <- apply(input, 1, function(x) {
  chr <- as.character(x[1]) 
  pos <- as.numeric(x[3]) 
  any(select_0.05$chrom == chr & pos >= select_0.05$start & pos <= select_0.05$end)
})

select_0.02<-read.table("Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12.merge.bed",header=F)
colnames(select_0.02)<-c("chrom","start","end")
input$SelectiveRegion0.02 <- apply(input, 1, function(x) {
  chr <- as.character(x[1])
  pos <- as.numeric(x[3])
  any(select_0.02$chrom == chr & pos >= select_0.02$start & pos <= select_0.02$end)
})

select_0.1<-read.table("Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12.merge.bed",header=F)
colnames(select_0.1)<-c("chrom","start","end")
input$SelectiveRegion0.1 <- apply(input, 1, function(x) {
  chr <- as.character(x[1])
  pos <- as.numeric(x[3])
  any(select_0.1$chrom == chr & pos >= select_0.1$start & pos <= select_0.1$end)
})
#fisher_Pvalue and p.adjust_fdr
Prop_LanDeleNum=input$Lan_DeleHome*2+input$Lan_DeleHeter
Prop_CanDeleNum=input$Can_DeleHome*2+input$Can_DeleHeter
Lan_GenoNum=input$Lan_Geno*2
Progenitor_GenoNum=input$Can_Geno*2
Pro_test=data.frame(Prop_LanDeleNum,Prop_CanDeleNum,Lan_GenoNum,Progenitor_GenoNum)
Pro_test$fisher_Pvalue=apply(Pro_test,1,function(x) fisher.test(matrix(c(x[1],x[2],x[3],x[4]),2))$p.value)
input$fisher_Pvalue=Pro_test$fisher_Pvalue
input$p.adjust_fdr=p.adjust(input$fisher_Pvalue, method = "fdr")
#rmBothMajor
input_rmBothMajor=input[!(input$Lan_DeleFreq > 0.5 & input$Can_DeleFreq > 0.5) & !is.na(input$Lan_DeleFreq) & !is.na(input$Can_DeleFreq) & (input$Lan_DeleFreq >0 | input$Can_DeleFreq >0),]
write.table(input_rmBothMajor, 
            file = "Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_SolMsaNonMajorDele_Genotype_GERP2.75.txt.filter.rmBothMajor.txt",
            sep = "\t",
            row.names = FALSE,          
            col.names = TRUE,           
            quote = FALSE)     
