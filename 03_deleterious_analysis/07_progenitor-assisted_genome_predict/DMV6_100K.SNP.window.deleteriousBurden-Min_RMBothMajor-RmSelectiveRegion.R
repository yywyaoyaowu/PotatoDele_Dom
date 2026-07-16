chr.len=read.table("Input/Solanum_tuberosumDM.fa.fai",sep="\t")
Remove_SelectRegion=read.table("/data/dzhst/Potato/20260622_tpc/select/result5/Wholegenome_NoSelect_filter.txt")
names(Remove_SelectRegion)=c("Chrom","Start","End")
#GroupInfo=read.table("/home/wuyaoyao/04_Potato/04_DomDele/01_VCF/sample_info_noEtu_final")
GroupInfo=read.table("Input/sample_info_502")

names(GroupInfo)=c("AccessionID","GroupType")
S.can_ID=GroupInfo$AccessionID[GroupInfo$GroupType=="Candolleanum"]
Inbred_ID=GroupInfo$AccessionID[GroupInfo$GroupType=="Inbread"]
Landrace_ID=GroupInfo$AccessionID[GroupInfo$GroupType=="Landrace"]

options(scipen = 200)
win.len=100000
#step=200000
All.win.infor=NULL
chrs=as.character(unique(chr.len[1:12,1]))
Homo_Dele_Genome=NULL
Heter_Dele_Genome=NULL

All_Dele_Genome=NULL
for (i in 1:length(chrs)){
   VCF_Freq=read.table(paste0("Input/DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_",chrs[i],"_Freq"),skip=1)

   deleterious.infor_ALL=read.table(paste0("Input/DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_",chrs[i],"_ConservationCutoff2.75_SolMsaNonMajorDele_Genotype.txt"),header=T)
   deleterious.infor=deleterious.infor_ALL[!(deleterious.infor_ALL$Land_DeleFreq > 0.5 & deleterious.infor_ALL$SCand_DeleFreq > 0.5) & !is.na(deleterious.infor_ALL$Land_DeleFreq) & !is.na(deleterious.infor_ALL$SCand_DeleFreq) & (deleterious.infor_ALL$Land_DeleFreq >0 | deleterious.infor_ALL$SCand_DeleFreq >0),]
   dim(deleterious.infor)# 1261004      24

    Remove_SelectRegion_Chr=subset(Remove_SelectRegion,Chrom==chrs[i])

   for (seg in 1:nrow(Remove_SelectRegion_Chr)){
     VCF_Freq_Region=VCF_Freq[VCF_Freq[,1]==chrs[i] & VCF_Freq[,2] > Remove_SelectRegion_Chr$Start[seg] & VCF_Freq[,2] <= Remove_SelectRegion_Chr$End[seg],]

   WinNum=floor(nrow(VCF_Freq_Region)/win.len)
   win.star_XL=seq(0,WinNum)*win.len+1
   win.end_XL=win.star_XL+win.len-1

   SNP_num=win.len
   win.star=VCF_Freq_Region[win.star_XL,2]
   win.end=VCF_Freq_Region[win.end_XL,2]
   chr=rep(chrs[i],length(win.star))

   win.infor=data.frame(chr,win.star,win.end,SNP_num)
   win.infor[nrow(win.infor),3]=VCF_Freq_Region[nrow(VCF_Freq_Region),2]
    win.infor$SNP_num[nrow(win.infor)]=nrow(VCF_Freq_Region)%%win.len

  for (m in 1:nrow(win.infor)){
    print(paste(i,m))
    deleterious_window =deleterious.infor[deleterious.infor$position>=win.infor$win.star[m] & deleterious.infor$position<=win.infor$win.end[m],-c(1:20)]
    w_Conserved=deleterious.infor$GERP[deleterious.infor$position>=win.infor$win.star[m] & deleterious.infor$position<=win.infor$win.end[m]]
    DeleBurdenMatrix=(w_Conserved * deleterious_window)
    DeleBurden=apply(DeleBurdenMatrix,2,function(x) sum(x,na.rm=T))
    #Homo_Dele=apply(deleterious_window,2,function(x) sum(x==2,na.rm=T))
    #Heter_Dele=apply(deleterious_window,2,function(x) sum(x==1,na.rm=T))
    #All_Dele=Homo_Dele+0.5*Heter_Dele
    All_Dele_Info=cbind(win.infor[m,],t(data.frame(DeleBurden)))
    #All_Dele_HomoInfo=cbind(win.infor[m,],t(data.frame(Homo_Dele)))
    #All_Dele_HeterInfo=cbind(win.infor[m,],t(data.frame(Heter_Dele)))
    All_Dele_Genome=rbind(All_Dele_Genome,All_Dele_Info)
   }
   }
   }

  S.can_index=match(S.can_ID,names(All_Dele_Genome))
  Inbred_index=match(Inbred_ID,names(All_Dele_Genome))
  Landrace_index=match(Landrace_ID,names(All_Dele_Genome))
  A626_index=match("A626",names(All_Dele_Genome))
  E463_index=match("E463",names(All_Dele_Genome))

All_Dele_Genome$Hotspot_AllDele_CNDMin=apply(All_Dele_Genome[,S.can_index],1,min)
All_Dele_Genome$Hotspot_AllDele_LandMin=apply(All_Dele_Genome[,Landrace_index],1,min)
All_Dele_Genome$Hotspot_AllDele_Min=apply(All_Dele_Genome[,c(Landrace_index,S.can_index)],1,min)

 sum(All_Dele_Genome$Hotspot_AllDele_CNDMin < All_Dele_Genome$Hotspot_AllDele_LandMin & All_Dele_Genome$Hotspot_AllDele_CNDMin  < All_Dele_Genome$A626)
 #503
sum(All_Dele_Genome$Hotspot_AllDele_LandMin < All_Dele_Genome$Hotspot_AllDele_CNDMin  & All_Dele_Genome$Hotspot_AllDele_LandMin < All_Dele_Genome$A626)
#91

All_Dele_Genome2=All_Dele_Genome[(All_Dele_Genome$A626 >0 | All_Dele_Genome$E463 >0),c(1:5,121,122,7,6,507:509,8:120,123:506)]

write.table(All_Dele_Genome2,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)
#840

  sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
#836
 All_Dele_Genome_CNDLessThanInbread=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
write.table(All_Dele_Genome_CNDLessThanInbread,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanInbread_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)



 All_Dele_Genome_CNDLessThanA626=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626),]
write.table(All_Dele_Genome_CNDLessThanA626,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanA626_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)

index_CND=match(S.can_ID,names(All_Dele_Genome2))
 All_Dele_Genome_CNDLessThanA626_CND= All_Dele_Genome_CNDLessThanA626[,c(1:12,index_CND)]
 write.table(All_Dele_Genome_CNDLessThanA626_CND,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanA626_RmSelectiveRegion_CND.deleInfo.bed",quote=F,sep="\t",row.name=F)
 CND_min=apply(All_Dele_Genome_CNDLessThanA626[,index_CND],1,min)
Reduced_A626=(All_Dele_Genome_CNDLessThanA626$A626 -All_Dele_Genome_CNDLessThanA626$Hotspot_AllDele_CNDMin )

Reduced_A626_100KSNP=(All_Dele_Genome_CNDLessThanA626$A626[All_Dele_Genome_CNDLessThanA626$SNP_num==100000] -All_Dele_Genome_CNDLessThanA626$Hotspot_AllDele_CNDMin[All_Dele_Genome_CNDLessThanA626$SNP_num==100000]  )


 CND_min_col=apply(All_Dele_Genome_CNDLessThanA626[,index_CND],1,which.min)
CND_min_col2=apply(All_Dele_Genome_CNDLessThanA626[All_Dele_Genome_CNDLessThanA626_CND$SNP_num!=100000,index_CND],1,which.min)

All_Dele_Genome_CNDLessThanA626


 All_Dele_Genome_CNDLessThanE463=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
write.table(All_Dele_Genome_CNDLessThanE463,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanE463_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)



index_cnd=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
sum(All_Dele_Genome2$A626[index_cnd]) #93974.43
sum(All_Dele_Genome2$A626[!index_cnd]) #95.03

sum(All_Dele_Genome2$A626) #94069.46  (non-selective region)        #105154.6 (all genome)
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd])+sum(All_Dele_Genome2$A626[!index_cnd])#40272.24


index_cnd_A626=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626)
sum(All_Dele_Genome2$A626[index_cnd_A626]) #93782.88
sum(All_Dele_Genome2$A626[!index_cnd_A626]) #286.58
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd_A626])+sum(All_Dele_Genome2$A626[!index_cnd_A626])
#40218.86

sum(All_Dele_Genome2$E463) ##93120.85 (non-selective region)                   #104544
index_cnd_E463=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
sum(All_Dele_Genome2$E463[index_cnd_E463])#92752.61
sum(All_Dele_Genome2$E463[!index_cnd_E463]) #368.24
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd_E463])+sum(All_Dele_Genome2$E463[!index_cnd_E463])
#40245.76

 All_Dele_Genome_CNDLessThanAll=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin) & (All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
#506
write.table(All_Dele_Genome_CNDLessThanAll,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanAll_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)

All_Dele_Genome_CNDLessThanAll_IncludeBothInbred=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin) & (All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 & All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
#501
write.table(All_Dele_Genome_CNDLessThanAll_IncludeBothInbred,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanAll_IncludeBothInbred_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)

###
  sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin & All_Dele_Genome2$Hotspot_AllDele_CNDMin  < All_Dele_Genome2$A626)
#503
 sum(All_Dele_Genome2$Hotspot_AllDele_LandMin < All_Dele_Genome2$Hotspot_AllDele_CNDMin  & All_Dele_Genome2$Hotspot_AllDele_LandMin < All_Dele_Genome2$A626)
#91
