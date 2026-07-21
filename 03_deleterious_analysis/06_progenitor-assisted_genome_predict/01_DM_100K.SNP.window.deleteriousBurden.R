chr.len=read.table("/home/wuyaoyao/03-Solanaceae/09_SolOmics/01_Assembly/01_genome/Solanum_tuberosumDM.fa.fai",sep="\t")
#GroupInfo=read.table("/home/wuyaoyao/04_Potato/04_DomDele/01_VCF/sample_info_noEtu_final")
GroupInfo=read.table("/home/wuyaoyao/04_Potato/04_DomDele/DMref_502Accessions/01_VCF/sample_info_502")

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
   VCF_Freq=read.table(paste0("DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_",chrs[i],"_Freq"),skip=1)
   WinNum=floor(nrow(VCF_Freq)/win.len)
   #chr_end=floor((chr.len[i,2]-win.len)/step)+1
   win.star_XL=seq(0,WinNum)*win.len+1
   win.end_XL=win.star_XL+win.len-1
   win.star=VCF_Freq[win.star_XL,2]
   win.end=VCF_Freq[win.end_XL,2]
  chr=rep(chrs[i],length(win.star))
  win.infor=data.frame(chr,win.star,win.end)
  win.infor[nrow(win.infor),3]=chr.len[i,2]
  deleterious.infor_ALL=read.table(paste0("DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_",chrs[i],"_ConservationCutoff2.75_SolMsaNonMajorDele_Genotype.txt"),header=T)
  deleterious.infor=deleterious.infor_ALL[!(deleterious.infor_ALL$Land_DeleFreq > 0.5 & deleterious.infor_ALL$SCand_DeleFreq > 0.5) & !is.na(deleterious.infor_ALL$Land_DeleFreq) & !is.na(deleterious.infor_ALL$SCand_DeleFreq) & (deleterious.infor_ALL$Land_DeleFreq >0 | deleterious.infor_ALL$SCand_DeleFreq >0),]
 dim(deleterious.infor)# 1261004      24

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

  S.can_index=match(S.can_ID,names(All_Dele_Genome))
  Inbred_index=match(Inbred_ID,names(All_Dele_Genome))
  Landrace_index=match(Landrace_ID,names(All_Dele_Genome))
  A626_index=match("A626",names(All_Dele_Genome))
  E463_index=match("E463",names(All_Dele_Genome))

All_Dele_Genome$Hotspot_AllDele_CNDMin=apply(All_Dele_Genome[,S.can_index],1,min)
All_Dele_Genome$Hotspot_AllDele_LandMin=apply(All_Dele_Genome[,Landrace_index],1,min)
All_Dele_Genome$Hotspot_AllDele_Min=apply(All_Dele_Genome[,c(Landrace_index,S.can_index)],1,min)

 sum(All_Dele_Genome$Hotspot_AllDele_CNDMin < All_Dele_Genome$Hotspot_AllDele_LandMin & All_Dele_Genome$Hotspot_AllDele_CNDMin  < All_Dele_Genome$A626)
 #500
sum(All_Dele_Genome$Hotspot_AllDele_LandMin < All_Dele_Genome$Hotspot_AllDele_CNDMin  & All_Dele_Genome$Hotspot_AllDele_LandMin < All_Dele_Genome$A626)
#94

All_Dele_Genome2=All_Dele_Genome[(All_Dele_Genome$A626 >0 | All_Dele_Genome$E463 >0 | All_Dele_Genome$Hotspot_AllDele_Min >0 ),c(1:4,120,121,6,5,506:508,7:119,122:505)]

write.table(All_Dele_Genome2,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne.bed",quote=F,sep="\t",row.name=F)
#816

  sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
#813
 All_Dele_Genome_CNDLessThanInbread=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
write.table(All_Dele_Genome_CNDLessThanInbread,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanInbread.bed",quote=F,sep="\t",row.name=F)

index_cnd=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
sum(All_Dele_Genome2$A626[index_cnd]) #104889.9
sum(All_Dele_Genome2$A626[!index_cnd]) #264.65

sum(All_Dele_Genome2$A626) #105154.6
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd])+sum(All_Dele_Genome2$A626[!index_cnd])# 46471.29


index_cnd_A626=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626)
sum(All_Dele_Genome2$A626[index_cnd_A626]) #104705.8
sum(All_Dele_Genome2$A626[!index_cnd_A626]) #448.76
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd_A626])+sum(All_Dele_Genome2$A626[!index_cnd_A626])
#46436.22

sum(All_Dele_Genome2$E463)#104544
index_cnd_E463=(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463)
sum(All_Dele_Genome2$E463[index_cnd_E463])#103934.6
sum(All_Dele_Genome2$E463[!index_cnd_E463]) #609.45
sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin[index_cnd_E463])+sum(All_Dele_Genome2$E463[!index_cnd_E463])
#46471.75

 All_Dele_Genome_CNDLessThanAll=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin) & (All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 | All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
#500
write.table(All_Dele_Genome_CNDLessThanAll,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanAll.bed",quote=F,sep="\t",row.name=F)

All_Dele_Genome_CNDLessThanAll_IncludeBothInbred=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin) & (All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626 & All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
#499
write.table(All_Dele_Genome_CNDLessThanAll_IncludeBothInbred,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanAll_IncludeBothInbred.bed",quote=F,sep="\t",row.name=F)

###
  sum(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$Hotspot_AllDele_LandMin & All_Dele_Genome2$Hotspot_AllDele_CNDMin  < All_Dele_Genome2$A626)
#500
 sum(All_Dele_Genome2$Hotspot_AllDele_LandMin < All_Dele_Genome2$Hotspot_AllDele_CNDMin  & All_Dele_Genome2$Hotspot_AllDele_LandMin < All_Dele_Genome2$A626)
#94












    win.infor$Hotspot_AllDele_S.canMedian[m]=median(All_Dele[S.can_index],na.rm=T)
    win.infor$Hotspot_AllDele_LandraceMedian[m]=median(All_Dele[Landrace_index],na.rm=T)

    win.infor$Hotspot_AllDele_S.canMean[m]=mean(All_Dele[S.can_index],na.rm=T)
    win.infor$Hotspot_AllDele_LandraceMean[m]=mean(All_Dele[Landrace_index],na.rm=T)

    win.infor$Hotspot_AllDele_S.canMin[m]=min(All_Dele[S.can_index],na.rm=T)
    win.infor$Hotspot_AllDele_LandraceMin[m]=min(All_Dele[Landrace_index],na.rm=T)
    win.infor$Hotspot_AllDele_Min[m]=min(All_Dele,na.rm=T)
    win.infor$Hotspot_AllDele_p[m]=t.test(All_Dele[S.can_index],All_Dele[Landrace_index])$p.value

      if ( All_Dele[A626_index] > 0){
          win.infor$CNDlessNum_AllDele[m]=sum((All_Dele[S.can_index] / All_Dele[A626_index]) <= 0.8)
          win.infor$LandracelessNum_AllDele[m]=sum((All_Dele[Landrace_index] / All_Dele[A626_index]) <= 0.8)
      }
     # win.infor$Hotspot_AllDele_fold[m]= win.infor$Hotspot_AllDele_p[m]=NA


    All_Dele_Genome=rbind(All_Dele_Genome,All_Dele_Info)
  }
         #All.win.infor=rbind(All.win.infor,win.infor)
}

All_Dele_Genome$Hotspot_AllDele_fold_median=All_Dele_Genome$Hotspot_AllDele_S.canMedian / All_Dele_Genome$Hotspot_AllDele_LandraceMedian
All_Dele_Genome$Hotspot_AllDele_fold_mean=All_Dele_Genome$Hotspot_AllDele_S.canMean / All_Dele_Genome$Hotspot_AllDele_LandraceMean

names(All_Dele_Genome[,8:9])
All_Dele_Genome2=All_Dele_Genome[,c(1:13,16,351:353,14,15,17:350)]
sum(All_Dele_Genome$Hotspot_AllDele_S.canMin < All_Dele_Genome$A626 | All_Dele_Genome$Hotspot_AllDele_LandraceMin <  All_Dele_Genome$A626 )
write.table(All_Dele_Genome,"DMRef_1Mb.window.200k.step_AllDeleNum_S.canLandrace.txt",quote=F,sep="\t",row.name=F)
write.table(All_Dele_Genome2,"DMRef_1Mb.window.200k.step_AllDeleNum_S.canLandrace.bed",quote=F,sep="\t",row.name=F)


ThreshouldTop5_A626= sort(All_Dele_Genome$A626,decreasing=T)[nrow(All_Dele_Genome)*0.05]
 ThreshouldA626=ThreshouldTop5_A626

 Landrace_SigMoreDele=All_Dele_Genome[(All_Dele_Genome$Hotspot_AllDele_fold_mean < 1 & All_Dele_Genome$Hotspot_AllDele_p < 0.01 ),]
 Landrace_SigLessDele=All_Dele_Genome[(All_Dele_Genome$Hotspot_AllDele_fold_mean > 1 & All_Dele_Genome$Hotspot_AllDele_p < 0.01 ),]

Landrace_SigMoreDele_InbredTop=All_Dele_Genome[(All_Dele_Genome$Hotspot_AllDele_fold_mean < 1 & All_Dele_Genome$Hotspot_AllDele_p < 0.01 & All_Dele_Genome$A626 >= ThreshouldA626),]
Landrace_SigLessDele_InbredTop=All_Dele_Genome[(All_Dele_Genome$Hotspot_AllDele_fold_mean > 1 & All_Dele_Genome$Hotspot_AllDele_p < 0.01 & All_Dele_Genome$A626 >= ThreshouldA626),]


sum(All_Dele_Genome$LandracelessNum_AllDele < All_Dele_Genome$CNDlessNum_AllDele & All_Dele_Genome$A626 >= 98)
 A626_DeleTopRegion=subset(All_Dele_Genome,A626 >= ThreshouldA626)

# sum (All_Dele_Genome$Hotspot_AllDele_AccumulatedPer < 0 & All_Dele_Genome$Hotspot_AllDele_p < 0.01,na.rm=T)
[1] 2074
# sum (All_Dele_Genome$Hotspot_AllDele_AccumulatedPer > 0 & All_Dele_Genome$Hotspot_AllDele_p < 0.01,na.rm=T)
[1] 363


ThreshouldTop5_A626= sort(All_Dele_Genome$A626,decreasing=T)[nrow(All_Dele_Genome)*0.05]
A626_DeleTopRegion=subset(All_Dele_Genome,A626 >= ThreshouldA626)
write.table(A626_DeleTopRegion,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTopRegion.bed",quote=F,sep="\t",row.name=F)
AtLeastOneCND_Less=All_Dele_Genome[(All_Dele_Genome$CNDlessNum_AllDele > 0 & All_Dele_Genome$A626 >= ThreshouldTop5_A626),]
AtLeastOneCND.More_Less=All_Dele_Genome[(All_Dele_Genome$CNDlessNum_AllDele > All_Dele_Genome$LandracelessNum_AllDele & All_Dele_Genome$A626 >= ThreshouldTop5_A626),]
#174
write.table(AtLeastOneCND_Less,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTopRegion_AtLeastOneCND_Less.bed",quote=F,sep="\t",row.name=F)
write.table(AtLeastOneCND.More_Less,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTopRegion_AtLeastOneCND.MoreLess.bed",quote=F,sep="\t",row.name=F)
#53





ThreshouldTop5_A626= sort(All_Dele_Genome$A626,decreasing=T)[nrow(All_Dele_Genome)*0.01]
A626_DeleTopRegion=subset(All_Dele_Genome,A626 >= ThreshouldA626)
write.table(A626_DeleTopRegion,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTopRegion.bed",quote=F,sep="\t",row.name=F)
AtLeastOneCND_Less=All_Dele_Genome[(All_Dele_Genome$CNDlessNum_AllDele > 0 & All_Dele_Genome$A626 >= ThreshouldTop5_A626),]
AtLeastOneCND.More_Less=All_Dele_Genome[(All_Dele_Genome$CNDlessNum_AllDele > All_Dele_Genome$LandracelessNum_AllDele & All_Dele_Genome$A626 >= ThreshouldTop5_A626),]
nrow(AtLeastOneCND_Less) #38
nrow(AtLeastOneCND.More_Less)#11
write.table(AtLeastOneCND_Less,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTop1Region_AtLeastOneCND_Less.bed",quote=F,sep="\t",row.name=F)
write.table(AtLeastOneCND.More_Less,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTop1Region_AtLeastOneCND.MoreLess.bed",quote=F,sep="\t",row.name=F)









write.table(A626_DeleTopRegion,"DMRef_1Mb.window.200k.step_AllDeleNum_A626_DeleTopRegion.bed",quote=F,sep="\t",row.name=F)

####sliding window merge
#All_Dele_Genome=read.table("graphRef_1Mb.window.200k.step_AllDeleNum_S.canLandrace.txt",header=T)
All_Dele_Genome_WithDele=All_Dele_Genome[All_Dele_Genome$Hotspot_AllDele_S.canMean > 0 | All_Dele_Genome$Hotspot_AllDele_LandraceMean >0,]
## accumulation
Landrace_SigMoreDele=All_Dele_Genome_WithDele[(All_Dele_Genome_WithDele$Hotspot_AllDele_fold_mean < 1 & All_Dele_Genome_WithDele$Hotspot_AllDele_p < 0.01),1:17]
#362 (p 0.01)
Landrace_LessDele=All_Dele_Genome_WithDele[(All_Dele_Genome_WithDele$Hotspot_AllDele_fold_mean > 1 & All_Dele_Genome_WithDele$Hotspot_AllDele_p < 0.01),1:17]
#2072

Landrace_SigMoreDeleFC15=All_Dele_Genome[All_Dele_Genome$Hotspot_AllDele_fold <= 0.6666667 & !is.na(All_Dele_Genome$Hotspot_AllDele_fold) & All_Dele_Genome$Hotspot_AllDele_p < 0.01,1:17]
#120 (p 0.05) ; 112 (p 0.01)

Landrace_SigMoreDele[,2]=Landrace_SigMoreDele[,2]-1
Landrace_SigMoreDeleFC15[,2]=Landrace_SigMoreDeleFC15[,2]-1

Landrace_SigLessDeleFC=All_Dele_Genome[All_Dele_Genome$Hotspot_AllDele_fold > 1 & !is.na(All_Dele_Genome$Hotspot_AllDele_fold) & All_Dele_Genome$Hotspot_AllDele_p < 0.01,1:17]
#1618(p 0.05) ; 1496 (p 0.01)
Landrace_SigLessDeleFC15=All_Dele_Genome[All_Dele_Genome$Hotspot_AllDele_fold >= 1.5 & !is.na(All_Dele_Genome$Hotspot_AllDele_fold) & All_Dele_Genome$Hotspot_AllDele_p < 0.01,1:17]
#174 (p 0.05); 165 (p 0.01)
Landrace_SigLessDeleFC[,2]=Landrace_SigLessDeleFC[,2]-1
Landrace_SigLessDeleFC15[,2]=Landrace_SigLessDeleFC15[,2]-1


write.table(Landrace_SigMoreDele,"DMRef_1Mb.window.200k.step_AllDeleNum_LandraceSigAccumDele.bed",quote=F,sep="\t",row.name=F)
write.table(Landrace_SigMoreDeleFC15,"graphRef_1Mb.window.200k.step_AllDeleNum_LandraceSigAccumDele15FC.bed",quote=F,sep="\t",row.name=F)

write.table(Landrace_SigLessDeleFC,"DMRef_1Mb.window.200k.step_AllDeleNum_LandraceSigPurgeDele.bed",quote=F,sep="\t",row.name=F)
write.table(Landrace_SigLessDeleFC15,"DMRef_1Mb.window.200k.step_AllDeleNum_LandraceSigPurgeDele15FC.bed",quote=F,sep="\t",row.name=F)
