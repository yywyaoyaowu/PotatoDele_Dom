library(SNPRelate)
library(regress)
args<-commandArgs(TRUE)
chr=args[1]
GERP_file=args[2]
vcf_file=args[3]
geno_prefix=args[4]
pop1_file=args[5]
pop2_file=args[6]
pop3_file=args[7]
ConservationCutoff=as.numeric(args[8])

vcf_freq_file=paste0(geno_prefix,"_Freq")
gbs_file=paste(geno_prefix,"_RefCount.gds",sep='')

print(geno_prefix)
print(GERP_file)
snpgdsVCF2GDS(vcf_file, gbs_file, method="biallelic.only")
geno <- snpgdsOpen(gbs_file) ##if gbs from vcf, ref allele count:##Alt homo : 0; heter 1,Ref homo:2 ; missing NA
snp_infor_combine=snpgdsSNPList(geno, sample.id=NULL)
snp_infor_combine[1:3,]
sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
sample.id[1:3]
pop1=read.table(pop1_file)
pop2=read.table(pop2_file)
pop3=read.table(pop3_file)
print (paste("pop1 size is ", nrow(pop1)))
print (paste("pop2 size is ", nrow(pop2)))
print (paste("pop3 size is ", nrow(pop3)))

X_all_combine <- snpgdsGetGeno(geno)
X_all_combine[1:5,1:5]
##get the burden weight
GERP_OneChr=read.table(GERP_file,header=T)
names(GERP_OneChr)[1:6]=c("Chrom","Pos_start","Pos","GERP","Neutral","Depth")
GERP_FilterConserved=subset(GERP_OneChr,GERP>=ConservationCutoff)
# heter, home 2,deleter1, weight by GERP_Z,calculated by depth >=20, conserved GERP score:
# get the constrainted snp: dele snp matrix and gerp score #######################
index_Conserved=match(snp_infor_combine$position,GERP_FilterConserved$Pos)
snp_infor_Constrainted=snp_infor_combine[!is.na(index_Conserved),]
snp_infor_Constrainted[1:3,]
X_Conserved=X_all_combine[,!is.na(index_Conserved)]
w_Conserved=GERP_FilterConserved$GERP[match(snp_infor_Constrainted$position,GERP_FilterConserved$Pos)]

SnpNum_Combine=nrow(snp_infor_combine)
Conserved_SnpNnum=sum(!is.na(index_Conserved))

Geno_0_LineNum=apply(X_Conserved,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_Conserved,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum))
print(sum(Geno_2_LineNum))
print(sum(Geno_0_LineNum<Geno_2_LineNum))

vcf_freq=read.table(vcf_freq_file,skip=1)
names(vcf_freq)[1:8]=c("Chrom","position","N_ALLELES","N_ch","DM_Allele","RefFreq","Alt_Allele","AltFreq")
vcf_freq[1:3,]

snp_infor_Constrainted$RefAllele=vcf_freq$DM_Allele[match(snp_infor_Constrainted$position,vcf_freq$position)]
snp_infor_Constrainted$AltAllele=vcf_freq$Alt_Allele[match(snp_infor_Constrainted$position,vcf_freq$position)]

snp_infor_Constrainted$Sol100Major=GERP_FilterConserved$Sol100_MajorAllel[match(snp_infor_Constrainted$position,GERP_FilterConserved$Pos)]
sum(is.na(snp_infor_Constrainted$Sol100Major))
snp_infor_Constrainted$GERP=w_Conserved
X_DeleNum=(X_Conserved-2)*(-1)##X_DeleNum, Ref is 0; Alt is 2

index_concis=which(snp_infor_Constrainted$Sol100Major==snp_infor_Constrainted$RefAllele)  #ref is ancestor allele, that is neutral/benefit allele
index_turnOver=which(snp_infor_Constrainted$Sol100Major==snp_infor_Constrainted$AltAllele)  #ref is derived allele, that is dele
index_bothDele=which(snp_infor_Constrainted$Sol100Major!=snp_infor_Constrainted$RefAllele & snp_infor_Constrainted$Sol100Major!=snp_infor_Constrainted$AltAllele)
#both ref and alt is derived
snp_infor_Constrainted$DeleAllele=NA
snp_infor_Constrainted$DeleAllele[index_concis]="AltAlele"
snp_infor_Constrainted$DeleAllele[index_turnOver]="RefAlele"
snp_infor_Constrainted$DeleAllele[index_bothDele]="BothRefAlt"

X_DeleNum2=X_DeleNum
X_DeleNum2[,index_turnOver]=(X_DeleNum2[,index_turnOver]-2)*(-1)

X_DeleNum3=X_DeleNum2
X_DeleNum3[,index_bothDele]=(X_DeleNum2[,index_bothDele]*0-2)*(-1)

justDeleGeno=data.frame(t(X_DeleNum3))
names(justDeleGeno)=sample.id
#Dele_Genotype=cbind(snp_infor_Constrainted,justDeleGeno)
#write.table(Dele_Genotype,paste(geno_prefix,"_ConservationCutoff",ConservationCutoff,"_SolMsaNonMajorDele_Genotype.txt",sep=""),quote=F,row.name=F,sep='\t')

###STAT dele info for each accession
HomeDeleNum=apply(X_DeleNum3,1,function(x) sum(x==2,na.rm=T))
HeterDeleNum=apply(X_DeleNum3,1,function(x) sum(x==1,na.rm=T))

Dele_matrix_weight <- t(w_Conserved * t(X_DeleNum3))
Burden_index_EachiLine=apply(Dele_matrix_weight,1,function(x) sum(x,na.rm=T))

X_Dele_Heter=X_DeleNum3
X_Dele_Heter[X_DeleNum3==2]=0

Burden_index_HeterMatrix=t(w_Conserved * t(X_Dele_Heter))
Burden_index_Heter=apply(Burden_index_HeterMatrix,1,function(x) sum(x,na.rm=T))
Burden_index_Homo=(Burden_index_EachiLine-Burden_index_Heter)/2

BurdenInfor=data.frame(chr,SnpNum_Combine,Conserved_SnpNnum,sample.id,HeterDeleNum,HomeDeleNum,Burden_index_EachiLine,Burden_index_Heter,Burden_index_Homo)
write.table(BurdenInfor,paste(geno_prefix,"_SolMsaNonMajorDele_BurdenInfo_GERP2.75.txt",sep=""),quote=F,row.name=F,sep='\t')

##sub group stat
index_pop1=NA
Can_deleGenotype=NA
index_pop1=!is.na(match(sample.id,pop1[,1]))
Can_deleGenotype=justDeleGeno[,index_pop1]

Can_DeleHome=apply(Can_deleGenotype,1,function(x) sum(x==2,na.rm=T))
Can_DeleHeter=apply(Can_deleGenotype,1,function(x) sum(x==1,na.rm=T))
Can_NeuHomo=apply(Can_deleGenotype,1,function(x) sum(x==0,na.rm=T))
Can_Geno=apply(Can_deleGenotype,1,function(x) sum(!is.na(x),na.rm=T))
Can_DeleFreq=(Can_DeleHome*2+Can_DeleHeter)/(2*Can_Geno)

index_pop2=NA
Lan_deleGenotype=NA
index_pop2=!is.na(match(sample.id,pop2[,1]))
Lan_deleGenotype=justDeleGeno[,index_pop2]

Lan_DeleHome=apply(Lan_deleGenotype,1,function(x) sum(x==2,na.rm=T))
Lan_DeleHeter=apply(Lan_deleGenotype,1,function(x) sum(x==1,na.rm=T))
Lan_NeuHomo=apply(Lan_deleGenotype,1,function(x) sum(x==0,na.rm=T))
Lan_Geno=apply(Lan_deleGenotype,1,function(x) sum(!is.na(x),na.rm=T))
Lan_DeleFreq=(Lan_DeleHome*2+Lan_DeleHeter)/(2*Lan_Geno)

index_pop3=NA
Inbred_deleGenotype=NA
index_pop3=!is.na(match(sample.id,pop3[,1]))
Inbred_deleGenotype=justDeleGeno[,index_pop3]

subgroupDeleFreq=data.frame(snp_infor_Constrainted,Can_DeleFreq,Can_DeleHome,Can_DeleHeter,Can_NeuHomo,Can_Geno,Lan_DeleFreq,Lan_DeleHome,Lan_DeleHeter,Lan_NeuHomo,Lan_Geno)
write.table(subgroupDeleFreq,paste(geno_prefix,"_SolMsaNonMajorDele_subgroupDeleFreq_GERP2.75.txt",sep=""),quote=F,row.name=F,sep='\t')
Dele_Genotype.f=cbind(subgroupDeleFreq,Can_deleGenotype,Lan_deleGenotype,Inbred_deleGenotype)
write.table(Dele_Genotype.f,paste(geno_prefix,"_SolMsaNonMajorDele_Genotype_GERP2.75.txt",sep=""),quote=F,row.name=F,sep='\t')
