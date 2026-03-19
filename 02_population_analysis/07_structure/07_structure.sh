#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J admixture.sh
#SBATCH --error=err_%J_admixture.sh
#SBATCH --output=out_%J_admixture.sh
##################################################################
# @Author: huyong
# @Created Time : Fri May 17 22:47:00 2024

# @File Name: admixture.sh
# @Description:
##################################################################

ls `pwd`/*.vcf > vcf.list
gatk MergeVcfs -I vcf.list -O DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf
perl /home/huyong/software/VCF_add_id-master/VCF_add_id.pl DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf DP4_100.GQ10.Q30.MR0.5.maf0.001.addID.vcf

## LD pruning and make bed file
/home/lihongbo/bin/plink --vcf DP4_100.GQ10.Q30.MR0.5.maf0.001.addID.vcf --maf 0.05 --indep-pairwise 50 10 0.2 --double-id -out DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD --make-bed --allow-extra-chr
/home/lihongbo/bin/plink --bfile DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD --extract DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.prune.in --out DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture --make-bed --allow-extra-chr

for i in {2..15}
do
	echo """#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 1

/public/software/env01/bin/admixture -j10 --cv DP4_100.GQ10.Q30.MR0.5.maf0.001.addID_LD.admixture.bed ${i} | tee log${i}.out
""" > ${i}.sh
sbatch ${i}.sh
done



