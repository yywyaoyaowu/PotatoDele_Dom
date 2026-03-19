#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Jul 30 14:54:21 2025

# @File Name: work.sh
# @Description:
##################################################################

#for i in {01..12}
#do
#    echo """~/software/snpEff/scripts/snpEff ann Solanum_tuberosumDM chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf > chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf""" >> split_anno
#done
#/public/software/env01/bin/parallel -j 12 < split_anno

cut -f1-8 chr01.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf > DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info
for i in {02..12}
do
   cut -f1-8 chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf | grep -v "^#" >> DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info
done


paste <(grep -v "^#" DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info  | awk -F '\t' '{print $1"\t"$2-1"\t"$2}') <(grep -v "^#" DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info  | awk -F '|' '{print $2}') > DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info.bed


