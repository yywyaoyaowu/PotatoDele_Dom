#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Fri Jul 19 21:31:51 2024

# @File Name: work.sh
# @Description:
##################################################################

/home/huyong/software/PopLDdecay-3.42/bin/PopLDdecay -InVCF ../DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.vcf -MAF 0.05 -OutStat Candolleanum_LDdecay -SubPop Candolleanum_list


