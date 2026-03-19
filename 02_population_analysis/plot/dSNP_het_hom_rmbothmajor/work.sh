#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J vcftools_stat.sh
#SBATCH --error=err_%J_vcftools_stat.sh
#SBATCH --output=out_%J_vcftools_stat.sh
##################################################################
# @Author: huyong
# @Created Time : Mon Aug 18 09:15:13 2025

# @File Name: vcftools_stat.sh
# @Description:
##################################################################

input=DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_Conservation2_SolMsaNonMajorDele_subgroupDeleFreq_AllChrs_FisherPvalue_RecomTop5.Centremere5Mb.10Mb_rmbothmajor_header.txt

awk '$17 > 0 && $10 >= 2 && $10 < 2.75' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation2_2.75.bed
awk '$12 > 0 && $10 >= 2 && $10 < 2.75' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation2_2.75.bed

awk '$17 > 0 && $10 >= 2.75 && $10 < 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation2.75_3.5.bed
awk '$12 > 0 && $10 >= 2.75 && $10 < 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation2.75_3.5.bed

awk '$17 > 0 && $10 >= 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation3.5.bed
awk '$12 > 0 && $10 >= 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation3.5.bed

Rscript hetDelRatio.R Lan_del_Conservation2_2.75.bed Lan_del_Conservation2.75_3.5.bed Lan_del_Conservation3.5.bed Can_del_Conservation2_2.75.bed Can_del_Conservation2.75_3.5.bed Can_del_Conservation3.5.bed dSNP_het_ratio_rmbothmajor.pdf



####### 
awk '$17 > 0 && $10 >= 2' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation2.bed
awk '$12 > 0 && $10 >= 2' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation2.bed

awk '$17 > 0 && $10 >= 2.75' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation2.75.bed
awk '$12 > 0 && $10 >= 2.75' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation2.75.bed

awk '$17 > 0 && $10 >= 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$19/($19+$18)}' > Lan_del_Conservation3.5.bed
awk '$12 > 0 && $10 >= 3.5' ${input} | grep -v "chrom" | awk '{print $1"\t"$3"\t"$14/($14+$13)}' > Can_del_Conservation3.5.bed

Rscript hetDelRatio.R Lan_del_Conservation2.bed Lan_del_Conservation2.75.bed Lan_del_Conservation3.5.bed Can_del_Conservation2.bed Can_del_Conservation2.75.bed Can_del_Conservation3.5.bed dSNP_het_ratio_rmbothmajor_2.pdf



