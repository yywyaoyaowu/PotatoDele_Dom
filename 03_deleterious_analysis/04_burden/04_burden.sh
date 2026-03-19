#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 3

for i in {01..12}
do
    chr=chr${i}
    vcf_file=${chr}.DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.vcf
    geno_prefix=${chr}.DP4_100.GQ10.Q30.MR0.5.maf0.0001
    GERP_file=/home/wuyaoyao/Solab/Forweizhi/00_Constrainted/Sol_msa_ConstraintedGERP2_withDepth_withMajorAllele_${chr}.txt
    subGroup_Can=Candolleanum_list
    subGroup_Lan=Landrace_list
    subGroup_Inbred=Inbred_list
    Rscript burdenstat.R ${chr} ${GERP_file} ${vcf_file} ${geno_prefix} $subGroup_Can $subGroup_Lan $subGroup_Inbred 2.75
done
