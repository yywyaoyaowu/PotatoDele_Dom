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

#ls ../perSite/*bed | awk -F '/' '{print $3}' | sed 's/.bed//g' > list
#grep Lan list > list_Lan
#grep Can list > list_Can

for i in $(cat list_Lan)
do
	cut -f1,3 ../perSite/${i}.bed > ${i}.pos
	nohup /public/software/env01/bin/vcftools --vcf Lan.snp.vcf --positions ${i}.pos --het --out ${i}.perSample &
done

for i in $(cat list_Can)
do
	cut -f1,3 ../perSite/${i}.bed > ${i}.pos
	nohup /public/software/env01/bin/vcftools --vcf Can.snp.vcf --positions ${i}.pos --het --out ${i}.perSample &
done



