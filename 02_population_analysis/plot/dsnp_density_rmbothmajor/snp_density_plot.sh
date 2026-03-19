#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Sun Nov  9 11:57:26 2025

# @File Name: work.sh
# @Description:
##################################################################

input=$1
prefix=$2

awk '$12 > 0 && $10 >= 2.75' ${input} | sed '1d' | cut -f1,3 | sed 's/chr/Can_chr/g' > Can_dsnp_density.csv
awk '$17 > 0 && $10 >= 2.75' ${input} | sed '1d' | cut -f1,3 | sed 's/chr/Lan_chr/g' > Lan_dsnp_density.csv

for i in {01..12}
do
    grep chr${i} Lan_dsnp_density.csv >> dsnp_density.csv
    grep chr${i} Can_dsnp_density.csv >> dsnp_density.csv
done
awk 'BEGIN{i=1}{print "dsnp_"i"\t"$0; i++}' dsnp_density.csv > multitrack_plot.csv

source activate R_4.1
Rscript snp_density_plot.R multitrack_plot.csv
mv Marker_Density.Trait1_Trait0.pdf ${prefix}.pdf

