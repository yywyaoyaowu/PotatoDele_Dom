#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Thu Oct  9 11:08:01 2025

# @File Name: work.sh
# @Description:
##################################################################

sed '1d' dat_same.txt | sort -k4 > dat_same.txt.sort
mv dat_same.txt.sort dat_same.txt
awk '{print $4}' dat_same.txt | uniq -c

source activate /public/agis/huangsanwen_group/huyong/software/anaconda3/envs/R4.1.1
Rscript generate_points.R Bolivia 41 | grep -v "Latitude" > points.txt
Rscript generate_points.R Colombia 32 | grep -v "Latitude" >> points.txt
Rscript generate_points.R Ecuador 3 | grep -v "Latitude" >> points.txt
Rscript generate_points.R Peru 95 | grep -v "Latitude" >> points.txt

paste dat_same.txt points.txt | awk '{print $1,$2,$3,$4,$5,$6}' OFS='\t' > dat_same_points.txt

grep -f final_data dat_same_points.txt > dat_same_points.txt.final
mv dat_same_points.txt.final dat_same_points.txt

