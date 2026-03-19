#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -J work.sh
#SBATCH --qos=queue1
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Sat Aug 27 23:51:37 2022

# @File Name: work.sh
# @Description:
##################################################################

for i in $(cat list)
do
	mkdir ${i}
	raw=($(ls ../00_raw_data/${i}/*gz))
	echo """#!/bin/bash
#SBATCH --partition=low,big
#SBATCH -N 1
#SBATCH -c 4
source activate annotation_2
fastp --detect_adapter_for_pe -w 4 -i ../${raw[0]} -I ../${raw[1]} -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz --json ${i}_fastp.json --html ${i}_fastp.html""" >> ${i}/${i}_fastp.sh

cd ${i}
sbatch ${i}_fastp.sh
cd ..
done

