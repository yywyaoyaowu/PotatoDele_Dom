#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Nov  5 15:19:19 2025

# @File Name: work.sh
# @Description:
##################################################################

bash 1.sh
bash 2.sh
bash 3.sh

paste <(ls `pwd`/*gz) <(ls `pwd`/*gz | awk -F '/' '{print $NF}' | sed 's/_LDdecay.stat.gz//g') > LDdecayResult.list
perl /home/huyong/software/PopLDdecay-3.42/bin/Plot_MultiPop.pl -inList LDdecayResult.list -output LDdecayResult


