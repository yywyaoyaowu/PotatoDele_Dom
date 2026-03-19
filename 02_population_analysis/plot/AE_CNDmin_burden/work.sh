#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Dec 17 13:13:51 2025

# @File Name: work.sh
# @Description:
##################################################################

cut -f1,2,3,4,7,9 DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.bed > DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.plot.bed
Rscript plot.R DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.plot.bed DMRef_1MbSNP.window_AllDeleBurden_AtLeastOne.pdf

cut -f1,2,3,4 DMRef_1MbSNP.window_HomoDeleBurden_AtLeastOne_F1HomoDele.bed > F1HomoDele.plot.bed
sed -i 's/A626/F1HomoDele/g' F1HomoDele.plot.bed



Rscript plot.R DMRef_1MbSNP.window_HomoDeleBurden_AtLeastOne_F1HomoDele.plot.bed DMRef_1MbSNP.window_HomoDeleBurden_AtLeastOne_F1HomoDele.pdf


