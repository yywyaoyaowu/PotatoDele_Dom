#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Mon Jan 20 17:39:39 2025

# @File Name: work.sh
# @Description:
##################################################################

source activate R_topGO
Rscript GOenrichment.R DM_1-3_516_R44_potato.v6.1.working_models.iprscan_go_terms.txt gene_in_dom_top0.05_merge_geneID

conda deactivate
source activate R_4.1
Rscript GOenrichment_plot.R gene_in_dom_top0_GoInfo_GoEnrichment.txt

