#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Fri Dec 26 17:40:18 2025

# @File Name: work.sh
# @Description:
##################################################################

ln -s ~/tmp/TE/Solanum_tuberosumDM.mod.EDTA.TEanno.gff3 .
ln -s ../Solanum_tuberosumDM.gene.bed .
#cat /home/huyong/dom_del/circro/DM_v6.1_karyotype.txt | cut -f4 | xargs -n1000 | sed 's/ /,/g'

#### gene
/public/software/env01/bin/bedtools makewindows -g ../Solanum_tuberosumDM.len -w 500000 > Solanum_tuberosumDM.win500k
/public/software/env01/bin/bedtools intersect -b Solanum_tuberosumDM.gene.bed -a Solanum_tuberosumDM.win500k -c > Solanum_tuberosumDM.gene_density.win500k

#### TE
grep -v "#" Solanum_tuberosumDM.mod.EDTA.TEanno.gff3 | cut -f1,4-5 > TE
/public/software/env01/bin/bedtools merge -i TE > TE_cor
/public/software/env01/bin/bedtools intersect -b TE_cor -a Solanum_tuberosumDM.win500k -c > Solanum_tuberosumDM_chr_TE_density_500k


awk '{print $2"\t"$3-1"\t"$3}' ../Results_stat/genome_snp_density/snp_density.csv > snp_density
/public/software/env01/bin/bedtools intersect -b snp_density -a Solanum_tuberosumDM.win500k -c > Solanum_tuberosumDM_chr_snp_density_500k

/public/software/env01/bin/bedtools intersect -b maf0.05.bed -a Solanum_tuberosumDM.win500k -c > Solanum_tuberosumDM_chr_snp_density_maf0.05_500k

#### dele
awk '$11 > 0 && $9 >= 2.75' Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_SolMsaNonMajorDele_subgroupDeleFreq_GERP2.75.txt.filter.rmBothMajor.txt | cut -f1-3 > dele
/public/software/env01/bin/bedtools intersect -b dele -a Solanum_tuberosumDM.win500k -c > Solanum_tuberosumDM_dele_density_500k

~/software/circos-0.69-9/bin/circos -conf ./circos.conf -outputfile DM_dele -outputdir ./






