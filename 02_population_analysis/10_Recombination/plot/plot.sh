#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Sep 17 15:55:32 2025

# @File Name: work.sh
# @Description:
##################################################################

ln -s ../Lan_chr*output/*.PREDICT.txt . 
ln -s ../../Solanum_tuberosumDM.bed .
cp ~/dom_del/Recombination/plot/Solanum_tuberosumDM.length.txt .
ln -s /home/huyong/dom_del/09_diversity_final/Lan_Can.windowed.pi.dom .
ln -s /home/huyong/dom_del/09_diversity_final/Lan.windowed.pi .

for i in {01..12}
do
	grep -v "^#" /home/huyong/dom_del/99_results/Lan_chr${i}.snp.recode.vcf | cut -f1,2 | awk '{print $1"\t"$2-1"\t"$2}' >> Lan_snp_density.bed
done

cp Lan_chr01.snp.recode.PREDICT.txt Lan.snp.recode.PREDICT.txt
for i in {02..12}
do
        grep -v "chrom" Lan_chr${i}.snp.recode.PREDICT.txt >> Lan.snp.recode.PREDICT.txt
done
sed -i 's/b'\''//g' Lan.snp.recode.PREDICT.txt
sed -i 's/'\''//g' Lan.snp.recode.PREDICT.txt

source activate /public/agis/huangsanwen_group/huyong/software/anaconda3/envs/r-plot
Rscript recombination_line.R Lan.snp.recode.PREDICT.txt Solanum_tuberosumDM.length.txt
Rscript recombination_line_2.R Lan.snp.recode.PREDICT.txt Solanum_tuberosumDM.length.txt Lan 2.15
Rscript Recomb_Pi.R Lan.windowed.pi Lan.snp.recode.PREDICT.txt Landrace
Rscript Recomb_geneNum.R Lan.windowed.pi Lan.snp.recode.PREDICT.txt Solanum_tuberosumDM.bed Landrace
Rscript Recomb_PiRatio.R Lan_Can.windowed.pi.dom Lan.snp.recode.PREDICT.txt Landrace

DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_Conservation2_SolMsaNonMajorDele_subgroupDeleFreq_AllChrs_FisherPvalue_RecomTop5.Centremere5Mb.10Mb_rmbothmajor_header.txt

cut -f1-3 Lan.snp.recode.PREDICT.txt | grep -v chrom > df_windows.txt
/public/software/env01/bin/bedtools intersect -a df_windows.txt -b Lan_snp_density.bed -c > snp_Recomb_df_windows.csv
Rscript Recomb_del.R DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_Conservation2_SolMsaNonMajorDele_subgroupDeleFreq_AllChrs_FisherPvalue_RecomTop5.Centremere5Mb.10Mb_rmbothmajor_header.txt snp_Recomb_df_windows.csv Lan.snp.recode.PREDICT.txt

cut -f1-3 Lan_Can.windowed.pi.dom | grep -v CHROM > df_windows_dom.txt
/public/software/env01/bin/bedtools intersect -a df_windows_dom.txt -b Lan_snp_density.bed -c > snp_dom_df_windows.csv
Rscript del_PiRatio.R DMV6_DP4.100.GQ10.Q30.MR0.5.maf0.001_Conservation2_SolMsaNonMajorDele_subgroupDeleFreq_AllChrs_FisherPvalue_RecomTop5.Centremere5Mb.10Mb_rmbothmajor_header.txt Lan_Can.windowed.pi.dom snp_dom_df_windows.csv
/public/software/env01/bin/bedtools intersect -a Lan.snp.recode.PREDICT.txt -b Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.merge > Lan.Recomb.dom.txt


Rscript sort.R Lan.snp.recode.PREDICT.txt Lan.snp.recode.PREDICT.sort.txt
head -n 2383 Lan.snp.recode.PREDICT.sort.txt | sed '1d' | cut -f1-3 | sort -k1,1 -k2n,2 | /public/software/env01/bin/bedtools merge > Lan.Recombination.top5%
awk '{i+=($3-$2+1)}END{print i}' Lan.Recombination.top5%


head -n 7148 Lan.snp.recode.PREDICT.sort.txt | sed '1d' | cut -f1-3 | sort -k1,1 -k2n,2 | /public/software/env01/bin/bedtools merge > Lan.Recombination.top15percentage
awk '{i+=($3-$2+1)}END{print i}' Lan.Recombination.top15percentage

grep -v chr05 Lan.snp.recode.PREDICT.sort.txt | head -n6607 | sed '1d' | cut -f1-3 | sort -k1,1 -k2n,2 | /public/software/env01/bin/bedtools merge > Lan.Recombination.rmChr05.top15percentage


