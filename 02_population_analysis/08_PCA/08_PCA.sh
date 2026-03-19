#!/bin/bash
#SBATCH --partition=big
#SBATCH -N 1
#SBATCH -c 15
#SBATCH -J t.sh
#SBATCH --error=err_%J_t.sh
#SBATCH --output=out_%J_t.sh
##################################################################
# @Author: huyong
# @Created Time : Fri Jan 24 15:34:05 2025

# @File Name: t.sh
# @Description:
##################################################################

/public/lihongbo/bin/plink --threads 15 --noweb --vcf DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf --maf 0.05 --make-bed --out DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile
/public/lihongbo/bin/plink --threads 15 --bfile DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile --pca 338 --out DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile.pca338

cut -d" " -f2 DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile.pca338.eigenvec > f2
cat sample_info | while read a b other; do sed -i 's/'${a}'/'${b}'/g' f2; done
paste <(cut -d" " -f1 DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile.pca338.eigenvec) <(cat f2) <(cut -d" " -f3- DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.snp.pca.bfile.pca338.eigenvec) > eigenvec_pca338
sed -i 's/ /\t/g' eigenvec_pca338
paste <(echo -e "name\tgroup") <(for i in {1..338}; do echo PC${i}; done | xargs -n10000 | sed 's/ /\t/g') > eigenvec_header_pca338
cat eigenvec_header_pca338 eigenvec_pca338 > e
mv e eigenvec_pca338




