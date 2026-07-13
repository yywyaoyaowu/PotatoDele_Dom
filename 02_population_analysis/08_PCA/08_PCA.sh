#!/bin/bash
#SBATCH --partition=big
#SBATCH -N 1
#SBATCH -c 15
##################################################################
# @Author: huyong
# @Created Time : Fri Jan 24 15:34:05 2025

# @File Name: t.sh
# @Description:
##################################################################


/public/lihongbo/bin/plink --threads 40 --noweb --vcf DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.vcf --maf 0.05 --make-bed --out DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile
/public/lihongbo/bin/plink --threads 40 --bfile DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile --pca 502 --out DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile.pca502

cut -d" " -f2 DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile.pca502.eigenvec > f2
cat sample_info | while read a b other; do sed -i 's/'${a}'/'${b}'/g' f2; done
paste <(cut -d" " -f1 DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile.pca502.eigenvec) <(cat f2) <(cut -d" " -f3- DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.snp.pca.bfile.pca502.eigenvec) > eigenvec_pca502
sed -i 's/ /\t/g' eigenvec_pca502
paste <(echo -e "name\tgroup") <(for i in {1..502}; do echo PC${i}; done | xargs -n10000 | sed 's/ /\t/g') > eigenvec_header_pca502
cat eigenvec_header_pca502 eigenvec_pca502 > e
mv e eigenvec_pca502




