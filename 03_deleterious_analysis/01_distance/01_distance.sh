#!/bin/bash
#SBATCH --partition=smp01
#SBATCH -N 1
#SBATCH -c 5

inputs=(
    "DP4_100.GQ10.Q30.MR0.5.maf0.0001.recode.vcf"
    "dom.Filter0.0001.vcf.gz"
    "NotDom.Filter0.0001.vcf.gz"
)

for input in "${inputs[@]}"; do
    output_name=$(echo ${input} | sed 's/\.vcf\.gz$//;s/\.vcf$//')

    ## make bed file
    plink --vcf ${input} \
          --make-bed --geno 0.2 --maf 0.01 \
          --out ${output_name}.maf0.01.geno0.2.filter

    ## genetic distance
    plink --bfile ${output_name}.maf0.01.geno0.2.filter \
          --distance-matrix \
          --out ${output_name}.maf0.01.geno0.2.filter.distance
done










