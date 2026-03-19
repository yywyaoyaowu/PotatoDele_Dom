#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 5

##calculate allele frequency in different panels
panels=(
    "all:" 
    "can:Candolleanum_list" 
    "lan:Landrace_list" 
    "precan:CandolleanumPre_list"
)

for i in {01..12}; do
    for panel in "${panels[@]}"; do
        panel_name=${panel%:*}
        list_file=${panel#*:}
        
        input_pre="chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.0001"
        
        if [ "$panel_name" = "all" ]; then
            /public/software/env01/bin/vcftools --vcf ${input_pre}.recode.vcf --freq --out ${input_pre}
            sed 's/:/\t/g' ${input_pre}.frq > ${input_pre}_Freq
            awk '{print $1,$2-1,$2,$5,$6,$7,$8}' OFS='\t' ${input_pre}_Freq | sed '1d' > ${input_pre}_bed
        else
            /public/software/env01/bin/vcftools --vcf ${input_pre}.recode.vcf --freq --keep ${list_file} --out ${panel_name}_${input_pre}
            sed 's/:/\t/g' ${panel_name}_${input_pre}.frq > ${panel_name}_${input_pre}_Freq
            awk '{print $1,$2-1,$2,$5,$6,$7,$8}' OFS='\t' ${panel_name}_${input_pre}_Freq | sed '1d' > ${panel_name}_${input_pre}_bed
        fi
    done
done

##filter and calculate maf
panels=("Can" "Lan" "PreCan")
input_name="Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001"

for panel in "${panels[@]}"; do
    cat ${panel}_chr{01..12}.DP4_100.GQ10.Q30.MR0.5.maf0.0001_bed > ${panel}_${input_name}_bed
    grep -v "nan" ${panel}_${input_name}_bed > ${panel}_${input_name}_bed.filter
    awk '{MAF=($5 < $7) ? $5 : $7; print $0,MAF;}' OFS='\t' ${panel}_${input_name}_bed.filter > ${panel}_${input_name}_bed.filter.maf
done















