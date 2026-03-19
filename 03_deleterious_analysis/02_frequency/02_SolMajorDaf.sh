#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 5

#SolMajor
panels=("Can" "Lan")
input_path="/public/home/yuqing/home/dzh/Potato/NewGeno"
input_name="DP4_100.GQ10.Q30.MR0.5.maf0.0001"

for i in {01..12}; do
    for panel in "${panels[@]}"; do
        bedtools intersect -a ${panel}_chr${i}.${input_name}_bed \
            -b ${input_path}/SolMsa_100Species_chr${i}.fa_BiAllelefreq.bed \
            -wa -wb | cut -f 1,2,3,4,5,6,7,15,16 > ${panel}_chr${i}.${input_name}_bed_SolMajor
    done
done

for panel in "${panels[@]}"; do
    cat ${panel}_chr{01..12}.${input_name}_bed_SolMajor > ${panel}_Allchrs.${input_name}_bed_SolMajor
    awk '{Major=($5 > $7) ? $4 : $6; print $0,Major;}' OFS='\t' \
        ${panel}_Allchrs.${input_name}_bed_SolMajor > ${panel}_Allchrs.${input_name}_bed_SolMajor.major
    wc -l ${panel}_Allchrs.${input_name}_bed_SolMajor.major
done

paste Can_Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_bed_SolMajor.major Lan_Allchrs.DP4_100.GQ10.Q30.MR0.5.maf0.0001_bed_SolMajor.major > Can.Lan.SolMajor.Major.bed

##all snp daf
awk '{
    RefAllele = $4
    RefAlleleFreq = $5
    AltAllele = $6
    AltAlleleFreq = $7
    Sol100Major = $8
    Sol100MajorCount = $9
    
    if (Sol100MajorCount >= 3) {
        if (Sol100Major == RefAllele) {
            daf = AltAlleleFreq
        } else if (Sol100Major == AltAllele) {
            daf = RefAlleleFreq
        } else {
            daf = 1
        }
    }
    
    else {
        daf = (RefAlleleFreq < AltAlleleFreq) ? RefAlleleFreq : AltAlleleFreq
    }
   
    print $0 "\t" daf
}' Can.Lan.SolMajor.Major.bed > Can.Lan.SolMajor.Major.Candaf.bed

awk '{
    CanMajor = $10    
    RefAllele = $14    
    RefAlleleFreq = $15    
    AltAllele = $16    
    AltAlleleFreq = $17    
    Sol100Major = $18    
    Sol100MajorCount = $19    
    LanMajor = $20    

    if (Sol100MajorCount >= 3) {
        if (Sol100Major == RefAllele) {
            daf = AltAlleleFreq
        } else if (Sol100Major == AltAllele) {
            daf = RefAlleleFreq
        } else {
            daf = 1
        }
    }
    
    else if (Sol100MajorCount < 3 && CanMajor == LanMajor) {
        daf = (RefAlleleFreq < AltAlleleFreq) ? RefAlleleFreq : AltAlleleFreq
    }
    
    else if (Sol100MajorCount < 3 && CanMajor != LanMajor) {
        daf = (RefAlleleFreq > AltAlleleFreq) ? RefAlleleFreq : AltAlleleFreq
    }
    
    print $0 "\t" daf
}' Can.Lan.SolMajor.Major.Candaf.bed > Can.Lan.SolMajor.Major.Candaf.Landaf.bed



