#!/bin/bash
#SBATCH --partition=low,big
#SBATCH -N 1
#SBATCH -c 5

CDS=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_CDS.sort.merg.bed
Intron=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_intron.sort.merge.bed
UTR=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_prime_UTR.sort.merge.bed
UpDown5K=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_UpDown5K.sort.merge.bed
promoter=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_Promoter_TSS1K.bed_merge
GeneticRegion=/public/home/yuqing/home/dzh/Potato/Threshold_Test/Input/Solanum_tuberosumDM_v6_repr_hc_genePlusUpDown5k.sort.merge.bed

for panel in can lan; do
    input_file=snpEff.anno.info.bed.${panel}
    gerp_file=${input_file}_GERP2
    
    ### different region number
    bedtools intersect -a $input_file -b $CDS > ${input_file}_CDS
    bedtools intersect -a $input_file -b $Intron > ${input_file}_Intron
    bedtools intersect -a $input_file -b $UTR > ${input_file}_UTR
    bedtools intersect -a $input_file -b $promoter > ${input_file}_Promoter_TSS1K
    bedtools intersect -a $input_file -b $GeneticRegion > ${input_file}_Genetic
    bedtools intersect -a $input_file -b $UpDown5K > ${input_file}_UpDown5K
    bedtools subtract -a $input_file -b $GeneticRegion > ${input_file}_InterGenetic
    
    bedtools intersect -a $gerp_file -b $CDS > ${gerp_file}_CDS
    bedtools intersect -a $gerp_file -b $Intron > ${gerp_file}_Intron
    bedtools intersect -a $gerp_file -b $UTR > ${gerp_file}_UTR
    bedtools intersect -a $gerp_file -b $promoter > ${gerp_file}_Promoter_TSS1K
    bedtools intersect -a $gerp_file -b $GeneticRegion > ${gerp_file}_Genetic
    bedtools intersect -a $gerp_file -b $UpDown5K > ${gerp_file}_UpDown5K
    bedtools subtract -a $gerp_file -b $GeneticRegion > ${gerp_file}_InterGenetic
    
    Poly_num=`cat $input_file | wc -l`
    CDS_Poly_num=`cat ${input_file}_CDS | wc -l`
    Intron_Poly_num=`cat ${input_file}_Intron | wc -l`
    UTR_Poly_num=`cat ${input_file}_UTR | wc -l`
    Promoter_TSS1K_Poly_num=`cat ${input_file}_Promoter_TSS1K | wc -l`
    Genetic_Poly_num=`cat ${input_file}_Genetic | wc -l`
    UpDown5K_Poly_num=`cat ${input_file}_UpDown5K | wc -l`
    InterGenetic_Poly_num=`cat ${input_file}_InterGenetic | wc -l`
    NonSynonymous_Site_num=`cat ${input_file}| egrep 'MODERATE|HIGH' | wc -l`
    Synonymous_Site_num=`grep synonymous_variant ${input_file} | wc -l`
    Synonymous_Poly_num=$Synonymous_Site_num
    NonSynonymous_Poly_num=$NonSynonymous_Site_num
    
    Poly_num_GERP2=`cat $gerp_file | wc -l`
    CDS_Poly_num_GERP2=`cat ${gerp_file}_CDS | wc -l`
    Intron_Poly_num_GERP2=`cat ${gerp_file}_Intron | wc -l`
    UTR_Poly_num_GERP2=`cat ${gerp_file}_UTR | wc -l`
    Promoter_TSS1K_Poly_num_GERP2=`cat ${gerp_file}_Promoter_TSS1K | wc -l`
    Genetic_Poly_num_GERP2=`cat ${gerp_file}_Genetic | wc -l`
    UpDown5K_Poly_num_GERP2=`cat ${gerp_file}_UpDown5K | wc -l`
    InterGenetic_Poly_num_GERP2=`cat ${gerp_file}_InterGenetic | wc -l`
    NonSynonymous_Site_num_GERP2=`cat ${gerp_file}| egrep 'MODERATE|HIGH' | wc -l`
    Synonymous_Site_num_GERP2=`grep synonymous_variant ${gerp_file} | wc -l`
    Synonymous_Poly_num_GERP2=$Synonymous_Site_num_GERP2
    NonSynonymous_Poly_num_GERP2=$NonSynonymous_Site_num_GERP2
    
    Poly_num_GERP2_75=`awk '($8>=2.75){print}' $gerp_file | wc -l`
    CDS_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_CDS | wc -l`
    Intron_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_Intron | wc -l`
    UTR_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_UTR | wc -l`
    Promoter_TSS1K_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_Promoter_TSS1K | wc -l`
    Genetic_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_Genetic | wc -l`
    UpDown5K_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_UpDown5K | wc -l`
    InterGenetic_Poly_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}_InterGenetic | wc -l`
    NonSynonymous_Site_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file}| egrep 'MODERATE|HIGH' | wc -l`
    Synonymous_Site_num_GERP2_75=`awk '($8>=2.75){print}' ${gerp_file} | grep synonymous_variant | wc -l`
    Synonymous_Poly_num_GERP2_75=$Synonymous_Site_num_GERP2_75
    NonSynonymous_Poly_num_GERP2_75=$NonSynonymous_Site_num_GERP2_75
    
    Poly_num_GERP3_5=`awk '($8>=3.5){print}' $gerp_file | wc -l`
    CDS_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_CDS | wc -l`
    Intron_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_Intron | wc -l`
    UTR_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_UTR | wc -l`
    Promoter_TSS1K_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_Promoter_TSS1K | wc -l`
    Genetic_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_Genetic | wc -l`
    UpDown5K_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_UpDown5K | wc -l`
    InterGenetic_Poly_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}_InterGenetic | wc -l`
    NonSynonymous_Site_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file}| egrep 'MODERATE|HIGH' | wc -l`
    Synonymous_Site_num_GERP3_5=`awk '($8>=3.5){print}' ${gerp_file} | grep synonymous_variant | wc -l`
    Synonymous_Poly_num_GERP3_5=$Synonymous_Site_num_GERP3_5
    NonSynonymous_Poly_num_GERP3_5=$NonSynonymous_Site_num_GERP3_5
    
    DeleStat_output=${panel}_DeleStat_DiffCutoff.txt
    
    echo All $Poly_num > ${DeleStat_output}
    echo CDS $CDS_Poly_num >> ${DeleStat_output}
    echo Intron $Intron_Poly_num >> ${DeleStat_output}
    echo UTR $UTR_Poly_num >> ${DeleStat_output}
    echo Promoter $Promoter_TSS1K_Poly_num >> ${DeleStat_output}
    echo Updown5K $UpDown5K_Poly_num >> ${DeleStat_output}
    echo InterGenetic $InterGenetic_Poly_num >> ${DeleStat_output}
    echo Synonymous $Synonymous_Poly_num >> ${DeleStat_output}
    echo NonSynonymous $NonSynonymous_Site_num >> ${DeleStat_output}
    
    echo All_GERP2 $Poly_num_GERP2 >> ${DeleStat_output}
    echo CDS_GERP2 $CDS_Poly_num_GERP2 >> ${DeleStat_output}
    echo Intron_GERP2 $Intron_Poly_num_GERP2 >> ${DeleStat_output}
    echo UTR_GERP2 $UTR_Poly_num_GERP2 >> ${DeleStat_output}
    echo Promoter_GERP2 $Promoter_TSS1K_Poly_num_GERP2 >> ${DeleStat_output}
    echo Updown5K_GERP2 $UpDown5K_Poly_num_GERP2 >> ${DeleStat_output}
    echo InterGenetic_GERP2 $InterGenetic_Poly_num_GERP2 >> ${DeleStat_output}
    echo Synonymous_GERP2 $Synonymous_Poly_num_GERP2 >> ${DeleStat_output}
    echo NonSynonymous_GERP2 $NonSynonymous_Site_num_GERP2 >> ${DeleStat_output}
    
    echo All_GERP2_75 $Poly_num_GERP2_75 >> ${DeleStat_output}
    echo CDS_GERP2_75 $CDS_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo Intron_GERP2_75 $Intron_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo UTR_GERP2_75 $UTR_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo Promoter_GERP2_75 $Promoter_TSS1K_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo Updown5K_GERP2_75 $UpDown5K_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo InterGenetic_GERP2_75 $InterGenetic_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo Synonymous_GERP2_75 $Synonymous_Poly_num_GERP2_75 >> ${DeleStat_output}
    echo NonSynonymous_GERP2_75 $NonSynonymous_Site_num_GERP2_75 >> ${DeleStat_output}
    
    echo All_GERP3_5 $Poly_num_GERP3_5 >> ${DeleStat_output}
    echo CDS_GERP3_5 $CDS_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo Intron_GERP3_5 $Intron_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo UTR_GERP3_5 $UTR_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo Promoter_GERP3_5 $Promoter_TSS1K_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo Updown5K_GERP3_5 $UpDown5K_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo InterGenetic_GERP3_5 $InterGenetic_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo Synonymous_GERP3_5 $Synonymous_Poly_num_GERP3_5 >> ${DeleStat_output}
    echo NonSynonymous_GERP3_5 $NonSynonymous_Site_num_GERP3_5 >> ${DeleStat_output}
done
