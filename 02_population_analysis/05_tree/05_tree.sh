#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Wed Jul 30 14:54:21 2025

# @File Name: work.sh
# @Description:
##################################################################

for i in {01..12}
do
    echo """~/software/snpEff/scripts/snpEff ann Solanum_tuberosumDM chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf > chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf""" >> split_anno
done
/public/software/env01/bin/parallel -j 12 < split_anno

for i in {01..12}
do
    echo """python calc_4dTv_in_eff_vcf.py chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf /public/agis/huangsanwen_group/huyong/Cando_del/Solanum_tuberosumDM.fa""" >> split_calc_4dTv
done
/public/software/env01/bin/parallel -j 12 < split_calc_4dTv

cp chr01.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf
for i in {02..12}
do
       grep -v "^#" chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf >> DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf
done
python vcf2phylip.py -i DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.vcf --output-folder .
/public/software/env01/bin/iqtree -s DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy -nt AUTO -m GTR -bb 1000

cut -f1-8 chr01.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf > DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info
for i in {02..12}
do
	cut -f1-8 chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.vcf | grep -v "^#" >> DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.info
done

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep Candolleanum_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree Candolleanum_tree

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep Old_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree previous_Candolleanum_Landrace_tree

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep New_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree NewCandolleanum_tree

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep Landrace_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree Landrace_tree

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep CandolleanumPre_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree CandolleanumPre_tree

java -jar /home/wuyaoyao/software/PareTree1.0.2.jar -keep potato_list -t O -f DP4_100.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree
mv DP4_100_pared.GQ10.Q30.MR0.5.maf0.001.anno.4dTv.min4.phy.contree potato_tree


grep "New_candolleanum" ../sample_info_final | cut -f1 | awk '{print $0",branch,node,#F0A419,18,bold-italic,#F0A419"}' > cando_style
grep "Previous_candolleanum" ../sample_info_final | cut -f1 | awk '{print $0",branch,node,#F3C97F,18,bold-italic,#F3C97F"}' >> cando_style
grep -v "candolleanum" ../sample_info_final | grep "stenotomum_subsp_stenotomum" | cut -f1 | awk '{print $0",branch,node,#8CC63E,18,bold-italic,#8CC63E"}' > stenotomum_style
grep -v "candolleanum" ../sample_info_final | grep "stenotomum_subsp.goniocalyx" | cut -f1 | awk '{print $0",branch,node,#8CA0CB,18,bold-italic,#8CA0CB"}' > goniocalyx_style
grep -v "candolleanum" ../sample_info_final | grep "phureja" | cut -f1 | awk '{print $0",branch,node,#56B4E9,18,bold-italic,#56B4E9"}' > phureja_style
grep -v "candolleanum" ../sample_info_final | grep "ajanhuiri" | cut -f1 | awk '{print $0",branch,node,#009E73,18,bold-italic,#009E73"}' > ajanhuiri_style
egrep "etuberosum|founder|inbred" ../sample_info_final | cut -f1 | awk '{print $0",branch,node,#000000,18,bold-italic,#000000"}' > other_style
cat styles cando_style stenotomum_style goniocalyx_style phureja_style ajanhuiri_style other_style > style.txt

grep "New_candolleanum" ../sample_info_final | cut -f1 | awk '{print $0" #F0A419 COL#F0A419"}' > cando_srtip
grep "Previous_candolleanum" ../sample_info_final | cut -f1 | awk '{print $0" #F3C97F COL#F3C97F"}' >> cando_srtip
grep "stenotomum_subsp_stenotomum" ../sample_info_final | cut -f1 | awk '{print $0" #8CC63E COL#8CC63E"}' > stenotomum_srtip
grep "stenotomum_subsp.goniocalyx" ../sample_info_final | cut -f1 | awk '{print $0" #8CA0CB COL#8CA0CB"}' > goniocalyx_srtip
grep "phureja" ../sample_info_final | cut -f1 | awk '{print $0" #56B4E9 COL#56B4E9"}' > phureja_srtip
grep "ajanhuiri" ../sample_info_final | cut -f1 | awk '{print $0" #009E73 COL#009E73"}' > ajanhuiri_srtip
egrep "etuberosum|founder|inbred" ../sample_info_final | cut -f1 | awk '{print $0" #000000 COL#000000"}' > other_srtip
cat color_strip cando_srtip stenotomum_srtip goniocalyx_srtip phureja_srtip ajanhuiri_srtip other_srtip > color_strip.txt


# only split cando and Landrace
egrep "stenotomum|goniocalyx|phureja|ajanhuiri" ../sample_info_final | cut -f1 | awk '{print $0",branch,node,#548B54,18,bold-italic,#548B54"}' > Landrace_style
egrep "etuberosum" ../sample_info_final | cut -f1 | awk '{print $0",branch,node,#000000,18,bold-italic,#000000"}' > other_style
cat styles_split_CL cando_style Landrace_style other_style > style_2.txt

egrep "stenotomum|goniocalyx|phureja|ajanhuiri" ../sample_info_final | cut -f1 | awk '{print $0" #548B54 COL#548B54"}' > Landrace_srtip
egrep "etuberosum" ../sample_info_final | cut -f1 | awk '{print $0" #000000 COL#000000"}' > other_srtip
cat color_strip cando_srtip Landrace_srtip other_srtip > color_strip_2.txt




