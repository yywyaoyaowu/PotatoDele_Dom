#!/bin/bash
#SBATCH --partition=low,big,amd
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J work.sh
#SBATCH --error=err_%J_work.sh
#SBATCH --output=out_%J_work.sh
##################################################################
# @Author: huyong
# @Created Time : Sun Jul 21 11:47:04 2024

# @File Name: work.sh
# @Description:
##################################################################

for i in {01..12}
do
	echo "/public/software/env01/bin/vcftools --vcf chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf --window-pi 100000 --window-pi-step 10000 --out Can_chr${i} --keep Candolleanum_list" >> p
	echo "/public/software/env01/bin/vcftools --vcf chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf --window-pi 100000 --window-pi-step 10000 --out CanPre_chr${i} --keep CandolleanumPre_list" >> p2
	echo "/public/software/env01/bin/vcftools --vcf chr${i}.DP4_100.GQ10.Q30.MR0.5.maf0.001.recode.vcf --window-pi 100000 --window-pi-step 10000 --out Lan_chr${i} --keep Landrace_list" >> p
done
/public/software/env01/bin/parallel -j 24 < p
/public/software/env01/bin/parallel -j 12 < p2

for i in {02..12}; do grep -v "CHROM" Lan_chr${i}.windowed.pi >> Lan_chr01.windowed.pi; done
for i in {02..12}; do grep -v "CHROM" Can_chr${i}.windowed.pi >> Can_chr01.windowed.pi; done
for i in {02..12}; do grep -v "CHROM" CanPre_chr${i}.windowed.pi >> CanPre_chr01.windowed.pi; done

mv Can_chr01.windowed.pi Can.windowed.pi
mv Lan_chr01.windowed.pi Lan.windowed.pi
mv CanPre_chr01.windowed.pi CanPre.windowed.pi

##### remove other of diff
awk '{print $1"_"$2"_"$3}' Can.windowed.pi > Can.pi.win
awk '{print $1"_"$2"_"$3}' Lan.windowed.pi > Lan.pi.win
diff Can.pi.win Lan.pi.win | grep -v "d" | grep -v "a" | grep -v "-" | sed 's/< //g' | sed 's/> //g' > diff
awk '{print $1"_"$2"_"$3"\t"$0}' Can.windowed.pi | grep -vf diff > Can.windowed.pi1
awk '{print $1"_"$2"_"$3"\t"$0}' Lan.windowed.pi | grep -vf diff > Lan.windowed.pi1
mv Can.windowed.pi1 Can.windowed.pi
mv Lan.windowed.pi1 Lan.windowed.pi
cat <(echo -e "CHROM\tBIN_START\tBIN_END\tLan_PI\tCan_PI") <(paste <(cut -f2,3,4,6 Lan.windowed.pi) <(cut -f6 Can.windowed.pi) | grep -v CHROM) > Lan_Can.windowed.pi

Rscript stat.R
awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted > Lan_Can.windowed.pi.dom.sorted.filter
awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted | sed '1d' | wc -l  # 72679
awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted | head -n3634 > Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001
sort -k1,1 -k2n,2 Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001 | cut -f1,2,3 | sed '1d' > Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12
awk '{print $1"\t"$2-1"\t"$3}' Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12 > Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.bed
/public/software/env01/bin/bedtools merge -i Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.bed > Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.merge.bed
awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.merge.bed

# get gene in top 5% regions
/public/software/env01/bin/bedtools intersect -a Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.merge.bed -b ../Solanum_tuberosumDM.bed -wa -wb > gene_in_top0.05
awk '{print $1"_"$2"_"$3"\t"$0}' gene_in_top0.05 > gene_in_top0.05_2
python Tstat.py
cat <(cat header) <(awk -F '\t|_' '{print $1,$2,$3,$4,$5}' OFS='\t' gene_in_top0.05_merge.txt) > gene_in_dom_top0.05_merge.txt

cut -f4 gene_in_dom_top0.05_merge.txt | sed '1d' | sed 's/,/\n/g' | sort | uniq | wc -l
cut -f4 gene_in_dom_top0.05_merge.txt | sed '1d' | sed 's/,/\n/g' | sort | uniq > gene_in_dom_top0.05_merge_geneID

awk '{print $1"_"$2"_"$3"\t"$0}' CanPre.windowed.pi > CanPre.windowed.pi1                            
mv CanPre.windowed.pi1 CanPre.windowed.pi 

awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted > Lan_Can.windowed.pi.dom.sorted.filter

awk '{print $2,$3,$4,$1}' OFS='\t' ori.marker | sort -k1,1 -k2n,2 > marker
/public/software/env01/bin/bedtools intersect -b marker -a Lan_Can.windowed.pi.dom.sorted.top0.05.Filter0.0001.f12.merge.bed -wa -wb > marker_dom.sorted.top0.05.merge
cut -f4- marker_dom.sorted.top0.05.merge | sort | uniq > marker_in_dom_region

cat max | while read a b; do grep ${a} Lan.windowed.pi | grep ${b}; done
cat max | while read a b; do grep ${a} Can.windowed.pi | grep ${b}; done
awk '$5<=0.0001' Lan_Can.windowed.pi | cut -f1,2 > filter_out

###########
awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted | head -n7268 > Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001
sort -k1,1 -k2n,2 Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001 | cut -f1,2,3 | sed '1d' > Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12
awk '{print $1"\t"$2-1"\t"$3}' Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12 > Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12.bed
/public/software/env01/bin/bedtools merge -i Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12.bed > Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12.merge.bed
awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' Lan_Can.windowed.pi.dom.sorted.top0.10.Filter0.0001.f12.merge.bed

##########
awk '$5>0.0001' Lan_Can.windowed.pi.dom.sorted | head -n1454 > Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001
sort -k1,1 -k2n,2 Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001 | cut -f1,2,3 | sed '1d' > Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12
awk '{print $1"\t"$2-1"\t"$3}' Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12 > Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12.bed
/public/software/env01/bin/bedtools merge -i Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12.bed > Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12.merge.bed
awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' Lan_Can.windowed.pi.dom.sorted.top0.02.Filter0.0001.f12.merge.bed














