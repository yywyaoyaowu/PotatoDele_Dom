#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J work.sh
##################################################################
# @Author: huyong
# @Created Time : Sun Aug 28 11:26:41 2022

# @File Name: work.sh
# @Description:
##################################################################

head -n12 ../Solanum_tuberosumDM.fa.fai | while read chr len orther
do
   mkdir  ${chr}
   cd ${chr}

  for i in `seq 1 10000000 ${len}`
   do
       let j=${i}+9999999
       
       if (( ${j} > ${len} ))
       then
           j=${len}
       fi

       echo """#!/bin/bash
#SBATCH --partition=low,big
#SBATCH -N 1
#SBATCH -c 2
/home/huyong/software/gatk-4.1.9.0/gatk GenomicsDBImport -R ../../Solanum_tuberosumDM.fa -L ${chr}:${i}-${j} --sample-name-map ../gvcf_path --genomicsdb-workspace-path ${chr}_${i}.db
/home/huyong/software/gatk-4.1.9.0/gatk GenotypeGVCFs -R ../../Solanum_tuberosumDM.fa --sample-ploidy 2 --max-genotype-count 8192 -O ${chr}_${i}.combined.vcf -V gendb://${chr}_${i}.db -L ${chr}:${i}-${j}""" >>  ${chr}_${i}.sh
       
       sbatch ${chr}_${i}.sh
   done

   cd ..
done

