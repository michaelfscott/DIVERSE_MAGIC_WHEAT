#!/bin/bash -l
#$ -cwd 
#$ -N STEP0_join_chr_halves_pos_ref
#$ -o STEP0_join_chr_halves_pos_ref.log
#$ -l mem=5G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y

set -a

#to deal with bam file specification errors in GATK caused by long chromosomes, each chromosome on the reference genome was split in half before alignment. 
#this script joins the variants called on each half of each chromosome so that we run imputation on the full chromosome. 

source ../Input/data_and_variables.sh
#load R
module unload compilers
module unload mpi
module load r/recommended

Rscript merge_SNPlist.R ${founder_SNPlist_dir} ${founder_reference_panel}

for index in {1..21}; do
get_wheat_chr_genome $index
#gzip legend output from above rscript
gzip -f ${founder_reference_panel}/Founders_combined_chr${chr}${genome}.legend
#concatonate hap files
zcat ${founder_reference_panel}/Founders_combined_chr${chr}${genome}_1.hap.gz \
${founder_reference_panel}/Founders_combined_chr${chr}${genome}_2.hap.gz | \
gzip -c > ${founder_reference_panel}/Founders_combined_chr${chr}${genome}.hap.gz
cp ${founder_reference_panel}/Founders_combined_chr${chr}${genome}_1.sample ${founder_reference_panel}/Founders_combined_chr${chr}${genome}.sample
done

