#!/bin/bash -l
#$ -cwd
#$ -N STEP2_info_filter
#$ -o STEP2_info_filter.log
#memory and runtime options 
#$ -l mem=2G
#Number of threads per job, memory given above is per-thread
#$ -l h_rt=12:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=10G

#This script uses the Rscript infoScoreFilter.R to determine the sites that don't pass the input score threshold for STITCH

set -a
source ../Input/data_and_variables.sh
#source ../Input/software.sh
#for this, we will load r:
module unload compilers
module unload mpi
module load r/recommended

for index in {1..21}; do 

get_wheat_chr_genome $index

chromosome_name=chr${chr}${genome}

chromosome_dir=${STITCH_output_dir}/${chromosome_name}
input_RData=${chromosome_dir}/RData/EM.all.${chr}${genome}.RData
output_file_name=${chromosome_dir}/Exclude_Info_${infoScoreFilter}.txt

echo $chromosome_name
Rscript infoScoreFilter.R $input_RData $infoScoreFilter $output_file_name

done

cat ${STITCH_output_dir}/chr*/Exclude_Info_${infoScoreFilter}.txt > ${STITCH_output_dir}/Exclude_Info_${infoScoreFilter}_allChr.txt

