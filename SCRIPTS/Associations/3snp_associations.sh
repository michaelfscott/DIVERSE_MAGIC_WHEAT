#!/bin/bash -l
#$ -cwd 
#$ -N STEP3_snp_assocs
#$ -o STEP3_snp_assocs.log
#$ -l mem=5G
#$ -l h_rt=12:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=10G
#$ -t 1-73

set -a 
source ../Input/data_and_variables.sh
module unload compilers
module unload mpi 
module load r/recommended
source ../Input/software.sh

#for counter in {1..5}; do 

#TEMP=$(( $SGE_TASK_ID - 1 ))
#TASK_ID=$(( $TEMP * 5 + $counter ))
TASK_ID=$(( $SGE_TASK_ID + 2))

phenotype_name=$(head -n 1 ../phenotypes_2020/founder_and_MAGIC_phenotypes_with_header.tsv | cut -f $TASK_ID)

Rscript perform.mixed.model.snps.r ${phenotype_name} > log_snp/${TASK_ID}.${phenotype_name}.log

#done


