#!/bin/bash -l
#$ -cwd
#$ -N STEP2_generate_haps
#$ -o STEP2_generate_haps.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=50G

set -a 
source ../Input/data_and_variables.sh
module unload compilers
module unload mpi
module load r/recommended

mkdir ${STITCH_haplotype_dir}/LD_pruned_haplotypes
perl vcf2df.pl 1
Rscript import.stitch.dosages.R ${STITCH_haplotype_dir}/LD_pruned_haplotypes

