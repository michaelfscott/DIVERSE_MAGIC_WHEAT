#!/bin/bash -l
#$ -cwd 
#$ -N run_STITCH_join_halves
#$ -o run_STITCH_join_halves.log
#$ -l mem=5G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=50G
#$ -t 1-21

set -a

source ../Input/data_and_variables.sh
#for this, we will load r:
module unload compilers
module unload mpi
#load R 3.5.1
module load r/recommended
#bgzip also required for STITCH
PATH=/home/ucbtmf1/packages/STITCH:$PATH

get_wheat_chr_genome $SGE_TASK_ID

#create bamlist for MAGIC lines
bamlist_file=${TMPDIR}/MAGIC_bamlist.txt
ls ${MAGIC_alignments_dir}/*/*.bam > ${bamlist_file}
#before running this script, we create a list of lines to remove by comparing the genotypes from the MAGIC lines with the corresponding array genotype calls and removing all lines with a concordance of <0.96. 
#here, we remove these lines, found in ${array_data_dir}/lines_to_ignore_0.96.txt
while read line; do 
grep -v $line $bamlist_file > ${bamlist_file}.tmp
mv ${bamlist_file}.tmp $bamlist_file
done < ${array_data_dir}/lines_to_ignore_0.95.txt

#remove lines not to be included at all
while read line; do
grep -v $line $bamlist_file > ${bamlist_file}.tmp
mv ${bamlist_file}.tmp $bamlist_file
done < ${array_data_dir}/lines_to_ignore_nonMAGIC.txt

tmpdir=${TMPDIR}/tmpdir_${SGE_TASK_ID}
mkdir -p ${tmpdir}
outdir_STITCH=${TMPDIR}/output
mkdir -p ${outdir_STITCH}

K_param=16
G_param=3
chromosome=$chr$genome
chromosome_name=chr$chr$genome
bqFilter_param=30
minRate_param=0.001
bamlist_file=${bamlist_file}
posfile_dir=${founder_SNPlist_dir}
reference_panel_prefix=${founder_reference_panel}/Founders_combined_${chromosome_name}
reference_phred=30
reference_iterations=0
niterations=40
tmpdir_STITCH=${tmpdir}
outdir_STITCH=${outdir_STITCH}

log=${project_dir}/scripts/3STITCH_downsample/log_STITCH
mkdir -p ${log}
mkdir -p ${log}_input

### First run STITCH in generateInputOnly mode on each chromosome half

for half in _1 _2; do 

chromosome_half_name=$(grep ${chromosome_name}${half} $chr_names_lengths_file | awk '{print $2}')

Rscript run_stitch_generateInputOnly.R $K_param \
$G_param \
$chromosome_half_name \
${chromosome_name}${half} \
$bqFilter_param \
$minRate_param \
$bamlist_file \
$posfile_dir \
${reference_panel_prefix}${half} \
$reference_phred \
$reference_iterations \
$niterations \
$tmpdir_STITCH \
$outdir_STITCH > ${log}_input/stitch.${chromosome_name}${half}.log 2>&1 

done

sampleReads_1_path=$(ls ${outdir_STITCH}/${chromosome_name}_1/input/sample*.RData | head -n 1)
sampleReads_2_path=$(ls ${outdir_STITCH}/${chromosome_name}_2/input/sample*.RData | head -n 1)

### Then use the merge_posfile.r script to join the chromosome halves in the STITCH RData files

mkdir -p ${outdir_STITCH}/${chromosome_name}/RData
mkdir -p ${outdir_STITCH}/${chromosome_name}/input
sampleNames=$(ls ${outdir_STITCH}/${chromosome_name}_1/RData/sampleNames*.RData)
posfile_1=$(ls ${outdir_STITCH}/${chromosome_name}_1/RData/pos*.RData)
posfile_2=$(ls ${outdir_STITCH}/${chromosome_name}_2/RData/pos*.RData)
Rscript merge_posfile.r $sampleNames $posfile_1 $posfile_2 $sampleReads_1_path $sampleReads_2_path $chromosome_name $outdir_STITCH/${chromosome_name}

### Then run STITCH on the full chromosomes. 

Rscript run_stitch_ref.input.R $K_param \
$G_param \
$chromosome \
$chromosome_name \
$bqFilter_param \
$minRate_param \
$bamlist_file \
$posfile_dir \
$reference_panel_prefix \
$reference_phred \
$reference_iterations \
$niterations \
$tmpdir_STITCH \
$outdir_STITCH > ${log}/stitch.${chromosome_name}.log 2>&1

mkdir -p ${STITCH_output_dir}
cp -r ${outdir_STITCH}/${chromosome_name} --target-directory=${STITCH_output_dir}

echo "$SGE_TASK_ID ended at $(date)"

