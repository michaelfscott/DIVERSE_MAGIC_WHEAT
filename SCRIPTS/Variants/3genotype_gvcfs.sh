#!/bin/bash -l
#$ -cwd 
#$ -N STEP2_genotype_gvcfs
#$ -o STEP2_genotype_gvcfs.log
#$ -l mem=20G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -t 1-43
#$ -l tmpfs=70G

set -a
source ../Input/data_and_variables.sh
source ../Input/software.sh

#haplotypeCaller options
interval_padding_param=100
java_memory_param=12
minimum_mapping_quality_param=30
sample_ploidy_param=2
verbosity_param=ERROR

#first convert to get chr, genome and half
#get chr, genome, and half
get_wheat_chr_genome_half $SGE_TASK_ID

#log and directory for temporary files
log=${project_dir}/scripts/2Variants/log_Founders_combine_gvcfs_MQ30 
tmpdir=${TMPDIR}/tmp_job$SGE_TASK_ID 
mkdir -p ${tmpdir}
mkdir -p ${log}

#have to make an arguments file that takes in all the relevant g.vcf's
ls ${founder_gvcfs_dir}/*/*chr${chr}${genome}${half}.g.vcf > ${TMPDIR}/tmp.txt
input_files_arguments=${TMPDIR}/Founders_chr${chr}${genome}${half}.arguments.txt
while read line ; do
echo "-V $line"
done < ${TMPDIR}/tmp.txt > ${input_files_arguments}

#make a target bed specific to the chr half
target_bed=${TMPDIR}/gene_and_promoter_targets_HC_and_5UTR.edit.chr${chr}${genome}${half}.bed
if [[ $half == "_1" ]]; then
  grep ^${chr}${genome}:1- ${gene_and_prom_and_5UTR_bed} > ${target_bed}
else
  grep ^${chr}${genome} ${gene_and_prom_and_5UTR_bed} | \
  grep -v ^${chr}${genome}:1- > ${target_bed}
fi

samplename=Founders_combined

echo "$SGE_TASK_ID combine GVCFs for ${samplename} started for chr${chr}${genome} at $(date)"
gatk CombineGVCFs --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=$tmpdir \
--arguments_file $input_files_arguments \
--reference $ref \
--intervals $target_bed \
--interval-padding $interval_padding_param \
--verbosity $verbosity_param \
--output ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.combineGVCFs.log

echo "copying combined GVCF for ${samplename} chr${chr}${genome}${half} at $(date)"
mkdir -p ${founder_gvcfs_dir}/$samplename
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf* --target-directory=${founder_gvcfs_dir}/${samplename}

input_GVCF=${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf

echo "$TASK_ID genotype GVCFs for ${samplename} started for chr${chr}${genome}${half} at $(date)"
gatk GenotypeGVCFs --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=$tmpdir \
--variant $input_GVCF \
--reference $ref \
--intervals $target_bed \
--interval-padding $interval_padding_param \
--verbosity $verbosity_param \
--output ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.genotypeGVCFs.log

echo "copying combined VCF for ${samplename} chr${chr}${genome}${half} at $(date)"
mkdir -p ${founder_vcf_dir}/raw_vcfs
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.vcf* --target-directory=${founder_vcf_dir}/raw_vcfs

