#!/bin/bash -l
#$ -cwd 
#$ -N STEP1_call_gvcfs_Founders
#$ -o STEP1_call_gvcfs_Founders.log
#$ -l mem=4G
#$ -l h_rt=10:00:00 
#$ -S /bin/bash
#$ -j y
#$ -t 1-688
#$ -l tmpfs=70G

set -a
source ../Input/data_and_variables.sh
source ../Input/software.sh

#haplotypeCaller options
interval_padding_param=100
minimum_mapping_quality_param=30
sample_ploidy_param=2
java_memory=2
verbosity_param=ERROR

#for each of the chromosomes and each of the founders, we will run a separate job
echo "$SGE_TASK_ID $(date) on $(hostname)" #so that you know what time each task began and finished

#first convert SGE_TASK_ID to get a unique ID repeated across the founders
founder_index=$((($SGE_TASK_ID+15)%16 + 1))
chr_index=$((($SGE_TASK_ID +15)/16))

#get chr, genome, and half
get_wheat_chr_genome_half $chr_index
#get samplename for this job
samplename=$(get_founder_samplename $founder_index)
echo "$SGE_TASK_ID founder ${founder_index} $samplename chr $chr subgenome ${genome} half ${half}"

#log and directory for temporary files
log=${project_dir}/scripts/2Variants/log_Founders_gvcfs_MQ30 
tmpdir=${TMPDIR}/tmp_job$SGE_TASK_ID 
mkdir -p ${tmpdir}
mkdir -p ${log}

target_bed=${TMPDIR}/gene_and_promoter_targets_HC_and_5UTR.edit.chr${chr}${genome}${half}.bed
if [[ $half == "_1" ]]; then
  grep ^${chr}${genome}:1- ${gene_and_prom_and_5UTR_bed} > ${target_bed}
else
  grep ^${chr}${genome} ${gene_and_prom_and_5UTR_bed} | \
  grep -v ^${chr}${genome}:1- > ${target_bed}
fi

#only one input file for Holdfast but there is are two input bams for other founders (one from each capture)
if [[ $samplename == "Holdfast" ]]; then

input_bam=$(ls ${founder_alignments_dir}/${samplename}/*.dedup.bam | sed "1q;d")

echo "$SGE_TASK_ID HaplotypeCaller started for $samplename chr${chr}${genome}${half} at $(date)"
gatk HaplotypeCaller --java-options "-Xmx${java_memory}G" \
--TMP_DIR=$tmpdir/${samplename} \
--reference $ref \
--input $input_bam \
--intervals $target_bed \
--sample-ploidy $sample_ploidy_param \
--interval-padding $interval_padding_param \
--minimum-mapping-quality $minimum_mapping_quality_param \
--emit-ref-confidence GVCF \
--verbosity $verbosity_param \
--output ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.HaplotypeCaller.gvcf.log

else

input_bam1=$(ls ${founder_alignments_dir}/${samplename}_L1/*.dedup.2.bam | sed "1q;d")
input_bam2=$(ls ${founder_alignments_dir}/${samplename}_L2/*.dedup.2.bam | sed "1q;d")

echo "$SGE_TASK_ID HaplotypeCaller started for $samplename chr${chr}${genome}${half} at $(date)"
gatk HaplotypeCaller --java-options "-Xmx${java_memory}G" \
--TMP_DIR=$tmpdir/${samplename} \
--reference $ref \
--input $input_bam1 \
--input $input_bam2 \
--intervals $target_bed \
--sample-ploidy $sample_ploidy_param \
--interval-padding $interval_padding_param \
--minimum-mapping-quality $minimum_mapping_quality_param \
--emit-ref-confidence GVCF \
--verbosity $verbosity_param \
--output ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.HaplotypeCaller.gvcf.log

fi

echo "$SGE_TASK_ID indexing vcf for $samplename chr ${chr}${genome}${half} at $(date)" 
gatk IndexFeatureFile --java-options "-Xmx${java_memory}G" \
--feature-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.index.vcf.log

mkdir -p ${founder_gvcfs_dir}/${samplename}
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}.g.vcf* --target-directory=${founder_gvcfs_dir}/${samplename}

