#!/bin/bash -l
#$ -cwd 
#$ -N STEP4_genotype_MAGIC
#$ -o STEP4_genotype_MAGIC.log
#$ -l mem=20G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -t 1-43
#request TMPDIR space
#$ -l tmpfs=70G

set -a

source ../Input/data_and_variables.sh
source ../Input/software.sh

#haplotypeCaller options
java_memory_param=12
verbosity_param=ERROR
interval_padding_param=100
minimum_mapping_quality_param=30

#first convert to get chr, genome and half
#get chr, genome, and half
get_wheat_chr_genome_half $SGE_TASK_ID

#log and directory for temporary files
log=${project_dir}/scripts/2Variants/log_MAGIC_genotype_given_alleles
tmpdir=${TMPDIR}/tmp_job$SGE_TASK_ID 
mkdir -p ${tmpdir}
mkdir -p ${log}
mkdir -p $MAGIC_vcf_dir

samplename=MAGIC_combined

#create bamlist for MAGIC lines
bamlist_file=${TMPDIR}/MAGIC_bamlist.txt
ls ${MAGIC_alignments_dir}/*/*.bam > ${bamlist_file}

#we remove lines with low concordance or nonMAGIC lines, lists found in ${array_data_dir}
while read line; do
grep -v $line $bamlist_file > ${bamlist_file}.tmp
mv ${bamlist_file}.tmp $bamlist_file
done < ${array_data_dir}/lines_to_ignore_0.95.txt
#remove lines not to be included at all
while read line; do
grep -v $line $bamlist_file > ${bamlist_file}.tmp
mv ${bamlist_file}.tmp $bamlist_file
done < ${array_data_dir}/lines_to_ignore_nonMAGIC.txt

#have to make an arguments file that takes in all the relevant g.vcf's
input_files_arguments=${TMPDIR}/${samplename}.arguments.txt
while read line ; do
echo "--input $line"
done < $bamlist_file > ${input_files_arguments}

#get vcf for chromosome from founder calls
known_variants=${founder_vcf_dir}/nohet_biallelic_polymorphic_filter/Founders_combined_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf

#must first index the vcf at which variants will be called
gatk IndexFeatureFile --java-options "-Xmx${java_memory_param}G" \
--feature-file $known_variants 

mkdir -p ${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter
#haplotype caller at known sites
gatk HaplotypeCaller --java-options "-Xmx${java_memory_param}G" \
--reference ${ref} \
--TMP_DIR=${tmpdir} \
--arguments_file $input_files_arguments \
--alleles $known_variants \
--intervals ${known_variants} \
--interval-padding $interval_padding_param \
--minimum-mapping-quality ${minimum_mapping_quality_param} \
--verbosity $verbosity_param \
--genotyping-mode GENOTYPE_GIVEN_ALLELES \
--genotype-filtered-alleles true \
--output-mode EMIT_ALL_SITES \
--standard-min-confidence-threshold-for-calling 0 \
--output ${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter/${samplename}_chr${chr}${genome}${half}.vcf 2> ${log}/${samplename}_chr${chr}${genome}${half}.HaplotypeCaller.log

