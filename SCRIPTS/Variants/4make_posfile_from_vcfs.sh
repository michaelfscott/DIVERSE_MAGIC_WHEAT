#!/bin/bash -l
#$ -cwd 
#$ -N STEP3_make_posfile_vcfs_filter
#$ -o STEP3_make_posfile_vcfs_filter.log
#$ -l mem=5G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=10G

set -a
source ../Input/data_and_variables.sh
source ../Input/software.sh

#filtering parameters
samplename=Founders_combined
min_mean_DP_param=5
max_mean_DP_param=60
min_DP_param=2
max_DP_param=120
minQ_param=100
max_missing_param=1

for index in {1..43}; do 
#first convert to get chr, genome and half
#get chr, genome, and half
get_wheat_chr_genome_half $index

input_VCF=${founder_vcf_dir}/raw_vcfs/${samplename}_chr${chr}${genome}${half}.vcf

echo "$SGE_TASK_ID vcftools for ${samplename} started for chr${chr}${genome}${half} at $(date)"
vcftools --vcf $input_VCF \
--recode --recode-INFO-all \
--remove-indels \
--min-alleles 2 --max-alleles 2 \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--max-missing ${max_missing_param} \
--maf 0.01 \
--out ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs

bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
--output-type v \
--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.filter.vcf

mkdir -p ${founder_vcf_dir}/biallelic_polymorphic
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf* --target-directory=${founder_vcf_dir}/biallelic_polymorphic
mkdir -p ${founder_vcf_dir}/biallelic_polymorphic_filter
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.filter.vcf* --target-directory=${founder_vcf_dir}/biallelic_polymorphic_filter

#remove calls with heterozygotes

bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
--genotype ^het \
--output-type v \
--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.vcf

bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
--genotype ^het \
--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
--output-type v \
--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf

echo "copying nohet vcf files for ${samplename} chr${chr}${genome}${half} at $(date)"
mkdir -p ${founder_vcf_dir}/nohet_biallelic_polymorphic
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.vcf* --target-directory=${founder_vcf_dir}/nohet_biallelic_polymorphic
mkdir -p ${founder_vcf_dir}/nohet_biallelic_polymorphic_filter
cp ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf* --target-directory=${founder_vcf_dir}/nohet_biallelic_polymorphic_filter

grep -v "^#" ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf |
awk -v OFS='\t' '{print $1, $2, $4, $5}' > ${TMPDIR}/chr${chr}${genome}${half}.SNPs.txt

echo "copying SNPlist files for ${samplename} chr${chr}${genome}${half} at $(date)"
mkdir -p ${founder_SNPlist_dir}
cp ${TMPDIR}/chr${chr}${genome}${half}.SNPs.txt --target-directory=${founder_SNPlist_dir}

done

