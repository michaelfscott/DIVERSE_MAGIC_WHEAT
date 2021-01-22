#!/bin/bash -l
#$ -cwd
#$ -N STEP1_LD_prune
#$ -o STEP1_LD_prune.log
#$ -l mem=3G
#$ -l h_rt=05:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=10G

set -a 
source ../Input/data_and_variables.sh
source ../Input/software.sh

input_PLINK_prefix=${STITCH_PLINK_dir}/ALLchr.SNPs
reference_alleles_ALLchr=${founder_plink_dir_merged}/alternate_allele_maps/Founders_combined_ALLchr.txt

plink \
--file $input_PLINK_prefix \
--indep-pairwise ${window_size} ${window_step} ${LD_threshold} \
--chr-set 21 \
--reference-allele $reference_alleles_ALLchr \
--out ${TMPDIR}/ALLchr.SNPs

plink \
--file ${input_PLINK_prefix} \
--make-bed \
--extract ${TMPDIR}/ALLchr.SNPs.prune.in \
--chr-set 21 \
--reference-allele $reference_alleles_ALLchr \
--out ${TMPDIR}/ALLchr.SNPs.pruned
plink \
--bfile ${TMPDIR}/ALLchr.SNPs.pruned \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 21 \
--out ${TMPDIR}/ALLchr.SNPs.pruned
plink \
--bfile ${TMPDIR}/ALLchr.SNPs.pruned \
--recode A \
--reference-allele $reference_alleles_ALLchr \
--chr-set 21 \
--out ${TMPDIR}/ALLchr.SNPs.pruned

echo "copying all chromosome pruned plink files at $(date)"
cp ${TMPDIR}/ALLchr.SNPs.prune* --target-directory=${STITCH_PLINK_dir}

#then take the prune files to filter the STITCH_vcf output
LD_pruned_bed=$TMPDIR/LD_pruned.bed
awk -v FS=":" -v OFS='\t' 'BEGIN{print "chrom", "chromStart", "chromEnd"}{print $1, $2, $2}' ${TMPDIR}/ALLchr.SNPs.prune.in > $LD_pruned_bed

for index in {1..21}; do
get_wheat_chr_genome $index
input_VCF=${STITCH_output_dir}/chr${chr}${genome}/stitch.${chr}${genome}.vcf.gz
echo "$index vcftools for chr${chr}${genome} at $(date)"
vcftools --gzvcf $input_VCF \
--recode --recode-INFO-all \
--bed $LD_pruned_bed \
--stdout | \
gzip -c > ${TMPDIR}/stitch.chr${chr}${genome}.vcf.gz
done

mkdir -p ${STITCH_haplotype_dir}/LD_pruned_vcfs
cp ${TMPDIR}/stitch.chr*.vcf.gz --target-directory=${STITCH_haplotype_dir}/LD_pruned_vcfs


