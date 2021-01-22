#!/bin/bash -l
#$ -cwd
#$ -N STEP3_filter_vcfs_PLINK_convert
#$ -o STEP3_filter_vcfs_PLINK_convert.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=50G

set -a 
source ../Input/data_and_variables.sh
source ../Input/software.sh

#for converting to plink correctly
infoScoreExclude=${STITCH_output_dir}/Exclude_Info_${infoScoreFilter}_allChr.txt
chromosome_map=${TMPDIR}/chromosome_map.txt
awk -v FS='\t' '{print $1}' ${ref}.fai | awk -v FS=":" '{print $1}' | sort -u | awk -v OFS='\t' '{print $1, NR}' > $chromosome_map
reference_alleles_ALLchr=${founder_plink_dir_merged}/alternate_allele_maps/Founders_combined_ALLchr.txt

#first filter SNPs with low info score
for index in {1..21}; do

get_wheat_chr_genome $index

input_VCF=${STITCH_output_dir}/chr${chr}${genome}/stitch.${chr}${genome}.vcf.gz
low_info_SNPs=${STITCH_output_dir}/chr${chr}${genome}/Exclude_Info_${infoScoreFilter}.txt
exclude_positions=${TMPDIR}/exclude_positions.txt
awk -v FS=":" -v OFS="\t" '{print $1, $2}' $low_info_SNPs > ${exclude_positions}

#apply filters to VCF

echo "$index vcftools for chr${chr}${genome} at $(date)"
vcftools --gzvcf $input_VCF \
--max-missing ${max_missing_param} \
--maf ${minMAF_param} \
--exclude-positions $exclude_positions \
--recode --recode-INFO-all \
--stdout | \
gzip -c > ${STITCH_vcfs_dir}/chr${chr}${genome}.vcf.gz

echo "$index vcftools for chr${chr}${genome} at $(date)"
vcftools --gzvcf $input_VCF \
--plink --chrom-map ${chromosome_map} \
--max-missing ${max_missing_param} \
--exclude-positions $exclude_positions \
--maf ${minMAF_param} \
--out ${TMPDIR}/chr${chr}${genome}.SNPs

done

#merge files, PLINK input requires that one file is specified and then a list of other plink files to merge is supplied.
pedfiles=$(ls ${TMPDIR}/chr*.SNPs.ped)
mapfiles=$(ls ${TMPDIR}/chr*.SNPs.map)
paste <(echo "$pedfiles") <(echo "$mapfiles") --delimiters '\t' > ${TMPDIR}/plink_files_to_merge.1.txt
first_file=$(sed "1q;d" ${TMPDIR}/plink_files_to_merge.1.txt | awk '{print $1}')
first_file_prefix="${first_file%.ped*}"
tail -n+2 ${TMPDIR}/plink_files_to_merge.1.txt > ${TMPDIR}/plink_files_to_merge.2.txt
plink \
--file $first_file_prefix \
--merge-list ${TMPDIR}/plink_files_to_merge.2.txt \
--recode \
--exclude ${infoScoreExclude} \
--chr-set 43 \
--out ${TMPDIR}/ALLchr.SNPs
plink \
--file ${TMPDIR}/ALLchr.SNPs \
--recode A-transpose \
--exclude ${infoScoreExclude} \
--reference-allele $reference_alleles_ALLchr \
--chr-set 43 \
--out ${TMPDIR}/ALLchr.SNPs

echo "copying all chromosome plink files at $(date)"
mkdir -p ${STITCH_PLINK_dir}
cp ${TMPDIR}/ALLchr.SNPs* --target-directory=${STITCH_PLINK_dir}

