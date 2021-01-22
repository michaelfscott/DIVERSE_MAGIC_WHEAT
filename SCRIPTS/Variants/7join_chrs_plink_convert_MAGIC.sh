#!/bin/bash -l
#$ -cwd
#$ -N STEP5_join_chrs_plink_convert_MAGIC
#$ -o STEP5_join_chrs_plink_convert_MAGIC.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=50G

set -a 
source ../Input/data_and_variables.sh
source ../Input/software.sh

samplename=MAGIC_combined
mkdir -p ${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter_merged #output dir

for index in {1..43}; do

#variables for converting position
line=$(sed -n ${index}p $chr_names_lengths_file) 
chromosome=$(echo $line | awk '{print $2}')
start_position=$(echo $chromosome | awk -v FS=":" '{print $2}' | awk -v FS="-" '{print $1}')
start_position_add=$(( $start_position - 1 ))
chr_name=$(echo $chromosome | awk -v FS=":" '{print $1}')
chr_name_half=$(echo $line | awk '{print $1}') 

#including only polymorphic SNPs without het or missing calls
input_vcf=${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter/${samplename}_${chr_name_half}*.vcf
mkdir -p ${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter_merged #output dir
#add the start position to the vcfs
awk -v OFS='\t' -v chromosome=$chromosome -v start_position_add=$start_position_add -v chr_name=$chr_name '{ if ($1 ~ chromosome && $1 !~ /#/) {$1=chr_name; $2=$2+start_position_add;} print $0 }' $input_vcf | \
gzip -c > ${TMPDIR}/${chr_name_half}.nohet.vcf.gz
vcf-concat ${TMPDIR}/chr${chr_name}*.nohet.vcf.gz | gzip -c > ${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter_merged/chr${chr_name}.vcf.gz

done

###convert to PLINK format
chromosome_map=${TMPDIR}/chromosome_map.txt
awk -v FS='\t' '{print $1}' ${ref}.fai | awk -v FS=":" '{print $1}' | sort -u | awk -v OFS='\t' '{print $1, NR}' > $chromosome_map
#specify reference allele for plink (based on vcf's)
reference_alleles_ALLchr=${founder_plink_dir_merged}/alternate_allele_maps/Founders_combined_ALLchr.txt

for index in {1..22}; do
get_wheat_chr_genome $index
input_VCF_nohet=${MAGIC_vcf_dir}/nohet_biallelic_polymorphic_filter_merged/chr${chr}${genome}.vcf.gz
vcftools --gzvcf ${input_VCF_nohet} \
--plink --chrom-map ${chromosome_map} \
--out ${TMPDIR}/chr${chr}${genome}.nohet
#repeat with only depth 4 calls
vcftools --gzvcf ${input_VCF_nohet} \
--plink --chrom-map ${chromosome_map} \
--minDP 4 \
--out ${TMPDIR}/chr${chr}${genome}.nohet.DP4
done

#merge all the plink files. Plink syntax requires one file and then a list of other files.
pedfiles=$(ls ${TMPDIR}/chr*.nohet.ped)
mapfiles=$(ls ${TMPDIR}/chr*.nohet.map)
paste <(echo "$pedfiles") <(echo "$mapfiles") --delimiters '\t' > ${TMPDIR}/plink_files_to_merge.1.txt
first_file=$(sed "1q;d" ${TMPDIR}/plink_files_to_merge.1.txt | awk '{print $1}')
first_file_prefix="${first_file%.ped*}"
tail -n+2 ${TMPDIR}/plink_files_to_merge.1.txt > ${TMPDIR}/plink_files_to_merge.2.txt
plink \
--file $first_file_prefix \
--merge-list ${TMPDIR}/plink_files_to_merge.2.txt \
--recode \
--reference-allele $reference_alleles_ALLchr \
--chr-set 22 \
--out ${TMPDIR}/ALLchr
plink \
--file ${TMPDIR}/ALLchr \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 22 \
--out ${TMPDIR}/ALLchr

echo "copying all chromosome plink files for ${samplename} at $(date)"
mkdir -p $MAGIC_plink_dir_merged
cp ${TMPDIR}/ALLchr.* --target-directory=${MAGIC_plink_dir_merged}

###repeat for the calls with at least DP4
#merge all the plink files. Plink syntax requires one file and then a list of other files.
pedfiles=$(ls ${TMPDIR}/chr*.nohet.DP4.ped)
mapfiles=$(ls ${TMPDIR}/chr*.nohet.DP4.map)
paste <(echo "$pedfiles") <(echo "$mapfiles") --delimiters '\t' > ${TMPDIR}/plink_files_to_merge.1.txt
first_file=$(sed "1q;d" ${TMPDIR}/plink_files_to_merge.1.txt | awk '{print $1}')
first_file_prefix="${first_file%.ped*}"
tail -n+2 ${TMPDIR}/plink_files_to_merge.1.txt > ${TMPDIR}/plink_files_to_merge.2.txt
plink \
--file $first_file_prefix \
--merge-list ${TMPDIR}/plink_files_to_merge.2.txt \
--recode \
--reference-allele $reference_alleles_ALLchr \
--chr-set 22 \
--out ${TMPDIR}/ALLchr.DP4
plink \
--file ${TMPDIR}/ALLchr.DP4 \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 22 \
--out ${TMPDIR}/ALLchr.DP4

echo "copying all chromosome plink files for ${samplename} at $(date)"
mkdir -p $MAGIC_plink_dir_merged
cp ${TMPDIR}/ALLchr.DP4* --target-directory=${MAGIC_plink_dir_merged}

