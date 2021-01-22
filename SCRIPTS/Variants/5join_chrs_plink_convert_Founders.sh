#!/bin/bash -l
#$ -cwd
#$ -N STEP4_join_chrs_plink_convert_Founders
#$ -o STEP4_join_chrs_plink_convert_Founders.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=50G

set -a 
source ../Input/data_and_variables.sh
source ../Input/software.sh

samplename=Founders_combined
mkdir -p ${founder_vcf_dir}/biallelic_polymorphic_filter_merged #output dir

for index in {1..43}; do

#variables for converting position
line=$(sed -n ${index}p $chr_names_lengths_file) 
chromosome=$(echo $line | awk '{print $2}')
start_position=$(echo $chromosome | awk -v FS=":" '{print $2}' | awk -v FS="-" '{print $1}')
start_position_add=$(( $start_position - 1 ))
chr_name=$(echo $chromosome | awk -v FS=":" '{print $1}')
chr_name_half=$(echo $line | awk '{print $1}') 
echo $chr_name_half

input_vcf=${founder_vcf_dir}/raw_vcfs/${samplename}_${chr_name_half}*.vcf
#add the start position to the vcfs
awk -v OFS='\t' -v chromosome=$chromosome -v start_position_add=$start_position_add -v chr_name=$chr_name '{ if ($1 ~ chromosome && $1 !~ /#/) {$1=chr_name; $2=$2+start_position_add;} print $0 }' $input_vcf | \
gzip -c > ${TMPDIR}/${chr_name_half}.raw.vcf.gz
#for the second half of each chromosome, this will concat both halves together
vcf-concat ${TMPDIR}/chr${chr_name}*.raw.vcf.gz | gzip -c > ${founder_vcf_dir}/raw_vcfs_merged/chr${chr_name}.vcf.gz

#repeat, but after filtering 
input_vcf=${founder_vcf_dir}/biallelic_polymorphic_filter/${samplename}_${chr_name_half}*.vcf
#add the start position to the vcfs
awk -v OFS='\t' -v chromosome=$chromosome -v start_position_add=$start_position_add -v chr_name=$chr_name '{ if ($1 ~ chromosome && $1 !~ /#/) {$1=chr_name; $2=$2+start_position_add;} print $0 }' $input_vcf | \
gzip -c > ${TMPDIR}/${chr_name_half}.all.vcf.gz
#for the second half of each chromosome, this will concat both halves together
vcf-concat ${TMPDIR}/chr${chr_name}*.all.vcf.gz | gzip -c > ${founder_vcf_dir}/biallelic_polymorphic_filter_merged/chr${chr_name}.vcf.gz

#repeat, including only polymorphic SNPs without het or missing calls
input_vcf=${founder_vcf_dir}/nohet_biallelic_polymorphic_filter/${samplename}_${chr_name_half}*.vcf
mkdir -p ${founder_vcf_dir}/nohet_biallelic_polymorphic_filter_merged #output dir
#add the start position to the vcfs
awk -v OFS='\t' -v chromosome=$chromosome -v start_position_add=$start_position_add -v chr_name=$chr_name '{ if ($1 ~ chromosome && $1 !~ /#/) {$1=chr_name; $2=$2+start_position_add;} print $0 }' $input_vcf | \
gzip -c > ${TMPDIR}/${chr_name_half}.nohet.vcf.gz
vcf-concat ${TMPDIR}/chr${chr_name}*.nohet.vcf.gz | gzip -c > ${founder_vcf_dir}/nohet_biallelic_polymorphic_filter_merged/chr${chr_name}.vcf.gz

done

### then convert all this to PLINK format
mkdir -p ${founder_plink_dir_merged}/alternate_allele_maps
###first create the list of alternate/reference alleles
for index in {1..22}; do
get_wheat_chr_genome $index
zcat ${founder_vcf_dir}/biallelic_polymorphic_filter_merged/chr${chr}${genome}.vcf.gz | \
awk -v OFS='\t' '/^[^#]/ {print $1":"$2,$5}' > ${founder_plink_dir_merged}/alternate_allele_maps/${samplename}_chr${chr}${genome}.txt
done
cat ${founder_plink_dir_merged}/alternate_allele_maps/${samplename}_chr*.txt > ${founder_plink_dir_merged}/alternate_allele_maps/${samplename}_ALLchr.txt

###convert to PLINK format
chromosome_map=${TMPDIR}/chromosome_map.txt
awk -v FS='\t' '{print $1}' ${ref}.fai | awk -v FS=":" '{print $1}' | sort -u | awk -v OFS='\t' '{print $1, NR}' > $chromosome_map
reference_alleles_ALLchr=${founder_plink_dir_merged}/alternate_allele_maps/${samplename}_ALLchr.txt
for index in {1..22}; do
get_wheat_chr_genome $index
#specify reference allele for plink (based on vcf's)
reference_alleles=${founder_plink_dir_merged}/alternate_allele_maps/${samplename}_chr${chr}${genome}.txt
input_VCF=${founder_vcf_dir}/raw_vcfs_merged/chr${chr}${genome}.vcf.gz
vcftools --gzvcf $input_VCF \
--plink --chrom-map ${chromosome_map} \
--remove-indels \
--min-alleles 2 --max-alleles 2 \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--max-missing ${max_missing_param} \
--out ${TMPDIR}/chr${chr}${genome}.raw
#repeat, but after filtering 
input_VCF=${founder_vcf_dir}/biallelic_polymorphic_filter_merged/chr${chr}${genome}.vcf.gz
vcftools --gzvcf $input_VCF \
--plink --chrom-map ${chromosome_map} \
--remove-indels \
--min-alleles 2 --max-alleles 2 \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--max-missing ${max_missing_param} \
--out ${TMPDIR}/chr${chr}${genome}.all
#repeat but use version with no heterozygotes 
input_VCF_nohet=${founder_vcf_dir}/nohet_biallelic_polymorphic_filter_merged/chr${chr}${genome}.vcf.gz
vcftools --gzvcf ${input_VCF_nohet} \
--plink --chrom-map ${chromosome_map} \
--remove-indels \
--min-alleles 2 --max-alleles 2 \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--max-missing ${max_missing_param} \
--out ${TMPDIR}/chr${chr}${genome}.nohet
done

#merge all the plink files. Plink syntax requires one file and then a list of other files.
pedfiles=$(ls ${TMPDIR}/chr*.raw.ped)
mapfiles=$(ls ${TMPDIR}/chr*.raw.map)
paste <(echo "$pedfiles") <(echo "$mapfiles") --delimiters '\t' > ${TMPDIR}/plink_files_to_merge.1.txt
first_file=$(sed "1q;d" ${TMPDIR}/plink_files_to_merge.1.txt | awk '{print $1}')
first_file_prefix="${first_file%.ped*}"
tail -n+2 ${TMPDIR}/plink_files_to_merge.1.txt > ${TMPDIR}/plink_files_to_merge.2.txt
plink \
--file $first_file_prefix \
--merge-list ${TMPDIR}/plink_files_to_merge.2.txt \
--recode \
--chr-set 43 \
--out ${TMPDIR}/ALLchr.raw
plink \
--file ${TMPDIR}/ALLchr \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 43 \
--out ${TMPDIR}/ALLchr.raw

#merge all the plink files. Plink syntax requires one file and then a list of other files.
pedfiles=$(ls ${TMPDIR}/chr*.all.ped)
mapfiles=$(ls ${TMPDIR}/chr*.all.map)
paste <(echo "$pedfiles") <(echo "$mapfiles") --delimiters '\t' > ${TMPDIR}/plink_files_to_merge.1.txt
first_file=$(sed "1q;d" ${TMPDIR}/plink_files_to_merge.1.txt | awk '{print $1}')
first_file_prefix="${first_file%.ped*}"
tail -n+2 ${TMPDIR}/plink_files_to_merge.1.txt > ${TMPDIR}/plink_files_to_merge.2.txt
plink \
--file $first_file_prefix \
--merge-list ${TMPDIR}/plink_files_to_merge.2.txt \
--recode \
--chr-set 43 \
--out ${TMPDIR}/ALLchr
plink \
--file ${TMPDIR}/ALLchr \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 43 \
--out ${TMPDIR}/ALLchr

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
--chr-set 43 \
--out ${TMPDIR}/ALLchr.nohet
plink \
--file ${TMPDIR}/ALLchr.nohet \
--recode A-transpose \
--reference-allele $reference_alleles_ALLchr \
--chr-set 43 \
--out ${TMPDIR}/ALLchr.nohet

echo "copying all chromosome plink files for ${samplename} at $(date)"
cp ${TMPDIR}/ALLchr.* --target-directory=${founder_plink_dir_merged}

