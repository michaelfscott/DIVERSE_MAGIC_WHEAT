#!/bin/bash -l
#$ -cwd 
#$ -N genotype_Balfourier_founders
#$ -o genotype_Balfourier_founders.log
#$ -l mem=20G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -t 1-43
#$ -l tmpfs=20G

set -a
source ../Input/data_and_variables.sh
source ../Input/software.sh

#haplotypeCaller options
java_memory_param=12
verbosity_param=ERROR
interval_padding_param=100
minimum_mapping_quality_param=30

get_wheat_chr_genome_half $SGE_TASK_ID
echo chr${chr}${genome}${half}

known_variants=${TMPDIR}/chr${chr}${genome}${half}.bed

if [ "$half" == _2 ] ; then
  line=$(grep ${chr}${genome} ${ref}.fai | head -n 2 | tail -n 1)
  chr_out=$(echo $line | awk '{print $1}')
  length_half=$(grep ${chr}${genome} ${ref}.fai | head -n 1 | awk '{print $2}')
  awk -F $'\t' \
    -v match_chr="chr$chr$genome" \
    -v length_out="$length_half" \
    -v chr_out="$chr_out" \
    '$3~match_chr && $4>length_out {OFS=FS; $4-=length_out; print chr_out, ($4 - 1), $4}' $Balfourier_known_variants \
    > ${known_variants}.tmp
else 
  line=$(grep ${chr}${genome} ${ref}.fai | head -n 1 )
  chr_out=$(echo $line | awk '{print $1}')
  length_half=$(echo $line | awk '{print $2}')
  awk -F $'\t' -v match_chr="chr$chr$genome" -v length_out="$length_half" -v chr_out="$chr_out" '$3~match_chr && $4<=length_out {OFS=FS; print chr_out, ($4 - 1), $4}' $Balfourier_known_variants > ${known_variants}.tmp
fi

module load bedtools
bedtools sort -i ${known_variants}.tmp > ${known_variants}

#log and directory for temporary files
log=${project_dir}/scripts/2Variants_genotype_diversity_populations/log_genotype_Balfourier
tmpdir=${TMPDIR}/tmp_job$SGE_TASK_ID
mkdir -p ${tmpdir}
mkdir -p ${log}
mkdir -p $founder_vcf_dir_Balfourier
mkdir -p ${founder_vcf_dir_Balfourier}/chr${chr}${genome}${half}

samplename=Founders_combined

#must first index the vcf at which variants will be called
gatk IndexFeatureFile --java-options "-Xmx${java_memory_param}G" \
--feature-file $known_variants

# in order to call reference alleles, must genotype each founder individually, then combine
for i in {1..16}; do
samp=$(get_founder_samplename $i)
bamlist_file=${TMPDIR}/${samp}_bamlist.txt
if [ $samp == "Holdfast" ]; then  
ls ${founder_alignments_dir}/Holdfast/*.dedup.bam > ${bamlist_file}
else 
ls ${founder_alignments_dir}/${samp}*/${samp}*.dedup.2.bam > ${bamlist_file}
fi

input_files_arguments=${TMPDIR}/${samp}.arguments.txt
while read line ; do
echo "--input $line"
done < ${bamlist_file} > ${input_files_arguments}

gatk HaplotypeCaller --java-options "-Xmx${java_memory_param}G" \
--reference ${ref} \
--TMP_DIR=${tmpdir} \
--arguments_file $input_files_arguments \
--intervals ${known_variants} \
--minimum-mapping-quality ${minimum_mapping_quality_param} \
--emit-ref-confidence BP_RESOLUTION \
--verbosity $verbosity_param \
--output-mode EMIT_ALL_SITES \
--standard-min-confidence-threshold-for-calling 0 \
--output ${TMPDIR}/${samp}_chr${chr}${genome}${half}.vcf 

bgzip -c ${TMPDIR}/${samp}_chr${chr}${genome}${half}.vcf > ${founder_vcf_dir_Balfourier}/chr${chr}${genome}${half}/${samp}_chr${chr}${genome}${half}.vcf.gz
bcftools index ${founder_vcf_dir_Balfourier}/chr${chr}${genome}${half}/${samp}_chr${chr}${genome}${half}.vcf.gz

done

#combine founders
${HOME}/packages/bcftools/bcftools merge --gvcf ${ref} ${founder_vcf_dir_Balfourier}/chr${chr}${genome}${half}/*_chr${chr}${genome}${half}.vcf.gz \
--output ${founder_vcf_dir_Balfourier}/${samplename}_chr${chr}${genome}${half}.vcf


