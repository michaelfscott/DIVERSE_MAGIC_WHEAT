#$ -N Split_Founder_Fastq
#$ -o Split_Founder_Fastq.log
#$ -cwd
#$ -l mem=10G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=30G
#$ -t 1-32

source ../Input/data_and_variables.sh
source ../Input/software.sh

#get all the fastqs, trim off extension that includes read pair, sort uniq to remove duplicate
for file in ${raw_Founder_fastq_dir}/*/*.fastq.gz ; do echo "${file%_R[1,2].fastq.gz}"; done | sort -u > ${TMPDIR}/experiment_list.txt

#select from sample_list the fastq file to work with
input_fastq_prefix=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/experiment_list.txt)
#get the basename
experimentname=$(basename $input_fastq_prefix)
#trim the trailing part of the filename, including the lane and experiment number
samplename=$(echo "${experimentname%_PRO183*}")

numthreads=1
num_lines_per_split=40000000
sample_founder_split_fastq_dir=${founder_split_fastq_dir}/${samplename}/
mkdir -p $sample_founder_split_fastq_dir

#split the fastq files using fastp 
fastp \
--in1 ${input_fastq_prefix}_R1.fastq.gz \
--in2 ${input_fastq_prefix}_R2.fastq.gz \
--out1 ${sample_founder_split_fastq_dir}/${experimentname}_R1.fastq.gz \
--out2 ${sample_founder_split_fastq_dir}/${experimentname}_R2.fastq.gz \
--thread $numthreads \
--split_by_lines $num_lines_per_split \
--disable_trim_poly_g \
--disable_adapter_trimming \
--disable_quality_filtering \
--disable_length_filtering \
--split_prefix_digits 4


