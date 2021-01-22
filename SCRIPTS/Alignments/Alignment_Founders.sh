#$ -N Align_Founders
#$ -o Align_Founders.log
#$ -cwd
#$ -l mem=30G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=30G
#$ -t 1-572

source ../Input/data_and_variables.sh
source ../Input/software.sh

#parameters 
numthreads=1

#get all the fastqs, trim off extension that includes read pair, sort uniq to remove duplicate
for file in ${founder_split_fastq_dir}/*/*.fastq.gz ; do echo "${file%_R[1,2].fastq.gz}"; done | sort -u > ${TMPDIR}/split_fastq_list.txt

#select from sample_list the fastq file to work with
input_fastq_prefix=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/split_fastq_list.txt)
#get the basename
experimentname=$(basename $input_fastq_prefix) #here, also split by sub-division of fastq file.
#trim the trailing part of the filename, including the lane and experiment number
samplename=$(echo "${experimentname%_[B,P][R,D][O,S]*}" | cut -c 6-) #two possible experiment extensions in these file names. 

fastq_1=${input_fastq_prefix}_R1.fastq.gz
fastq_2=${input_fastq_prefix}_R2.fastq.gz

### Get Read Group
library=LIB-${samplename}
string=$(gzip -cd $fastq_1 | head -n 2 | grep "^@") #get the read name from the head of fastq1
instrument_tmp=$(echo $string | cut -d":" -f 1) #first field is a unique instrument name
instrument=${instrument_tmp#"@"} #trim off leading @ from the fastq read string
flowcell=$(echo $string | cut -d":" -f 3) #flowcell ID
barcode_tmp=$(echo $string | cut -d":" -f 10) #the barcode is the last part of this string
barcode=$(echo $barcode_tmp | cut -d"+" -f 1) #the barcode has a +, I remove this because I don't want to introduce special characters
lane=$(echo $string | cut -d":" -f 4) #the lane is the fourth part.
samplename_for_RG=${samplename%_L[1-2]}

runRG="@RG\tID:${flowcell}.${lane}\tLB:${library}\tPL:ILLUMINA\tSM:${samplename_for_RG}\tPU:${flowcell}.${lane}.${barcode}"

mkdir -p ${TMPDIR}/${experimentname} #specific tmpdir for each sample

echo "$SGE_TASK_ID Aligning $SGE_TASK_ID ${experimentname} at $(date) on $(hostname)"
bwa mem -M -R $runRG -t $numthreads $ref $fastq_1 $fastq_2 2> $log/${experimentname}.bwamem.log | \
samtools sort -@ $numthreads -T ${TMPDIR}/${experimentname}/${experimentname} -o ${TMPDIR}/${experimentname}.sorted.bam - 2> ${log}/${experimentname}.sort.log #sort with specified tmp directory and output to .sorted.bam. 
echo "$SGE_TASK_ID indexing ${experimentname} before dedup at $(date)"
samtools index ${TMPDIR}/${experimentname}.sorted.bam

echo "$SGE_TASK_ID copying before dedup at $(date)"
mkdir -p ${founder_alignments_dir}/${samplename}
cp ${TMPDIR}/${experimentname}.sorted.bam* --target-directory=${founder_alignments_dir}/${samplename}/

