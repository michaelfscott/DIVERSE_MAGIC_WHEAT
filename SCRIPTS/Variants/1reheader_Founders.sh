#!/bin/bash -l
#$ -cwd 
#$ -N STEP0_reheader_Founders
#$ -o STEP0_reheader_Founders.log
#$ -l mem=10G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -t 1-16
#$ -l tmpfs=100G

set -a
source ../Input/data_and_variables.sh
source ../Input/software.sh

founder_index=$SGE_TASK_ID

#get samplename for this job
samplename=$(get_founder_samplename $founder_index)
echo "$SGE_TASK_ID founder ${founder_index} $samplename"

#only one input file for Holdfast but there is are two input bams for other founders (one from each capture)
if [[ $samplename == "Holdfast" ]]; then

echo "Holdfast"

else

input_bam1=$(ls ${founder_alignments_dir}/${samplename}_L1/*.dedup.bam | sed "1q;d")
input_bam2=$(ls ${founder_alignments_dir}/${samplename}_L2/*.dedup.bam | sed "1q;d")

samtools view -H $input_bam1 > ${TMPDIR}/header1.txt
sed 's/SM:'"$samplename"'_L1/SM:'"$samplename"''/g ${TMPDIR}/header1.txt > ${TMPDIR}/header1.edit.txt
samtools reheader ${TMPDIR}/header1.edit.txt $input_bam1 > ${founder_alignments_dir}/${samplename}_L1/${samplename}_L1.dedup.2.bam
samtools index ${founder_alignments_dir}/${samplename}_L1/${samplename}_L1.dedup.2.bam

samtools view -H $input_bam2 > ${TMPDIR}/header2.txt
sed 's/SM:'"$samplename"'_L2/SM:'"$samplename"''/g ${TMPDIR}/header2.txt > ${TMPDIR}/header2.edit.txt
samtools reheader ${TMPDIR}/header2.edit.txt $input_bam2 > ${founder_alignments_dir}/${samplename}_L2/${samplename}_L2.dedup.2.bam
samtools index ${founder_alignments_dir}/${samplename}_L2/${samplename}_L2.dedup.2.bam

fi

