#$ -cwd
#$ -N STEP3_Merge_Dedup_Founders
#$ -o STEP3_Merge_Dedup_Founders.log
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=500G
#$ -t 1-31

#parameters
numthreads=1
java_memory_param=18

source ../Input/data_and_variables.sh
source ../Input/software.sh

#for alignment statistics
mkdir -p ${founder_alignments_dir}/stats
mkdir -p ${founder_alignments_dir}/dedup_metrics

#match only the directories for founder data (either _L[1,2] or Holdfast
datadir=${TMPDIR}/datadirs.txt
ls -d ${founder_alignments_dir}/*{_L,Holdfast}* > $datadir
input_datadir=$(sed "${SGE_TASK_ID}q;d" $datadir )
samplename=$(echo $input_datadir | awk -F'/' '{print $(NF)}')

#make a list of bams to merge
bamlist=${TMPDIR}/bams.txt
ls ${input_datadir}/*.sorted.bam > $bamlist

#the c and p options are so that RG and PG tags with colliding IDs are merged rather than adding a suffix to differentiate them. 
echo "$SGE_TASK_ID merging ${samplename} at $(date)"
samtools merge -c -p -b $bamlist ${TMPDIR}/${samplename}.bam 
echo "$SGE_TASK_ID indexing ${samplename} at $(date)"
samtools index ${TMPDIR}/${samplename}.bam

#copy out alignment file in case dedup fails
echo "$SGE_TASK_ID copying at $(date)"
mkdir -p ${founder_alignments_dir}/${samplename}
cp ${TMPDIR}/${samplename}.bam* --target-directory=${founder_alignments_dir}/${samplename}/

#get alignment statistics before marking duplicates (and filtering unmapped reads)
echo "$SGE_TASK_ID alignment statistics ${samplename} at $(date)"
sample_stats_dir=${founder_alignments_dir}/stats/${samplename}
mkdir -p $sample_stats_dir
samtools stats ${TMPDIR}/${samplename}.bam > ${sample_stats_dir}/${samplename}.raw.stats
grep "^COV" ${sample_stats_dir}/${samplename}.raw.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.raw.cov.txt
samtools stats --target-regions $gene_bed ${TMPDIR}/${samplename}.bam > ${sample_stats_dir}/${samplename}.raw.gene.stats
grep "^COV" ${sample_stats_dir}/${samplename}.raw.gene.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.raw.gene.cov.txt
samtools stats --target-regions $prom_and_5UTR_bed ${TMPDIR}/${samplename}.bam > ${sample_stats_dir}/${samplename}.raw.prom_and_5UTR.stats
grep "^COV" ${sample_stats_dir}/${samplename}.raw.prom_and_5UTR.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.raw.prom_and_5UTR.cov.txt

#repeat but filter for mapQ30
min_mapQ=30
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.bam | samtools stats > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.stats
grep "^COV" ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.cov.txt
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.bam | samtools stats --target-regions $gene_bed > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.gene.stats
grep "^COV" ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.gene.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.gene.cov.txt
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.bam | samtools stats --target-regions $prom_and_5UTR_bed > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.prom_and_5UTR.stats
grep "^COV" ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.prom_and_5UTR.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.mapQ${min_mapQ}.prom_and_5UTR.cov.txt

#for picard markDuplicates, need to remove unmapped reads first. 
echo "filtering ${samplename} to mapped reads only at $(date)"
samtools view -b -h -F 4 ${TMPDIR}/${samplename}.bam > ${TMPDIR}/${samplename}.mapped.bam
samtools index ${TMPDIR}/${samplename}.mapped.bam

#make a temporary directory for dedup
mkdir -p ${TMPDIR}/${samplename}
#output dedup metrics to file in this dir
dedup_metrics_dir=${founder_alignments_dir}/dedup_metrics/${samplename}
mkdir -p ${dedup_metrics_dir}

#mark Duplicates and index
echo "$SGE_TASK_ID mark Duplicates for ${samplename} at $(date)"
gatk MarkDuplicates --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=${TMPDIR}/${samplename} \
--INPUT ${TMPDIR}/${samplename}.mapped.bam \
--METRICS_FILE ${dedup_metrics_dir}/${samplename}.dedup_metrics.txt \
--OUTPUT ${TMPDIR}/${samplename}.dedup.bam 
echo "$SGE_TASK_ID indexing merged dedup bam ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${samplename}.dedup.bam 

#copy out dedup file
echo "$SGE_TASK_ID copying at $(date)"
mkdir -p ${founder_alignments_dir}/${samplename}
cp ${TMPDIR}/${samplename}.dedup.bam* --target-directory=${founder_alignments_dir}/${samplename}/

#get alignment statistics 
echo "$SGE_TASK_ID alignment statistics ${samplename} at $(date)"
sample_stats_dir=${founder_alignments_dir}/stats/${samplename}
mkdir -p $sample_stats_dir
samtools stats ${TMPDIR}/${samplename}.dedup.bam > ${sample_stats_dir}/${samplename}.dedup.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.cov.txt
samtools stats --target-regions $gene_bed ${TMPDIR}/${samplename}.dedup.bam > ${sample_stats_dir}/${samplename}.dedup.gene.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.gene.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.gene.cov.txt
samtools stats --target-regions $prom_and_5UTR_bed ${TMPDIR}/${samplename}.dedup.bam > ${sample_stats_dir}/${samplename}.dedup.prom_and_5UTR.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.prom_and_5UTR.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.prom_and_5UTR.cov.txt

#repeat but filter for mapQ30
min_mapQ=30
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.dedup.bam | samtools stats > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.cov.txt
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.dedup.bam | samtools stats --target-regions $gene_bed > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.gene.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.gene.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.gene.cov.txt
samtools view -h -b -q $min_mapQ ${TMPDIR}/${samplename}.dedup.bam | samtools stats --target-regions $prom_and_5UTR_bed > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.prom_and_5UTR.stats
grep "^COV" ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.prom_and_5UTR.stats | cut -f 2- > ${sample_stats_dir}/${samplename}.dedup.mapQ${min_mapQ}.prom_and_5UTR.cov.txt


