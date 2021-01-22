#here we specify various standard bioinformatic tools for analysis of NGS data
module load bwa/0.7.12/gnu-4.9.2
module load samtools/1.3.1/gnu-4.9.2
module load java/1.8.0_92 gatk/4.0.8.0
module load perl/5.22.0
module load vcftools/0.1.15/gnu-4.9.2
module load bcftools/2.1/gnu-4.9.2
module load plink/1.90b3.40

#we use the fastp toolkit to split paired-end fastq files
PATH=/home/ucbtmf1/packages/fastp_v0.19.7/:$PATH
PATH=/home/ucbtmf1/packages/htslib/:$PATH
#FaSTLMM
PATH=/home/ucbtmf1/packages/FaSTLMM.207c.Linux/Linux_MKL/:$PATH
