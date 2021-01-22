#all data will be stored in the project dir
project_dir=/lustre/projects/MAGIC_WHEAT/
#log files
log=/home/ucbtmf1/Scratch/MAGIC_WHEAT/log

### fastqs
#the founder fastq dir
raw_Founder_fastq_dir=${project_dir}/1raw_Founder_fastq
#the MAGIC fastq dir
raw_MAGIC_fastq_dir=${project_dir}/1raw_MAGIC_fastq
new_MAGIC_fastqs=${raw_MAGIC_fastq_dir}/raw_data2

### Reference genome
reference_genome_dir=${project_dir}/2Reference_Genomes/triticum_aestivum_split
#the reference genome has been split to comply with SAM format specifications
ref=${reference_genome_dir}/Triticum_aestivum.IWGSC.dna.toplevel.split.fasta
#the gene and promoter capture bed files have also been modified to correspond to the split reference genome
gene_bed=${reference_genome_dir}/gene_capture_targets_HC.bed
prom_and_5UTR_bed=${reference_genome_dir}/promoter_capture_targets_HC_and_5UTR.edit.bed
gene_and_prom_and_5UTR_bed=${reference_genome_dir}/gene_and_promoter_targets_HC_and_5UTR.edit.bed
gff_features_file=${project_dir}/2Reference_Genomes/triticum_aestivum_gff/Triticum_aestivum.IWGSC.42.gff3
#known_variants
Balfourier_known_variants=${project_dir}/2Reference_Genomes/Balfourier_et_al_Array/Balfourier_et_al_Wheat_Phylogeography_DataS2.tab
He_known_variants=${project_dir}/2Reference_Genomes/He_et_al_Exomes/all.GP08_mm75_het3_publication01142019.vcf

### Alignments 
alignments_dir=${project_dir}/3Alignments
#where the founder files will be split into
founder_split_fastq_dir=${alignments_dir}/1Founder_split_fastq
#where the aligned founder files will go
founder_alignments_dir=${alignments_dir}/2Founder_Alignments
#where the MAGIC alignments are
MAGIC_split_fastq_dir=${alignments_dir}/1MAGIC_split_fastq
MAGIC_alignments_dir_before_dedup=${alignments_dir}/Alignments_before_dedup
MAGIC_alignments_dir=${alignments_dir}/Alignments

### Variant Calls
variants_dir=${project_dir}/4Variants
founder_gvcfs_dir=${variants_dir}/gvcfs_Founders
founder_vcf_dir=${variants_dir}/vcfs_Founders
founder_plink_dir=${variants_dir}/plink_Founders
founder_plink_dir_merged=${variants_dir}/plink_Founders_merged
founder_reference_panel=${variants_dir}/reference_haplotypes
founder_SNPlist_dir=${variants_dir}/SNPlists_Founders
MAGIC_vcf_dir=${variants_dir}/vcfs_MAGIC
MAGIC_vcf_dir_Balfourier=${variants_dir}/vcfs_MAGIC_Balfourier
founder_vcf_dir_Balfourier=${variants_dir}/vcfs_Founders_Balfourier
MAGIC_vcf_dir_He=${variants_dir}/vcfs_MAGIC_He
founder_vcf_dir_He=${variants_dir}/vcfs_Founders_He
MAGIC_plink_dir=${variants_dir}/plink_MAGIC
MAGIC_plink_dir_merged=${variants_dir}/plink_MAGIC_merged
array_data_dir=${variants_dir}/ARRAY_DATA

### for STITCH
#parameters
infoScoreFilter=0.4
max_missing_param=0.9
minMAF_param=0.025
#chromosome names and lengths file
chr_names_lengths_file=${project_dir}/scripts/Input/chr_names_lengths.txt
#MAGIC lines alignments dir
MAGIC_alignments_dir=${alignments_dir}/Alignments
STITCH_output_dir=${project_dir}/5STITCH/STITCH_output
STITCH_PLINK_dir=${project_dir}/5STITCH/STITCH_PLINK
STITCH_vcfs_dir=${project_dir}/5STITCH/STITCH_vcfs
STITCH_haplotype_dir=${project_dir}/5STITCH/STITCH_haplotypes

### for Association mapping
#pruning parameters
LD_pruning_type=pairwise
window_size=500
window_step=10
LD_threshold=0.99
#directories
FaSTLMM_dir=${project_dir}/6FaSTLMM
phenodir=${FaSTLMM_dir}/phenotypes_2019
FaSTLMM_output_dir=${FaSTLMM_dir}/FaSTLMM_output
FaSTLMM_plots_dir=${FaSTLMM_dir}/FaSTLMM_plots
FaSTLMM_output_topHit_dir=${FaSTLMM_dir}/FaSTLMM_output_topHit
FaSTLMM_plots_topHit_dir=${FaSTLMM_dir}/FaSTLMM_plots_topHit
#plot buffer
buffer_mb=15
buffer=$(( $buffer_mb * 1000000 ))
#LOD drop
tolerance=2

#function to get chromosome subgenome from input number
get_wheat_chr_genome () {
  local SGE_TASK_ID=$1
  #in order to cycle through the chromosomes and subgenomes
  local genomes="A B D"
  local genome_num=$((($SGE_TASK_ID+2)%3 + 1))
  genome=$(echo $genomes | cut -d" " -f$genome_num)
  local chr_tmp=$((($SGE_TASK_ID+2)/3))
  if [ $chr_tmp == 8 ]; then
    chr=Un
    genome=""
  else
    chr=$chr_tmp
  fi
}


#function to get chromosome subgenome and half from input number
get_wheat_chr_genome_half () {
  local SGE_TASK_ID=$1
  #in order to cycle through the chromosomes and subgenomes
  local genomes="A B D"
  local halves="_1 _2"
  local genome_num=$((($SGE_TASK_ID+2)%3 + 1))
  local tmp=$((($SGE_TASK_ID +2)/3))
  local half_num=$((($tmp +1)%2 + 1))
  half=$(echo $halves | cut -d" " -f$half_num)
  genome=$(echo $genomes | cut -d" " -f$genome_num)
  local chr_tmp=$((($SGE_TASK_ID+5)/6))
  if [ $chr_tmp == 8 ]; then
    chr=Un
    genome=""
    half=""
  else
    chr=$chr_tmp
  fi
}

get_founder_samplename () {
  local foundernames="Banco Bersee Brigadier Copain Cordiale Flamingo Gladiator Holdfast Kloka MarisFundin Robigus Slejpner Soissons Spark Steadfast Stetson"
  local founder_samplename=$(echo $foundernames | cut -d" " -f$1)
  echo $founder_samplename
}
