args = commandArgs(trailingOnly=TRUE)
    if (! length(args)==14) {
        stop("supply arguments, 1-K_param 2-G_param 3-chromosome 4-chromosome_name 5-bqFilter 6-minRate_param 7-bamlist_file 8-posfile_dir 9-reference_panel_prefix 10-reference_phred 11-reference_iterations 12-num_iterations 13-tempdir 14-outdir", call.=FALSE)
} else {

library(STITCH)

#parameters
K_param=as.numeric(args[1])
G_param=as.numeric(args[2])
chromosome=args[3]
chromosome_name=args[4]
bqFilter_param=as.numeric(args[5])
minRate_param=as.numeric(args[6])
bamlist_file=args[7]
posfile_dir=args[8]
reference_panel_prefix=args[9]
ref_phred=as.numeric(args[10])
ref_iterations=as.numeric(args[11])
num_iterations=as.numeric(args[12])
tempdir=args[13]
outdir=args[14]

Num_nCores = 1
nCores = 1
system(paste0("mkdir -p ", outdir))
setwd(outdir)
outputdir = paste0(outdir,"/",chromosome_name)

### STEP2: Start imputation using STITCH (NOTE: The chromosome name should stay the same between BAM files and SNP list!) ###

Num_nGen = G_param
Num_K = K_param
Method = "diploid-inbred"

posfile=paste0(posfile_dir,"/",chromosome_name,".SNPs.txt")

ref_hap=paste0(reference_panel_prefix,".hap.gz")
ref_sample=paste0(reference_panel_prefix,".sample")
ref_legend=paste0(reference_panel_prefix,".legend.gz")
#referece populations to use (from sample file)
ref_pop="FOUNDER"

system(paste0("mkdir -p ", outputdir))
system(paste0("mkdir -p ", tempdir))
bamlist = bamlist_file
options(bitmapType='cairo')
STITCH(outputdir = outputdir, 
	chr=chromosome, bqFilter=bqFilter_param, minRate=minRate_param,
	posfile = posfile, 
	reference_sample_file=ref_sample, reference_legend_file=ref_legend, reference_haplotype_file=ref_hap,
	reference_phred=ref_phred, reference_iterations=ref_iterations, reference_populations=ref_pop, 
	bamlist = bamlist, tempdir = tempdir, K = Num_K, outputHaplotypeProbabilities=TRUE, output_haplotype_dosages=TRUE,
	reference_shuffleHaplotypeIterations=NA, 
	shuffleHaplotypeIterations = NA,
	refillIterations = NA,
	generateInputOnly = TRUE, 
	nCores = Num_nCores, nGen = Num_nGen, method = Method)

}


