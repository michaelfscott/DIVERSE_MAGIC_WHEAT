# DIVERSE_MAGIC_WHEAT

Scripts associated with the analysis of the NIAB Diverse MAGIC wheat population. More information about this population can be found at http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/index.html.

# SCRIPTS

The SCRIPTS directory contains shell scripts used for alignment, variant calling, imputation, and association mapping. These scripts were tested on the UCL Research Computing HPC cluster, which uses a Sun Grid Engine scheduler. The flags for submission via SGE are at the head of each submission script. The required software is specified/loaded in "Input/software.sh". This will have to be modified to reflect the local install locations. The locations for accessing raw data and writing intermediate files is specified in Input/data_and_variables.sh.

# vcf2df.pl

Perl script for extracting haplotype dosages from vcfs outputted by STITCH (Davies et al., 2016). 

# mixed.model.functions.r

Functions for performing mixed model GWAS at the level of haplotypes or SNPs. 

Key functions are haplotype.dosage.mixed.model and mixed.model.snps. These require a file containing phenotype measurements (phenfile), the column name of the phenotype (phenotype) and the column name containing the sample IDs (id.column). A genetic relationship matrix (object "grm") is also loaded as an RData object using the grm.datafile option. n.perm specifies the number of phenotype permutations to perform to establish genomewide significance thresholds. For SNP-based mapping, SNP genotypes and locations can be read into R from a PLINK .traw file. For haplotype-based mapping, the function import.stitch.haplotype.dosages.wrapper should be performed on the output from vcf2df and the location specified in the haplotype.dosage.mixed.model function using the argument dosage.dir. The additional arguments chrs and thin specify the chromosome prefix for each chromosome dosage file and the level to which haplotypes should be thinned to reduce comuptational time (thin=1 to include all sites). 

# compare.founders.r

Contains founder.mosaics function that compares founder genotypes to determine pairwise similarity by using a dynamic programming algorithm. The function has three arguments (1) a traw file genotypes, outputted by PLINK (2) a penalty score for mismatching genotypes (3) a jumpcost penalty applied when inferring a transition between matching and non-matching states. 
