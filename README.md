# DIVERSE_MAGIC_WHEAT

Scripts associated with the analysis of the NIAB Diverse MAGIC wheat population. More information about this population can be found at http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/index.html.

# vcf2df.pl

Perl script for extracting haplotype dosages from vcfs outputted by STITCH (Davies et al., 2016). 

# mixed.model.functions.r

Functions for performing mixed model GWAS at the level of haplotypes or SNPs. 

Key functions are haplotype.dosage.mixed.model and mixed.model.snps. These require a file containing phenotype measurements (phenfile), the column name of the phenotype (phenotype) and the column name containing the sample IDs (id.column). A genetic relationship matrix (object "grm") is also loaded as an RData object using the grm.datafile option. n.perm specifies the number of phenotype permutations to perform to establish genomewide significance thresholds. For SNP-based mapping, SNP genotypes and locations can be read into R from a PLINK .traw file. For haplotype-based mapping, the function import.stitch.haplotype.dosages.wrapper should be performed on the output from vcf2df and the location specified in the haplotype.dosage.mixed.model function using the argument dosage.dir. The additional arguments chrs and thin specify the chromosome prefix for each chromosome dosage file and the level to which haplotypes should be thinned to reduce comuptational time (thin=1 to include all sites). 
