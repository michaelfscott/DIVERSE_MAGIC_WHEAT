args = commandArgs(trailingOnly=TRUE)
phenotype=args[1]

source("../0Rfunctions/mixed.model.functions.r")
source("../0Rfunctions/functions.r")

#load genotypes
traw_path_MAGIC_full="../../5STITCH/STITCH_PLINK/ALLchr.SNPs.traw"
traw_path_MAGIC_pruned="../../5STITCH/STITCH_PLINK/ALLchr.SNPs.pruned.traw"
traw_path_founders="../../4Variants/plink_Founders_merged/ALLchr.traw"
array_data_prefix="../../4Variants/ARRAY_DATA/"
lines_name_conversion_path=paste0(array_data_prefix,"line_name_conversion.csv")
traw_MAGIC_pruned<-read_tsv(traw_path_MAGIC_pruned) %>% rearrange_and_rename_samples() 
traw_MAGIC_full<-read_tsv(traw_path_MAGIC_full) %>% rearrange_and_rename_samples()
traw_founders<-read_tsv(traw_path_founders) %>% rearrange_and_rename_samples() %>% filter(!chromosome=="Un")
###

phenfile="../phenotypes_2020/founder_and_MAGIC_phenotypes_with_header.tsv"
#phenfile="../phenotypes_2020/founder_and_MAGIC_compound_phenotypes_with_header.tsv"
#phenfile="RData/recombination_pheno.tsv"
phenos<-read_tsv(phenfile)
#phenotype="GR"
#phenotype="PY_1"
#phenotype="PY_2"
#phenotype="PC11"
phenotypes=colnames(phenos)[-c(1:2)]
nperm=1000
alpha=c(0.05,0.01)

grm.datafile="RData/MAGIC.snp.grm.RData"
out_dir<-paste0("../../6SNPs/1SNP_associations/",phenotype,"/")
plot_dir<-paste0("../../6SNPs/2SNP_plots/",phenotype,"/")
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

#### SNP based

plots_and_tables<-function(snp.output, out_dir, plot_dir, phenotype, iteration=1, alpha=c(0.05,0.01), nperm=100, covar=NULL){
	#convert underflow p values to 314
	logP<-snp.output$logP
	logP$logP[logP$logP==Inf]<-314
	max.logP.perm= snp.output$max.logP.perm
	#output as tables
	write_tsv(tbl_df(logP), path=paste0(out_dir, phenotype, ".",iteration, ".logP.tsv"))
	write_tsv(tbl_df(data.frame(max.logP.perm=max.logP.perm)), path=paste0(out_dir, phenotype, ".",iteration, ".max.logP.perm.tsv"))
	threshold=sort(max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
	mhtn_plot<-plot_manhattan_snps(logP, threshold)
	ggsave(mhtn_plot, filename= paste0(plot_dir, phenotype, ".",iteration, ".manhattan.pdf"), height=3, width=7)
	#if there's a QTL that exceeds the threshold, produce plots and table
	if(logP[which.max(logP$logP),"logP"]>threshold[1]){
		hit_chr_name<-filter(logP, logP==max(logP, na.rm=T))$chromosome[1]
		snp.output.hit<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", filter(traw_MAGIC_full, chromosome==hit_chr_name), n.perm= 0, covar)
		snp.output.hit$logP[snp.output.hit$logP==Inf]<-314
		local<-plot_tophit_snps(snp.output.hit$logP, 1, phenotype, phenfile, id.column="IID", threshold, filter(traw_MAGIC_full, chromosome==hit_chr_name), filter(traw_founders, chromosome== hit_chr_name), logdrop=2, buffer_mb=10, max_QTL_extension=5000000, covar=covar)
		outplt<-local$plot
		local$info$heritability<-snp.output$vg
		local$info$iteration<-iteration
		write_tsv(local$info, path=paste0(plot_dir, phenotype, ".", iteration, ".info.tsv"))
		ggsave(outplt, filename=paste0(plot_dir, phenotype,".", iteration, ".QTL.pdf"), height=6, width=6)
	} else {
		hit_chr_name<-filter(logP, logP==max(logP))$chromosome[1]
		local<-plot_tophit_snps(snp.output$logP, 1, phenotype, phenfile, id.column="IID", threshold, filter(traw_MAGIC_full, chromosome==hit_chr_name), filter(traw_founders, chromosome== hit_chr_name), logdrop=2, buffer_mb=10)
		local$info[1,]<-NA
		local$info$phenotype=phenotype
		local$info$heritability<-snp.output$vg
		local$info$iteration<-iteration
		write_tsv(local$info, path=paste0(plot_dir, phenotype, ".", iteration, ".info.tsv"))
	}
}

#first iteration
iteration=1
snp.output<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw_MAGIC_pruned, n.perm= nperm)
save(snp.output, file=paste0(out_dir, phenotype, ".snp.output.RData"))
plots_and_tables(snp.output, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm)

#if there's a QTL that exceeds the threshold, fit it as a covariate
threshold1=sort(snp.output$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
if(snp.output$logP[which.max(snp.output$logP$logP),"logP"]> threshold1[1]){
	hit_locus_name<-as.character(snp.output$logP[which.max(snp.output$logP$logP),"SNP"])
	covar_tmp<-traw_MAGIC_full[traw_MAGIC_full$SNP==hit_locus_name,-(1:6)]
	covar<-t(covar_tmp)
	snp.output.covar1<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw_MAGIC_pruned, n.perm= nperm, covar=covar)
	iteration=2
	plots_and_tables(snp.output.covar1, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar)
	#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
	threshold2=sort(snp.output.covar1$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
	if(snp.output.covar1$logP[which.max(snp.output.covar1$logP$logP),"logP"]> threshold2[1]){
		hit_locus_name<-as.character(snp.output.covar1$logP[which.max(snp.output.covar1$logP$logP),"SNP"])
		covar_tmp<-t(traw_MAGIC_full[traw_MAGIC_full$SNP==hit_locus_name,-(1:6)])
		covar<-cbind(covar,covar_tmp)
		snp.output.covar2<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw_MAGIC_pruned, n.perm= nperm, covar=covar)
		iteration=3
		plots_and_tables(snp.output.covar2, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar)	
		#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
		threshold3=sort(snp.output.covar2$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
		if(snp.output.covar2$logP[which.max(snp.output.covar2$logP$logP),"logP"]> threshold3[1]){
			hit_locus_name<-as.character(snp.output.covar2$logP[which.max(snp.output.covar2$logP$logP),"SNP"])
			covar_tmp<-t(traw_MAGIC_full[traw_MAGIC_full$SNP==hit_locus_name,-(1:6)])
			covar<-cbind(covar,covar_tmp)
			snp.output.covar3<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw_MAGIC_pruned, n.perm= nperm, covar=covar)
			iteration=4
			plots_and_tables(snp.output.covar3, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar)
			#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
			threshold4=sort(snp.output.covar3$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
			if(snp.output.covar3$logP[which.max(snp.output.covar3$logP$logP),"logP"]> threshold4[1]){
				hit_locus_name<-as.character(snp.output.covar3$logP[which.max(snp.output.covar3$logP$logP),"SNP"])
				covar_tmp<-t(traw_MAGIC_full[traw_MAGIC_full$SNP==hit_locus_name,-(1:6)])
				covar<-cbind(covar,covar_tmp)
				snp.output.covar4<-mixed.model.snps(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw_MAGIC_pruned, n.perm= nperm, covar=covar)
				iteration=5
				plots_and_tables(snp.output.covar4, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar)
			}
		}
	}
}

### maximum of four covariates fit by this code