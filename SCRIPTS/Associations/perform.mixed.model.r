args = commandArgs(trailingOnly=TRUE)
phenotype=args[1]

source("../0Rfunctions/mixed.model.functions.r")
source("../0Rfunctions/functions.r")

phenfile="../phenotypes_2020/founder_and_MAGIC_phenotypes_with_header.tsv"
#phenfile="RData/recombination_pheno.tsv"
phenos<-read_tsv(phenfile)
#phenotype="GR"

grm.datafile="RData/MAGIC.snp.grm.RData"
out_dir<-paste0("../../6Haplotypes/1hap_associations/",phenotype,"/")
plot_dir<-paste0("../../6Haplotypes/2hap_plots/",phenotype,"/")
#load(paste0(out_dir, phenotype, ".hap.output.RData"))
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

nperm=1000
alpha=c(0.05,0.01)

plots_and_tables<-function(hap.output, out_dir, plot_dir, phenotype, iteration=1, alpha=c(0.05,0.01), nperm=100, covar=NULL){
	#convert column names
	haplo.logP<-hap.output$haplo.logP
	haplo.logP$logP[haplo.logP$logP==Inf]<-314
	colnames(haplo.logP)[grepl("beta", colnames(haplo.logP))]<-paste0(founderNames, "_beta")
	colnames(haplo.logP)[grepl("SE", colnames(haplo.logP))]<-paste0(founderNames, "_se")
	#for threshold
	max.logP.perm= hap.output$max.logP.perm
	#output as tables
	write_tsv(tbl_df(haplo.logP), path=paste0(out_dir, phenotype, ".", iteration, ".logP.tsv"))
	write_tsv(tbl_df(data.frame(max.logP.perm=max.logP.perm)), path=paste0(out_dir, phenotype, ".", iteration, ".max.logP.perm.tsv"))
	threshold=sort(max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
	mhtn_plot<-plot_manhattan_haplotypes(haplo.logP, threshold)
	ggsave(mhtn_plot, filename= paste0(plot_dir, phenotype, ".",iteration, ".manhattan.pdf"), height=3, width=7)
	#if there's a QTL that exceeds the threshold, produce plots and table
	if(haplo.logP[which.max(haplo.logP$logP),"logP"]>threshold[1]){
		hit_chr_name<-as.character(filter(haplo.logP, logP==max(logP))$chr[1])
		hap.output.hit<-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/full_haplotypes", chrs=hit_chr_name, thin=1,  n.perm=0, mc.cores=1, covar)
		haplo.logP.hit<-hap.output.hit$haplo.logP
		haplo.logP.hit$logP[haplo.logP.hit$logP==Inf]<-314
		colnames(haplo.logP.hit)[grepl("beta", colnames(haplo.logP.hit))]<-paste0(founderNames, "_beta")
		colnames(haplo.logP.hit)[grepl("SE", colnames(haplo.logP.hit))]<-paste0(founderNames, "_se")
		local<-plot_tophit_haplotypes_raw_pheno(haplo.logP.hit, 1, phenotype, phenfile, dosage.dir="../../5STITCH/STITCH_haplotypes/full_haplotypes", id.column="IID", threshold, logdrop=2, buffer_mb=10, covar=covar)
		outplt<-local$plot
		outplt2<-plot_tophit_haplotypes(haplo.logP.hit, 1, phenotype, threshold, logdrop=2, buffer_mb=10)
		local$info$heritability<-hap.output$vg
		local$info$iteration<-iteration
		write_tsv(local$info, path=paste0(plot_dir, phenotype, ".", iteration, ".info.tsv"))
		ggsave(outplt, filename=paste0(plot_dir, phenotype,".", iteration, ".QTL.pdf"), height=6, width=6)
		ggsave(outplt2, filename=paste0(plot_dir, phenotype,".", iteration, ".QTL.effects.pdf"), height=6, width=6)
	} else {
		local<-plot_tophit_haplotypes_raw_pheno(haplo.logP, 1, phenotype, phenfile, dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", id.column="IID", threshold, logdrop=2, buffer_mb=10)
		local$info[1,]<-NA
		local$info$phenotype=phenotype
		local$info$heritability<-hap.output$vg
		local$info$iteration<-iteration
		write_tsv(local$info, path=paste0(plot_dir, phenotype, ".", iteration, ".info.tsv"))
	}
}

#first iteration
iteration=1
hap.output <-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs=paste0(rep(seq(1:7), each=3), c("A","B","D")), thin=1,  n.perm=nperm, mc.cores=1)
save(hap.output, file=paste0(out_dir, phenotype, ".hap.output.RData"))
plots_and_tables(hap.output, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm)

#if there's a QTL that exceeds the threshold, fit it as a covariate
threshold1=sort(hap.output$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
if(hap.output$haplo.logP[which.max(hap.output$haplo.logP$logP),"logP"]> threshold1[1]){
	hit_locus_name<-as.character(hap.output$haplo.logP[which.max(hap.output$haplo.logP$logP),"locus"])
	hit_chr_name<-as.character(hap.output$haplo.logP[which.max(hap.output$haplo.logP$logP),"chr"])
	covar1.dose<-load_named_dose(dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", hit_chr_name, hit_locus_name)	
	hap.output.covar1<-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs=paste0(rep(seq(1:7), each=3), c("A","B","D")), thin=1,  n.perm=nperm, mc.cores=1, covar=covar1.dose)
	iteration=2
	plots_and_tables(hap.output.covar1, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar=covar1.dose)

	#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
	threshold2=sort(hap.output.covar1$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
	if(hap.output.covar1$haplo.logP[which.max(hap.output.covar1$haplo.logP$logP),"logP"]> threshold2[1]){
		iteration=3
		hit_locus_name<-as.character(hap.output.covar1$haplo.logP[which.max(hap.output.covar1$haplo.logP$logP),"locus"])
		hit_chr_name<-as.character(hap.output.covar1$haplo.logP[which.max(hap.output.covar1$haplo.logP$logP),"chr"])
		covar2.dose<-load_named_dose(dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", hit_chr_name, hit_locus_name)	
		hap.output.covar2<-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs=paste0(rep(seq(1:7), each=3), c("A","B","D")), thin=1,  n.perm=nperm, mc.cores=1, covar=cbind(covar1.dose,covar2.dose))
		plots_and_tables(hap.output.covar2, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar=cbind(covar1.dose,covar2.dose))
		
		#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
		threshold3=sort(hap.output.covar2$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
		if(hap.output.covar2$haplo.logP[which.max(hap.output.covar2$haplo.logP$logP),"logP"]> threshold3[1]){
			iteration=4
			hit_locus_name<-as.character(hap.output.covar2$haplo.logP[which.max(hap.output.covar2$haplo.logP$logP),"locus"])
			hit_chr_name<-as.character(hap.output.covar2$haplo.logP[which.max(hap.output.covar2$haplo.logP$logP),"chr"])
			covar3.dose<-load_named_dose(dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", hit_chr_name, hit_locus_name)	
			hap.output.covar3<-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs=paste0(rep(seq(1:7), each=3), c("A","B","D")), thin=1,  n.perm=nperm, mc.cores=1, covar=cbind(covar1.dose,covar2.dose, covar3.dose))
			plots_and_tables(hap.output.covar3, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar=cbind(covar1.dose,covar2.dose, covar3.dose))
			
				#now add another set of dosages as covariates if there's another QTL that exceeds the threshold
			threshold4=sort(hap.output.covar3$max.logP.perm, decreasing=TRUE)[ceiling(alpha*nperm)]
			if(hap.output.covar3$haplo.logP[which.max(hap.output.covar3$haplo.logP$logP),"logP"]> threshold4[1]){
				iteration=5
				hit_locus_name<-as.character(hap.output.covar3$haplo.logP[which.max(hap.output.covar3$haplo.logP$logP),"locus"])
				hit_chr_name<-as.character(hap.output.covar3$haplo.logP[which.max(hap.output.covar3$haplo.logP$logP),"chr"])
				covar4.dose<-load_named_dose(dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", hit_chr_name, hit_locus_name)	
				hap.output.covar4<-haplotype.dosage.mixed.model(phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs=paste0(rep(seq(1:7), each=3), c("A","B","D")), thin=1,  n.perm=nperm, mc.cores=1, covar=cbind(covar1.dose,covar2.dose, covar3.dose, covar4.dose))
				plots_and_tables(hap.output.covar4, out_dir, plot_dir, phenotype, iteration=iteration, alpha=alpha, nperm=nperm, covar=cbind(covar1.dose,covar2.dose, covar3.dose, covar4.dose))
			}
		}
	}
}

### maximum of four covariates fit by this code

