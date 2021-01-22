library(tidyr)
library(readr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(grid)
library(ggpubr)
library(viridis)
library(RColorBrewer)

#for conversion of alleles when blast alignments are to reverse strand
complement_map<-setNames(c("A","C","T","G"), c("T", "G", "A", "C"))
founderNames=c("Banco", "Bersee", "Brigadier", "Copain", "Cordiale", "Flamingo", "Gladiator", "Holdfast", "Kloka", "MarisFundin", "Robigus", "Slejpner", "Soissons", "Spark", "Steadfast", "Stetson")

#to add to plots
Alabel<-textGrob("a", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Blabel<-textGrob("b", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Clabel<-textGrob("c", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Dlabel<-textGrob("d", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Elabel<-textGrob("e", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Flabel<-textGrob("f", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))

#colour scheme
#cp<-c("SteelBlue", "MediumVioletRed", "#B79F00", "#17377A", "#FF7844", "grey30", "#E5C616")
cp<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02", "#A6761D", "#666666")
founder_call_number<-1131251

#phenotypes
phenology=c("GR","PHS","SH","FLA_1", "FLA_2","FLF","GLA_16_12_15","GLA_17_01_06", "GLA_17_03_17", "GLA_17_03_31","GLA_17_04_20","GLA_17_04_28", "GLA_17_11_01", "GLA_17_11_21", "GLA_17_12_13", "GLA_18_02_06","GLA_18_03_13","GLA_18_03_27","GLA_18_04_10","GLA_18_04_19","GLA_18_04_27","GLA_18_06_25","GS39_1","GS55_1","GS65_1","GS39_2","GS55_2","GS65_2","JGH_1","JGH_2","JGH_N","GY_1","GY_2")
size<-c("LOD","HET_1","HET_2","HFLB_1","HFLB_2","FLED_1", "FLED_2","EL_1", "EL_2","ETA","FLL_1", "FLL_2","FLW_1","FLW_2","GY_1", "GY_2")
yld<-c("BIS_1","TS","GPS","GPE","EW","GA_1","GA_2","GL_1","GL_2","GW_1","GW_2","TGW_1","TGW_2","SW","GPC_1", "GPC_2","GY_1","GY_2")
other_phens<-c("LGLAU", "GLAU_1","GLAU_2", "AWN", "PIG_1", "PIG_2", "SPIG", "YR_1",  "YR_2", "YR_N", "GY_1", "GY_2")

text_size=8
plot_theme1<-theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill=NA),
	legend.title=element_text(size=text_size),
    legend.text=element_text(size=text_size),
    axis.title= element_text(size=text_size), axis.text=element_text(size=text_size),
    strip.text = element_text(size=text_size)) 
plot_theme2<-plot_theme1+
  theme(panel.grid.minor.x = element_blank(), 
    panel.grid.major.x=element_blank(), 
    axis.text.x=element_text(angle=45, hjust=1))

 #add two halves of each chromosome together
get_chr_startpos<-function(traw){
  start_positions_per_chr<-traw %>% group_by(CHR) %>%
    summarize(
    startpos= as.numeric(word(word(SNP[1],2,sep=":"),1,sep="-"))-1,
    endpos= as.numeric(word(word(SNP[1],2,sep=":"),2,sep="-")),
    chromosome=word(SNP[1],1,sep=":"))
  return(start_positions_per_chr)
}

join_chromosome_halves_rename_samples<-function(traw){
  start_positions_per_chr<-get_chr_startpos(traw) 
  samplenames<-colnames(traw)[7:ncol(traw)]
  edited_samplenames<-word(samplenames,1,sep="_")
  traw_new<-left_join(traw,start_positions_per_chr, by="CHR") %>% ungroup() %>%
    mutate(position=POS+startpos)
  #re-order columns
  column_order<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", sort(samplenames))
  traw_new <-traw_new[,column_order]
  colnames(traw_new)<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", sort(edited_samplenames))
  return(traw_new)
}

rearrange_and_rename_samples<-function(traw){
  start_positions_per_chr<-traw %>% group_by(CHR) %>%
    summarize(
    startpos= 0,
    chromosome=word(SNP[1],1,sep=":"))
  samplenames<-colnames(traw)[7:ncol(traw)]
  edited_samplenames<-word(samplenames,1,sep="_")
  traw_new<-left_join(traw,start_positions_per_chr, by="CHR") %>% ungroup() %>%
    mutate(position=POS+startpos)
  #re-order columns
  column_order<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", sort(samplenames))
  traw_new <-traw_new[,column_order]
  colnames(traw_new)<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", sort(edited_samplenames))
  return(traw_new)
}

convert_line_names<-function(traw, lines_name_conversion_path){
	convert_names<-read_csv(lines_name_conversion_path) 
	conversion_map=setNames(convert_names$line_name, convert_names$line_code_A)
	replacement<-colnames(traw)[-(1:6)]
	replacement[]<-conversion_map[replacement]
	colnames(traw)<-c(colnames(traw)[1:6], replacement)
	return(traw)
}

array_calls_with_position_and_allele<-function(genotyping_calls_path, blast_path, axiom_info_path, founders=FALSE){
	if(founders){
		genos<-read.table(file=genotyping_calls_path, stringsAsFactors=FALSE)
		genos <-tbl_df(tibble::rownames_to_column(genos, var="probeset_id"))	
	} else {
		#get the genotypes and remove the extraneous parts of the samplenames
		genos<-read.table(file=genotyping_calls_path, stringsAsFactors=FALSE)
		colnames(genos)<-str_extract(colnames(genos), "[A-Z]\\.[0-9][0-9][0-9]\\.[0-9]+")
		genos <-tbl_df(tibble::rownames_to_column(genos, var="probeset_id"))
	}
	#read in the blast results
	blast_data<-read_tsv(blast_path, comment = "#", col_names=c("probeset_id", "subject_acc.ver", "perc_identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"))
	#extract the allele information
	axiom_array_info<-read_tsv(axiom_info_path) %>%
		mutate(alleles=str_extract(Sequence, "[A-Z]/[A-Z]"), 
			allele1=str_sub(alleles,1,1), 
			allele2=str_sub(alleles,3,3), 
			probeset_id=X35K.SNPId) %>%
		select(probeset_id, allele1, allele2)
	#combine data
	left_join(left_join(axiom_array_info, blast_data, by="probeset_id"), genos, by="probeset_id") %>%
		filter(gap_opens==0) %>%
		mutate(
			chromosome=str_replace(subject_acc.ver, "chr", ""),
			position= ifelse((s.end-s.start)>0, s.start+q.start+34, s.start-q.start-34),
			allele1= ifelse((s.end-s.start)>0, allele1, complement_map[allele1]),
			allele2= ifelse((s.end-s.start)>0, allele2, complement_map[allele2]))
}

compare_condordances_traw_array<-function(concordances){
	array_lines<-colnames(concordances)
	seq_lines<-rownames(concordances)
	sequence_line=str_replace(rownames(concordances), "\\.y", "")
	best_match=array_lines[apply(concordances, 1, which.max)]
	best_match_value=apply(concordances, 1, max)
	same_line_value<-sapply(1:length(seq_lines), function(x)(concordances[seq_lines[x],sequence_line[x]]))
	tbl_df(data.frame(sequence_line=sequence_line, best_match_array=best_match, best_match_value=best_match_value, same_line_value=same_line_value, stringsAsFactors=FALSE))
}

get_concordance_same_line_name_by_SNP <-function(overlap){
	line_names<-gsub("\\.x","",grep("\\.x",colnames(overlap), value=TRUE))
	tmp<-lapply(line_names,function(x)(overlap[,c(paste0(x,".x"), paste0(x,".y"))]))
	tmp2<-lapply(tmp,function(x)(x[,1]==x[,2]))
	apply(do.call(cbind,tmp2),1,sum, na.rm=TRUE)/apply(!is.na(do.call(cbind,tmp2)),1,sum)
}

get_concordance_same_line_name <-function(overlap){
	line_names<-gsub("\\.x","",grep("\\.x",colnames(overlap), value=TRUE))
	tmp<-lapply(line_names,function(x)(overlap[,c(paste0(x,".x"), paste0(x,".y"))]))
	tmp2<-lapply(tmp,function(x)(x[,1]==x[,2]))
	apply(do.call(cbind,tmp2),2,sum, na.rm=TRUE)/apply(!is.na(do.call(cbind,tmp2)),2,sum)
}

get_Rsq_same_line_name <-function(overlap){
	line_names<-gsub("\\.x","",grep("\\.x",colnames(overlap), value=TRUE))
	tmp<-lapply(line_names,function(x)(overlap[,c(paste0(x,".x"), paste0(x,".y"))]))
	tmp2<-lapply(tmp,function(x)((cor(x[,1],x[,2],use="complete.obs")^2)))
	apply(do.call(cbind,tmp2),2,sum, na.rm=TRUE)/apply(!is.na(do.call(cbind,tmp2)),2,sum)
}

compare_traw_HC_array<-function(traw, array){
	overlap<-inner_join(array, traw, by=c("chromosome", "position"))
	concordance<-get_concordance_same_line_name(overlap)
	Rsq<-get_Rsq_same_line_name(overlap)
	calls<-select(traw, matches("[0-9][0-9][0-9].[0-9]")) %>%
		summarise_all(funs(sum(!is.na(.)))) 
	output<-tbl_df(data.frame(line_name=colnames(calls), call_rate=as.numeric(calls)/founder_call_number, stringsAsFactors=F)) %>%
		left_join(tbl_df(data.frame(line_name=str_replace(names(Rsq), ".y",""), Rsq=as.numeric(Rsq), stringsAsFactors=F)) , by="line_name") %>%
		left_join(tbl_df(data.frame(line_name=str_replace(names(concordance), ".x",""), concordance=as.numeric(concordance), stringsAsFactors=F)) , by="line_name" ) 
	return(output)
}

update_pos<-function(pos){
	summarise(pos, 
		CHR=CHR[1], 
		chrom=word(CHR[1],1,sep=":"), 
		start_pos=as.numeric(word(word(CHR[1],2,sep=":"), 1,sep="-"))) %>%
	right_join(pos) %>% 
	mutate(CHR=chrom, POS=POS+start_pos-1) %>% select(-c("chrom", "start_pos"))
}

plot_manhattan_snps<-function(logP, threshold){
	logP %>%
	ggplot(aes(x=position, y=logP, colour = as.factor((as.numeric(as.factor(chromosome))-1) %% 3))) + 
		geom_point(size=0.4) + 
		geom_hline(yintercept=threshold, linetype="dashed") +
		scale_color_manual(values=cp[1:3]) +
		scale_x_continuous(labels = NULL, breaks=0) +
		facet_grid(.~chromosome, scale = "free_x", space="free_x", switch="x") +
		xlab("Chromosome") +
		ylab(expression(paste("-log"[10], plain(P)))) +
		theme(panel.spacing = unit(0, "lines")) +
		theme(panel.background = element_blank()) +
		theme(panel.grid.major = element_line(colour="grey85"), panel.grid.minor = element_line(colour="grey85")) +
		theme(legend.position="none") +
		theme(axis.line.y = element_line(colour = "black")) +
		theme(axis.ticks.x = element_blank()) +
		theme(strip.background = element_blank()) +
		NULL
}

plot_manhattan_haplotypes<-function(haplo.logP, threshold){
	tidyr::gather(tbl_df(haplo.logP), key="part", value="position", c("start", "end")) %>%
		ggplot(aes(x=position, y=logP, colour = as.factor((as.numeric(as.factor(chr))-1) %% 3))) + 
			geom_line(size=0.4) + 
			geom_hline(yintercept=threshold, linetype="dashed") +
			scale_color_manual(values=cp[1:3]) +
			scale_x_continuous(labels = NULL, breaks=0) +
			facet_grid(.~chr, scale = "free_x", space="free_x", switch="x") +
			xlab("Chromosome") +
			ylab(expression(paste("-log"[10], plain(P)))) +
			theme(panel.spacing = unit(0, "lines")) +
			theme(panel.background = element_blank()) +
			theme(panel.grid.major = element_line(colour="grey85"), panel.grid.minor = element_line(colour="grey85")) +
			theme(legend.position="none") +
			theme(axis.line.y = element_line(colour = "black")) +
			theme(axis.ticks.x = element_blank()) +
			theme(strip.background = element_blank()) +
			NULL
}

break_founder_names<-function(tmp0){
	tmp0_v<-str_split(tmp0, pattern=", ")[[1]]
	if(length(tmp0_v)<=4){tmp0<-paste(tmp0_v, collapse=", ")} else {
		if(length(tmp0_v)<=8){tmp0<-paste0(paste(tmp0_v[1:4][!is.na(tmp0_v[1:4])], collapse=", "),"\n", paste(tmp0_v[5:length(tmp0_v)], collapse=", "))} else {
			if(length(tmp0_v)<=12){tmp0<-paste0(paste(tmp0_v[1:4], collapse=", "),"\n", paste(tmp0_v[5:8], collapse=", "),"\n", paste(tmp0_v[9:length(tmp0_v)], collapse=", "))} else {
				tmp0<-paste0(paste(tmp0_v[1:4], collapse=", "),"\n", paste(tmp0_v[5:8], collapse=", "),"\n", paste(tmp0_v[9:12], collapse=", "),"\n", paste(tmp0_v[13:length(tmp0_v)], collapse=", "))}}}
	return(tmp0)
}

convert_SDP<-function(SDP){
	gsub("tmp","0",gsub("0","2",gsub("2","tmp", SDP)))
}

plot_tophit_snps<-function(logP, chr_hitnum=1, phenotype, phenfile, id.column="IID", threshold, traw_MAGIC, traw_founders, logdrop=2, buffer_mb=10, max_QTL_extension=5000000, covar=NULL){
	#get chromosome hitnum order
	tmp<-tbl_df(logP) %>% group_by(chromosome) %>% summarise(mx=max(logP)) %>% arrange(desc(mx)) %>% data.frame()
	chr_tophits<-as.character(tmp[,1])
	tophit_tbl<-tbl_df(logP) %>%
		filter(chromosome==chr_tophits[chr_hitnum]) %>%
		filter(logP==max(logP))
	lgdrp=max(c(logdrop, 0.1*max(tophit_tbl$logP)))
	#get SDP for SNPs
	logP_SDP<-filter(tbl_df(traw_founders), chromosome==chr_tophits[chr_hitnum]) %>%
		right_join(filter(tbl_df(logP), chromosome==chr_tophits[chr_hitnum]), by = c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT")) %>%
		mutate(SDP=paste0(Banco, Bersee, Brigadier, Copain, Cordiale, Flamingo, Gladiator, Holdfast, Kloka, MarisFundin, Robigus, Slejpner, Soissons, Spark, Steadfast, Stetson)) %>% select(-founderNames)
	#core range of highly associated SNPs
	core_tophit_range<-logP_SDP %>%
		filter(chromosome== as.character(tophit_tbl$chromosome[1]), logP>mean(tophit_tbl$logP)-lgdrp) %>%
		group_by(SDP) %>%
		summarise(start=min(position), end=max(position))
	#go to next SNP with the same SDP as those in the core range
	min_ext_range<-filter(logP_SDP, position<min(core_tophit_range$start), SDP %in% c(core_tophit_range$SDP, convert_SDP(core_tophit_range$SDP))) %>%
		filter(position==max(position)) %>%
		mutate(position=max(c(position, min(core_tophit_range$start)-max_QTL_extension)))
	if(nrow(min_ext_range)==0){min_ext_range<-filter(logP_SDP, position==min(core_tophit_range$start)) %>% mutate(position=max(0,position-max_QTL_extension))}
	max_ext_range<-filter(logP_SDP, position>max(core_tophit_range$end), SDP %in% c(core_tophit_range$SDP, convert_SDP(core_tophit_range$SDP))) %>%
		filter(position==min(position)) %>%
		mutate(position=min(c(position,max(core_tophit_range$end)+max_QTL_extension)))
	if(nrow(max_ext_range)==0){max_ext_range<-filter(logP_SDP, position==max(core_tophit_range$end)) %>% mutate(position=max(0,position+max_QTL_extension))}
	tophit_range<-data.frame(start=min_ext_range$position[1], end=max_ext_range$position[1])
	region <- data.frame(xmin=rep(tophit_range$start[1]/1000000,2), xmax=rep(tophit_range$end[1]/1000000,2), ymin=rep(-Inf,2), ymax=rep(Inf,2))
	logP_plt_data<-tbl_df(logP) %>%
		filter(chromosome==as.character(tophit_tbl$chromosome[1]), position>tophit_range$start-(buffer_mb*1000000), position<tophit_range$end+(buffer_mb*1000000)) 
	region_plot<-ggplot() + 
		    geom_rect(data= region, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.3) +
			geom_point(data=logP_plt_data, aes(x=position/1000000, y=logP), colour=cp[1]) + 
			theme_minimal() +
			geom_hline(yintercept=threshold, linetype="dashed") +
			ylab(expression(paste("-log"[10], plain(P)))) +
			xlab("position (mb)") +
			theme(panel.border = element_rect(colour = "black", fill=NA),
	    	panel.grid.minor=element_blank()) +
	    NULL
	# load phenotypes
	phenos = read.delim(phenfile)
	y = phenos[[phenotype]]
	names(y) = phenos[[id.column]]
	y = y[!is.na(y)] 
	#load topSNP genotypes
	geno_tmp<-filter(traw_MAGIC, SNP %in% tophit_tbl$SNP[1]) %>% data.frame()
	geno_tmp_fnd<-filter(traw_founders, SNP %in% tophit_tbl$SNP[1]) %>% data.frame()
	geno_tmp<-geno_tmp[,-c(1:6)]
	geno_tmp_fnd<-geno_tmp_fnd[,-c(1:6)]
	geno<-c(as.numeric(geno_tmp[1,]), as.numeric(geno_tmp_fnd[1,]))
	names(geno)<-c(colnames(geno_tmp), colnames(geno_tmp_fnd))
	# get phenotype in the same order as dosages
	use.ids = sort(intersect(names(y), names(geno)))
	#use.ids_fnd<-c(use.ids,founderNames)
	y_fnd = y[match(use.ids,names(y),nomatch=0)]
	#prepare for plotting
	phenos_plotter<-tbl_df(data.frame(name=names(y_fnd), y=y_fnd, geno=geno[names(y_fnd)])) 
	fnd0<-filter(phenos_plotter, geno==0, name %in% founderNames) %>% summarise(carriers=paste(name, collapse=", ")) %>% as.character()
	tmp0<-break_founder_names(fnd0)
	fnd2<-filter(phenos_plotter, geno==2, name %in% founderNames) %>% summarise(carriers=paste(name, collapse=", ")) %>% as.character()
	tmp2<-break_founder_names(fnd2)
	phenos_plotter<-mutate(phenos_plotter, type=ifelse(!name %in% founderNames, "MAGIC line", ifelse(geno==2, tmp2, tmp0)) )
	snp_effect_plot<-ggplot() +
		geom_jitter(data=phenos_plotter, aes(x=geno, y=y, colour=factor(type, levels=c("MAGIC line", tmp0, tmp2)), alpha= factor(type, levels=c("MAGIC line", tmp0, tmp2)), shape=factor(type, levels=c("MAGIC line", tmp0, tmp2))), width=0.1, height=0) +
		geom_smooth(data=filter(phenos_plotter, type=="MAGIC line"), aes(x=geno, y=y), method='lm', colour="grey30") +
		scale_alpha_manual(values=c(0.3,1,1)) +
		scale_colour_manual(values=cp[1:3]) +
		scale_shape_manual(values=c(16,19,19)) +
		scale_x_continuous(breaks=seq(0,2,by=1)) +
		xlab(paste0("SNP dosage (",tophit_tbl$SNP[1], ")")) +
    	ylab(paste0(str_replace_all(phenotype, "_"," "))) +
	    theme_minimal() + 
    	theme(
    		legend.title=element_blank(),
    		panel.border = element_rect(colour = "black", fill=NA),
	    	panel.grid.minor=element_blank()) +
	    NULL
	plt<-ggarrange(region_plot, snp_effect_plot, ncol=1, nrow=2, heights=c(5/10, 5/10), labels="auto")
	outplot<-annotate_figure(plt, top = text_grob(paste0(str_replace_all(phenotype, "_"," "), " chromosome ", as.character(tophit_tbl$chromosome[1])), face = "bold"))
	#get rsq
	# get phenotype in the same order as dosages
	use.ids_MAGIC = sort(intersect(names(y), names(geno_tmp)))
	y_MAGIC = y[match(use.ids_MAGIC,names(y),nomatch=0)]
	geno_MAGIC<-as.numeric(geno_tmp[match(use.ids_MAGIC,names(geno_tmp),nomatch=0)])
	f2<-lm(y_MAGIC~geno_MAGIC)
	snp_effect=f2$coef[2]
	#get full model fit if there are covariates
	if(is.null(covar)){
		rsq=data.frame(rsq=summary(f2)$r.squared, rsq.adj=summary(f2)$adj.r.squared, rsq_full=summary(f2)$r.squared, rsq_full.adj=summary(f2)$adj.r.squared)
	} else {
		covar_MAGIC<-covar[match(use.ids_MAGIC,rownames(covar),nomatch=0),]
		f3<-lm(y_MAGIC~ geno_MAGIC+covar_MAGIC)
		rsq=data.frame(rsq=summary(f2)$r.squared, rsq.adj=summary(f2)$adj.r.squared, rsq_full=summary(f3)$r.squared, rsq_full.adj=summary(f3)$adj.r.squared)
	}
	output_info<-data.frame(phenotype=phenotype, chromosome=as.character(tophit_tbl$chromosome[1]), min_position=tophit_range$start, max_position=tophit_range$end, core_start_position=min(core_tophit_range$start), core_end_position=max(core_tophit_range$end), peak_start= min(tophit_tbl$position), peak_end=max(tophit_tbl$position), peak_logP=max(tophit_tbl$logP), threshold05=threshold[1], threshold01=threshold[2], rsq, snp_effect, peak_carriers0=fnd0, peak_carriers2=fnd2, all_highly_associated_SDP= paste(names(table(core_tophit_range$SDP)), collapse=", "))
	return(list(plot=outplot, info=output_info))
}


plot_tophit_haplotypes<-function(haplo.logP, chr_hitnum, phenotype, threshold, logdrop=2, buffer_mb=10, covar=NULL){
	#get chromosome hitnum order
	tmp<-tbl_df(haplo.logP) %>% group_by(chr) %>% summarise(mx=max(logP)) %>% arrange(desc(mx)) %>% data.frame()
	chr_tophits<-as.character(tmp[,1])
	tophit_tbl<-tbl_df(haplo.logP) %>%
		filter(chr==chr_tophits[chr_hitnum]) %>%
		filter(logP==max(logP))
	lgdrp=max(c(logdrop, 0.1*max(tophit_tbl$logP)))
	tophit_range<-tbl_df(haplo.logP) %>%
		filter(chr== as.character(tophit_tbl$chr[1]), logP>mean(tophit_tbl$logP)-lgdrp) %>%
		summarise(start=min(start), end=max(end))
	region <- data.frame(xmin=rep(tophit_range$start[1]/1000000,2), xmax=rep(tophit_range$end[1]/1000000,2), ymin=rep(-Inf,2), ymax=rep(Inf,2))
	logP_plt_data<-tidyr::gather(tbl_df(haplo.logP), key="part", value="position", c("start", "end")) %>%
		filter(chr==as.character(tophit_tbl$chr[1]), position>tophit_range$start-(buffer_mb*1000000), position<tophit_range$end+(buffer_mb*1000000)) 
	region_plot<-ggplot() + 
		    geom_rect(data= region, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.3) +
			geom_line(data=logP_plt_data, aes(x=position/1000000, y=logP), colour=cp[1]) + 
			geom_hline(yintercept=threshold, linetype="dashed") +
			theme_minimal() +
			ylab(expression(paste("-log"[10], plain(P)))) +
			xlab("position (mb)") +
			theme(panel.border = element_rect(colour = "black", fill=NA),
	    	panel.grid.minor=element_blank()) +
	    NULL
	founder_effects_plot_tbl<-tophit_tbl %>%
		tidyr::gather(key="founder", value="beta", contains("beta")) %>% 
		mutate(founder=word(founder,1,sep="_")) %>%
		tidyr::gather(key="founderse", value="se", contains("se")) %>% 
		mutate(founderse=word(founderse,1,sep="_")) %>%
		filter(founder==founderse) %>%
		select(c("founder", "beta", "se"))
	founder_effects_plot<-ggplot(founder_effects_plot_tbl, aes(x=founder, y=beta)) + 
    	geom_errorbar(aes(ymin= beta-se,ymax=beta+se), width=.1, colour=cp[1]) +
	    geom_point(colour=cp[1]) +
    	ylab("effect") +
	    theme_minimal() + 
    	theme(axis.title.x=element_blank(), 
    		panel.border = element_rect(colour = "black", fill=NA),
	    	axis.text.x=element_text(angle=45, hjust=1),
	    	panel.grid.minor=element_blank()) +
	    NULL
	plt<-ggarrange(region_plot, founder_effects_plot, ncol=1, nrow=2, heights=c(5/10, 5/10), labels="auto")
	annotate_figure(plt, top = text_grob(paste0(phenotype, " chromosome ", as.character(tophit_tbl$chr[1])), face = "bold"))
}



plot_tophit_haplotypes_raw_pheno<-function(haplo.logP, chr_hitnum, phenotype, phenfile, dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", id.column="IID", threshold, logdrop=2, buffer_mb=10, covar=NULL){
	#get chromosome hitnum order
	tmp<-tbl_df(haplo.logP) %>% group_by(chr) %>% summarise(mx=max(logP)) %>% arrange(desc(mx)) %>% data.frame()
	chr_tophits<-as.character(tmp[,1])
	tophit_tbl<-tbl_df(haplo.logP) %>%
		filter(chr==chr_tophits[chr_hitnum]) %>%
		filter(logP==max(logP))
	lgdrp=max(c(logdrop, 0.1*max(tophit_tbl$logP)))
	tophit_range<-tbl_df(haplo.logP) %>%
		filter(chr== as.character(tophit_tbl$chr[1]), logP>mean(tophit_tbl$logP)-lgdrp) %>%
		summarise(start=min(start), end=max(end))
	region <- data.frame(xmin=rep(tophit_range$start[1]/1000000,2), xmax=rep(tophit_range$end[1]/1000000,2), ymin=rep(-Inf,2), ymax=rep(Inf,2))
	logP_plt_data<-tidyr::gather(tbl_df(haplo.logP), key="part", value="position", c("start", "end")) %>%
		filter(chr==as.character(tophit_tbl$chr[1]), position>tophit_range$start-(buffer_mb*1000000), position<tophit_range$end+(buffer_mb*1000000)) 
	region_plot<-ggplot() + 
		    geom_rect(data= region, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.3) +
			geom_line(data=logP_plt_data, aes(x=position/1000000, y=logP), colour=cp[1]) + 
			theme_minimal() +
			geom_hline(yintercept=threshold, linetype="dashed") +			
			ylab(expression(paste("-log"[10], plain(P)))) +
			xlab("position (mb)") +
			theme(panel.border = element_rect(colour = "black", fill=NA),
	    	panel.grid.minor=element_blank()) +
	    NULL
	### hard code dosages and estimate on untransformed phenotype
	#load dosages and hard call founders
	dosage.file=paste0(dosage.dir, "/", chr_tophits[chr_hitnum], ".hap.dosage.RData")
	load(dosage.file)
	dose<-hap.dosage.data$dosages[as.character(tophit_tbl$locus[1]),,]
	haps<-1*(dose==(apply(dose,1,max))) # hard call which line carries which founder haplotype
	colnames(haps)<-founderNames
	# load phenotypes
	phenos = read.delim(phenfile)
	y = phenos[[phenotype]]
	names(y) = phenos[[id.column]]
	y = y[!is.na(y)] 
	# get phenotype in the same order as hard called dosages
	use.ids = sort(intersect(names(y), dimnames(haps)[[1]]))
	y = y[match(use.ids,names(y),nomatch=0)]
	match.haps = match( use.ids, dimnames(haps)[[1]], nomatch=0)
    haps<-haps[match.haps,]
 	f = lm(y ~ haps -1 ,qr=TRUE)
	beta=coef(f)
	se=sqrt(diag(vcov(f)))
	est_df<-tbl_df(data.frame(founder=founderNames,beta,se))
	#also calculate from dosages to get rsq
	match.dose = match( use.ids, dimnames(dose)[[1]], nomatch=0)
    dose<-dose[match.dose,]
	f2 = lm(y ~ dose ,qr=TRUE)
		#get full model fit if there are covariates
	if(is.null(covar)){
		rsq=data.frame(rsq=summary(f2)$r.squared, rsq.adj=summary(f2)$adj.r.squared, rsq_full=summary(f2)$r.squared, rsq_full.adj=summary(f2)$adj.r.squared)
	} else {
		match.dose.covar = match( use.ids, dimnames(covar)[[1]], nomatch=0)
		covar_MAGIC<-covar[match.dose.covar,]
		f3<-lm(y ~ dose + covar_MAGIC)
		rsq=data.frame(rsq=summary(f2)$r.squared, rsq.adj=summary(f2)$adj.r.squared, rsq_full=summary(f3)$r.squared, rsq_full.adj=summary(f3)$adj.r.squared)
	}
	phenos_plotter<-tbl_df(cbind(y,haps)) %>% tidyr::gather(key="founder", value="indicator", -c(y)) %>% filter(indicator==1) %>% select(-c(indicator))
	founder_effects_plot<-ggplot() +
		geom_jitter(data=phenos_plotter, aes(x=founder, y=y), width=0.1,height=0, colour=cp[1], alpha=0.3, shape=16) +
		geom_errorbar(data=est_df, aes(x=founder, ymin= beta-se,ymax=beta+se), width=.1) +
	    geom_point(data=est_df, aes(x=founder, y=beta)) +
    	ylab(paste0(str_replace_all(phenotype, "_"," "))) +
	    theme_minimal() + 
    	theme(axis.title.x=element_blank(), 
    		panel.border = element_rect(colour = "black", fill=NA),
	    	axis.text.x=element_text(angle=45, hjust=1),
	    	panel.grid.minor=element_blank()) +
	    NULL
	plt<-ggarrange(region_plot, founder_effects_plot, ncol=1, nrow=2, heights=c(5/10, 5/10), labels="auto")
	outplot<-annotate_figure(plt, top = text_grob(paste0(str_replace_all(phenotype, "_"," "), " chromosome ", as.character(tophit_tbl$chr[1])), face = "bold"))
	output_info<-data.frame(phenotype=phenotype, chromosome=as.character(tophit_tbl$chr[1]), min_position=tophit_range$start, max_position=tophit_range$end, peak_start= min(tophit_tbl$start), peak_end=max(tophit_tbl$end), peak_logP=max(tophit_tbl$logP), threshold05=threshold[1], threshold01=threshold[2], rsq, select(tophit_tbl, contains("_beta"))[1,])
	return(list(plot=outplot, info=output_info))
}

### plot founder effects along chromosome
#tidyr::gather(tbl_df(haplo.logP), key="part", value="position", c("start", "end")) %>%
#		filter(chr==as.character(tophit_tbl$chr[1]), position>tophit_range$start-(buffer_mb*1000000), position<tophit_range$end+(buffer_mb*1000000)) %>% select(-contains("_se")) %>%
#		tidyr::gather(key="founder", value="effect", contains("beta")) %>% mutate(founder=str_replace(founder, "_beta","")) %>%
#		ggplot(aes(x=position/1000000, y=effect, group=founder, colour=founder)) + geom_line() + 
#			theme_minimal() +
#			ylim(c(-3,3)) +
#			xlab("position (mb)") +
#			theme(panel.border = element_rect(colour = "black", fill=NA),
#	    	panel.grid.minor=element_blank()) +
#	    NULL

qqplot_df <- function(ps, ci = 0.95) {
  n  <- length(ps)
  data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
}

gg_qqplot<-function(ps, ci=0.95){
df<-qqplot_df(ps, ci=ci)
ggplot(df) +
 geom_ribbon(
     mapping = aes(x = expected, ymin = clower, ymax = cupper),
     alpha = 0.1
   ) +
 geom_point(aes(expected, observed), size = 0.5, colour=cp[1]) +
 geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
 xlab(expression(paste("Expected -log"[10], plain(P)))) +
 ylab(expression(paste("Observed -log"[10], plain(P)))) +
 theme(panel.spacing = unit(0, "lines")) +
 theme(panel.background = element_blank()) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
 theme(legend.position="none") +
 theme(axis.line = element_line(colour = "black")) +
NULL
}

