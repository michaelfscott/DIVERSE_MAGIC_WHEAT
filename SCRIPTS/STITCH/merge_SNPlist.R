args = commandArgs(trailingOnly=TRUE)

SNPlist_path=args[1]
ref_haplotypes_path=args[2]

source("../0Rfunctions/functions.r")

for(chromosome in (1:7)){
	for(subgenome in c("A", "B", "D")){
		pos_1<-read_tsv(paste0(SNPlist_path,"/chr",chromosome, subgenome, "_1.SNPs.txt"), col_names=c("CHR", "POS", "REF", "ALT"))
		pos_2<-read_tsv(paste0(SNPlist_path,"/chr",chromosome, subgenome, "_2.SNPs.txt"), col_names=c("CHR", "POS", "REF", "ALT"))
		pos<-rbind(update_pos(pos_1), update_pos(pos_2))
		pos = mutate_if(pos, is.numeric, as.integer)
		write_tsv(pos, path=paste0(SNPlist_path,"/chr",chromosome, subgenome,".SNPs.txt"), col_names=FALSE)
		colnames(pos)<-c("id", "position", "a0", "a1")
		write_delim(pos, path=paste0(ref_haplotypes_path, "/Founders_combined_chr",chromosome, subgenome,".legend"), delim=" ")
	}
}


