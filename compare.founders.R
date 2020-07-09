library(Rfast)
library(Rcpp)

sourceCpp( code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector matching2segments( NumericVector geno1, NumericVector geno2, double mismatch, double jump ) {
	// mismatch and jump must be negative (eg -3, -10)
	long int i, nsnp = geno1.length();
	NumericMatrix state(2,nsnp+1);
	NumericMatrix score(2,nsnp+1);
	NumericVector path(nsnp);

	if ( mismatch > 0) mismatch = -mismatch;
	if ( jump > 0 ) jump = -jump;
	
	for(i=0;i<nsnp;i++) {
		state(0,i) = state(1,i)=0.0;
		score(0,i) = score(1,i)=0.0;
	}	
	
	score(0,0) = geno1(0)!=geno2(0) ? 1 : mismatch;
	score(1,0) = geno1(0)==geno2(0) ? 1 : mismatch;
	for(i=1;i<nsnp;i++){
		double mscore = geno1(i)==geno2(i) ? 1 : mismatch;
		double matchScore1 = score(1,i-1) + mscore;
		double matchScore2 = score(0,i-1) + mscore + jump;
		if ( matchScore1 > matchScore2) {
			score(1,i) = matchScore1;
			state(1,i) = 1;
		} else {
			score(1,i) = matchScore2;
			state(1,i) = 0;
		}
		mscore = geno1(i)!=geno2(i) ? 1 : mismatch;
		double mismatchScore1 = score(0,i-1) + mscore;
		double mismatchScore2 = score(1,i-1) + mscore + jump;
		if ( mismatchScore1 > mismatchScore2) {
			score(0,i) = mismatchScore1;
			state(0,i) = 0;
		} else {
			score(0,i) = mismatchScore2;
			state(0,i) = 1;
		}
	}
	
	long int k = score(0,nsnp-1) > score(1,nsnp-1) ? 0 : 1;
	for(i=nsnp-1;i>=0;i--){
		path(i) = k;
		k = state(k,i);
	}
	return(path);
}')

founder.mosaics <- function( traw , mismatch, jump) {
	chrs = unique(traw$chromosome)
	colours = c("white", "grey", "black", "red", "orange", "pink", "brown", "lightblue", "blue", "darkblue", "cyan", "lightgreen", "green", "darkgreen", "coral", "yellow")
	mosaics = list()
	for( chr in chrs) {
	
		genos = traw[traw$chromosome==chr,7:ncol(traw)]
		details = traw[traw$chromosome==chr,1:6]
		cat("chr ",chr,nrow(genos),"\n")
#		pdf(paste0(chr,"_founder_mosaic.pdf"),paper="a4r")
		genos = as.matrix(genos)
		nf = ncol(genos)
		
		for(i in 1:ncol(genos)) {
			genos[,i] = as.numeric(genos[,i])
		}
		
		paths = matrix(NA,nrow=nrow(genos),ncol=ncol(genos)*(ncol(genos)-1)/2)		
		k=1;
		for(i in 2:nf) {
			for(j in 1:(i-1)){
				paths[,k] = matching2segments( genos[,i], genos[,j], mismatch, jump)
				k = k+1;
			}
		} 
		
		groups = matrix(NA,nrow=nrow(genos),ncol=nf)
		colnames(groups) = colnames(genos)
		
		for(i in 1:ncol(genos)) {
			groups[,i] = i
		}
		
		k=1;
		for(i in 2:nf) {
			for(j in 1:(i-1)) {
				groups[,i] = ifelse( groups[,i] == i & paths[,k]==1, j, groups[,i])
				k = k+1
			}
		}
		mosaic = data.frame(cbind(data.frame(details),data.frame(groups)))
		mosaics[[chr]]=mosaic

#		par(mar=c(0.1, 4.1, 2, 4.1), xpd=TRUE)
#		par(xpd=TRUE)
#		layout(matrix(c(1,2), 2, 1, byrow = TRUE),widths=c(1), heights=c(5,2))
#		image(groups,col=colours, axes=F,main=paste(chr,nrow(genos)))
#		x = (0:(nf-1))/(nf-1) 
#		axis(2,at=x,labels=colnames(genos),las=2,cex.axis=0.5)
#		box()
#		legend(y=1,x=1.02, legend=rev(colnames(genos)), fill=rev(colours), border="black",  title="Founders",cex=0.5)
		
#		nhap = apply( groups, 1, function(y) { length(unique(y))})
#		x = (1:length(nhap))/length(nhap)
#		plot(x,nhap, ylab="#haps", xlab=chr,main=NULL, t="s",xlim=c(0,1),xaxs="i")

#		dev.off()
		
	}
	save(mosaics,file=paste0("RData/mosaics.mm", mismatch, ".j", jump,".RData"))
}
