#core functions used for genomic prediction

#Install and load the required library
library(glmnet)
phenotype="HET_2" #Column name of the phenotype to be analysed
id.column="line_name" #Column name with sample names (must match sample names of genotype data)
phenfile<-read.table("DATA/PHENOTYPES/NDM_phenotypes.tsv", header=TRUE, stringsAsFactors=FALSE)
pruned_traw_file<-read.table("DATA/MAGIC_IMPUTED_PRUNED/MAGIC_imputed.pruned.traw", header=TRUE, stringsAsFactors=FALSE) #imputed MAGIC RIL genotypes in traw format, after LD pruning with PLINK option --indep-pairwise 500 10 0.99
#Format datasets for genomic prediction
#Extract phenotype data for the MAGIC lines only
founderNames=c("Banco", "Bersee", "Brigadier", "Copain", "Cordiale", "Flamingo", "Gladiator", "Holdfast", "Kloka", "MarisFundin", "Robigus", "Slejpner", "Soissons", "Spark", "Steadfast", "Stetson")
phenodata <- phenfile[-which(phenfile$line_name %in% founderNames),]
genodata <- as.data.frame(t(pruned_traw_file[,-c(1:6)]))

#impute genotype and phenotype data with the mean. 
NAtomean <- function(x){replace(x, is.na(x), mean(x, na.rm = TRUE))}
imputed.phenodata<-replace(phenodata, TRUE, lapply(phenodata, NAtomean))
imputed.genodata<-data.matrix(replace(genodata, TRUE, lapply(genodata, NAtomean)))

#match the sample names in phenodata and genodata
rownames(imputed.genodata)<-sub("\\_.*", "", rownames(imputed.genodata))
imputed.genodata<-imputed.genodata[imputed.phenodata$line_name,]

#Split data into training and test set
foldid=sample(rep(seq(10),length=454))
subset<-sort(sample(1:504, 50, FALSE))
test.geno<-imputed.genodata[subset,]
train.geno<-imputed.genodata[-subset,]
test.pheno<-imputed.phenodata[subset,]
train.pheno<-imputed.phenodata[-subset,]

#Fit model with cross validation and make predictions
#Fit ridge regression for HET_2
        cv.fit.rr<-cv.glmnet(train.geno, train.pheno$HET_2, foldid=foldid, keep=TRUE, nfolds=10, alpha=0)
#Extract coefficients        
        coef.rr<-predict(cv.fit.rr, s = cv.fit.rr$lambda.min, type = "coefficients")
        fullfit.rr<- coef.rr[1] + (train.geno %*% coef.rr[-1])
        pred.rr<- coef.rr[1] + (test.geno %*% coef.rr[-1])
#correlate predictions to test and training phenotype values
        corr.pred.rr<-cor(test.pheno$HET_2, pred.rr[,1])
        corr.fullfit.rr<-cor(train.pheno$HET_2, fullfit.rr)

#Fit Elastic net        
cv.fit.elnet<-cv.glmnet(train.geno, train.pheno$HET_2, foldid=foldid, keep=TRUE, nfolds=10, alpha=0.5)
        coef.elnet<-predict(cv.fit.elnet, s = cv.fit.elnet$lambda.min, type = "coefficients")
        fullfit.elnet<- coef.elnet[1] + (train.geno %*% coef.elnet[-1])
        pred.elnet<- coef.elnet[1] + (test.geno %*% coef.elnet[-1])
        corr.pred.elnet<-cor(test.pheno$HET_2, pred.elnet[,1])
        corr.fullfit.elnet<-cor(train.pheno$HET_2, fullfit.elnet[,1])
        
#Fit LASSO
cv.fit.lasso<-cv.glmnet(train.geno, train.pheno$HET_2, foldid=foldid, keep=TRUE, nfolds=10, alpha=1)
        coef.lasso<-predict(cv.fit.lasso, s = cv.fit.lasso$lambda.min, type = "coefficients")
        fullfit.lasso<- coef.lasso[1] + (train.geno %*% coef.lasso[-1])
        pred.lasso<- coef.lasso[1] + (test.geno %*% coef.lasso[-1])
        corr.pred.lasso<- cor(test.pheno$HET_2, pred.lasso[,1])
        corr.fullfit.lasso<- cor(train.pheno$HET_2, fullfit.lasso[,1])
  
