library(parallel)
library(abind)

estimate.mixed.model <- function( y, kinship, make.positive=TRUE ) {
    y = y[!is.na(y)]
    if ( length(y) < 1 ) return(NULL)
    
    use.ids = intersect( names(y), colnames(kinship))
    match.kinship = match( use.ids, colnames(kinship), nomatch=0)
    K = kinship[match.kinship,match.kinship]
    K.eigen.trunc = K.eigen = eigen(K,symmetric=TRUE)
    if ( make.positive ) {
        w.eigen = which( K.eigen$values/K.eigen$values[1] > 1.0e-8 )
        neigen = length(w.eigen) # number of non-trivial principal components
        K.eigen.trunc = K.eigen
        K.eigen.trunc$values = K.eigen.trunc$values[w.eigen]
        K.eigen.trunc$vectors = K.eigen.trunc$vectors[,w.eigen]
    }
    match.phen = match( use.ids, names(y), nomatch=0)
    y = y[match.phen]
    y = scale(y)

    z = t(K.eigen.trunc$vectors) %*% y # transformed phenotype
    zz = z*z
    lambda = K.eigen.trunc$values

    opt.theta1 = optimise( interval=c( 0,1 ),
        f =function(theta, zz, lambda ) { # log likelihood for mvn 
            u = theta[1]*lambda+1-theta[1]
            u = ifelse ( abs(u)<1.0e-10, 1.0e-10 , u )
            return(sum(zz/u) + sum( log(u)))
        }, zz=zz, lambda=lambda)

    vg = opt.theta1$minimum[1]
    ve = 1-vg
    cat("estimated heritability", vg,"\n")
    E = K.eigen.trunc$vectors
    Lam = vg*K.eigen.trunc$values +ve
    V = Lam^(-0.5)    
    inv.sqrt.sigma = ( E %*% diag(V) ) %*% t(E)

    mixed.model.data = list( y=y, K=K, vg=vg, ve=ve, inv.sqrt.sigma=inv.sqrt.sigma, eigen=K.eigen.trunc )
    return(mixed.model.data)
}

mixed.model.snps <- function( phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", traw, n.perm=0, covar=NULL) {
		
	y=load_phenotypes(phenfile, phenotype, id.column="IID")


	load(grm.datafile) # loads grm.data
	kinship = grm

	use.ids = sort(intersect(names(y), rownames(kinship)))
	use.ids = intersect(use.ids, rownames(kinship))
	y = y[match(use.ids,names(y),nomatch=0)]
	if(!is.null(covar)){
      	covar<-covar[match(use.ids,rownames(covar),nomatch=0),]
      	y=lm(y~covar)$residuals
    }

    mmd = estimate.mixed.model( y, kinship )

    use.ids = rownames(mmd$y)
    genos<-t(traw[,7:ncol(traw)])
    genos = genos[match( use.ids, rownames(genos), nomatch=0),]
    
    map<-traw[,1:6]
    
    if ( nrow(genos) != length(use.ids)) {
        cat( "ERROR sample ids in genotypes do not match phenotypes\n")
        return(NULL);
    }
    
    genos = apply( genos, 2, function( g ) {
        s = sd(g, na.rm=TRUE)
        mu = mean(g, na.rm=TRUE)
        if ( s > 0 ) {
            g = ifelse ( is.na(g), 0.0, (g-mu)/s )
        } else {
            g = rep(0, length(g));
        }
        return(g)
    })

    mm.transformation = mmd$inv.sqrt.sigma    
    mm.transformed.y = mm.transformation %*% mmd$y
    mm.transformed.g = mm.transformation %*% genos

    mm.gwas.cor =  cor( mm.transformed.y, mm.transformed.g )
    n = length(mm.transformed.y)
    df = n-2
    t.stat = mm.gwas.cor*sqrt((df-2)/(1.0e-10+sqrt(1-mm.gwas.cor*mm.gwas.cor)))
    pval = 2*pt( abs(t.stat), df, low=FALSE )
    logP<-data.frame(map,logP=(-log10(pval[1,])))

    if ( n.perm > 0 ) {
        perms = matrix( NA, ncol=n.perm, nrow=n)
        for( i in 1:n.perm ) {
            perms[,i] = sample(mm.transformed.y,replace=FALSE)
        }
        mm.gwas.cor.perm = cor( perms, mm.transformed.g )
        t.stat.perm = mm.gwas.cor.perm*sqrt((df-2)/(1.0e-10+sqrt(1-mm.gwas.cor.perm*mm.gwas.cor.perm)))
        pval.perm = 2*pt( abs(t.stat.perm), df, low=FALSE)
        logP.perm=t(-log10(pval.perm))
        logP.max = apply( logP.perm, 2, max, na.rm=TRUE)
        logP$perm.adj.pvalue = sapply( 1:nrow(logP), function( i, logP ) {mean(logP$logP[i] < logP.max )}, logP)
        results = list( logP=logP, vg=mmd$vg, max.logP.perm=logP.max)
    } else {
        results = list( logP=logP, vg=mmd$vg)
    }
    return(results)
}

load_phenotypes<-function(phenfile, phenotype, id.column="IID"){
	phenos = read.delim(phenfile)
	if ( is.null(phenos[[phenotype]])) { 
		cat("ERROR phenotype ", phenotype, "missing from ", phenfile, "\n")
		return(NULL)
	}
	if ( is.null(phenos[[id.column]])) {
		cat("ERROR subject id column ", id.column, "missing from ", phenfile, "\n")
		return(NULL)
	}
	y = phenos[[phenotype]]
	names(y) = phenos[[id.column]]
	y = y[!is.na(y)]  
}

load_named_dose<-function(dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chr_name, locus_name){
	dosage.file=paste0(dosage.dir, "/", chr_name, ".hap.dosage.RData")
  	load(dosage.file)
  	dose=hap.dosage.data$dosages[locus_name,,]
  	return(dose)
}

haplotype.dosage.mixed.model <- function( phenfile, phenotype, id.column="IID", grm.datafile="RData/MAGIC.snp.grm.RData", dosage.dir="../../5STITCH/STITCH_haplotypes/LD_pruned_haplotypes", chrs, thin=1,  n.perm=0, mc.cores=1, covar=NULL) {
  
  y=load_phenotypes(phenfile, phenotype, id.column="IID")
 
  load(grm.datafile) # loads grm.data
  kinship = grm
  haplo.logP = NULL
  haplo.logP.perm = NULL
  mmd = NULL

  for ( chr in chrs ) {
    dosage.file=paste0(dosage.dir, "/", chr, ".hap.dosage.RData")
    load(dosage.file)
    hmap = hap.dosage.data$map
    hmap$start=hmap$bp
    hmap$end = c( hmap$bp[2:nrow(hmap)]-1,hmap$bp[nrow(hmap)]+100 )
    sites = 1:nrow(hmap)
    use.sites = sites %% thin == 0
    hmap = hmap[use.sites,]
    
#    hmap = data.frame( chr=rep(chr,length(d$bp)), start=d$bp, end=c(d$bp.delta,100+d$bp.delta[length(d$bp.delta)])+d$bp-1 )
    if ( is.null(mmd)) {
      use.ids = sort(intersect(names(y), dimnames(hap.dosage.data$dosages)[[2]]))
      use.ids = intersect(use.ids, rownames(kinship))
      y = y[match(use.ids,names(y),nomatch=0)]
      if(!is.null(covar)){
      	covar<-covar[match(use.ids,rownames(covar),nomatch=0),]
      	y=lm(y~covar-1)$residuals
      }
      k.ids = match(use.ids,rownames(kinship),nomatch=0)
      kinship = kinship[k.ids,k.ids]
      mmd = estimate.mixed.model( y, kinship )
    }
    dosages=hap.dosage.data$dosages[use.sites,use.ids,]
	results = haplotype.dosage.mixed.model.chr( mmd, dosages, hmap, nperm=n.perm, mc.cores=mc.cores)
    
    haplo.logP = rbind(haplo.logP, cbind(results$logP, results$beta, results$se))
    haplo.logP.perm = rbind(haplo.logP.perm, results$logP.perm)
  }
  
  logP.max=NULL
  if ( n.perm>0) {
    logP.max = apply( haplo.logP.perm, 2, max, na.rm=TRUE)
    haplo.logP$perm.adj.pvalue = sapply( 1:nrow(haplo.logP), 
                              function( i, haplo.logP, haplo.logP.perm ) { 
                                mean(haplo.logP$logP[i] < logP.max )
                              }, haplo.logP, haplo.logP.perm )
  }
  return(list(haplo.logP=haplo.logP, vg=mmd$vg, max.logP.perm=logP.max))
}

haplotype.dosage.mixed.model.chr <- function( mmd, haplotype.dosage, haplo.map, nperm=10, mc.cores=1, covar=NULL) {
  
  # haplotypes is array with dims n.loci*n.individuals*n.founders
  
  mm.transformation = mmd$inv.sqrt.sigma
  y.transformed = mm.transformation %*% mmd$y
    
  n.loci = dim(haplotype.dosage)[1]
  n.id = dim(haplotype.dosage)[2]
  n.founders = dim(haplotype.dosage)[3]
  founders = 1:n.founders
  cat( "n.id=", n.id, " n.loci=", n.loci, " n.founders=", n.founders, "\n")
  if ( n.id==0 ) {
    cat( "ERROR sample ids in genotypes do not match phenotypes\n")
    return(NULL);
  }
  
  if ( nperm > 0 ) {
    perms = matrix( NA, ncol=nperm, nrow=n.id)
    for( i in 1:nperm ) {
      perms[,i] = sample(y.transformed,replace=FALSE)
    }
    #    perms[,1] = y.transformed
    logP.perm = matrix(0,nrow=n.loci,ncol=nperm)
  } else {
    perms = NULL
    logP.perm = NULL
  }
  
  logP = mclapply( 1:n.loci,  function( locus, haplotype.dosage, y.transformed, mm.transformation, perms ) { 
    
    dose = haplotype.dosage[locus,,]
    dose_norm = 2*(dose/apply(dose,1,sum)) #make sure dosage for one MAGIC line sums to exactly 2
    dose.transformed = mm.transformation %*% dose_norm
    
    f = lm(y.transformed ~ dose.transformed -1 ,qr=TRUE)
    a = anova(f)
    logP = -log10(a[1,5])
    if ( !is.null(perms)) {
      ss = apply(perms, 2, function( y.p, f.qr) {
        res=qr.resid(f.qr, y.p)
        rss = sum(res*res)
        fit = y.p-res
        fss = sum(fit*fit)
        return(c(rss,fss))
      }, f$qr )
      rss = ss[1,]
      fss = ss[2,]
      F = (fss/f$rank)/(rss/(n.id-f$rank))
      lp.perm = -pf(F,f$rank,n.id-f$rank,low=FALSE,log.p=TRUE)/log(10)
    } else {
      lp.perm = FALSE
    }
    betas=coef(f)
    names(betas)=paste0("beta",seq(1:n.founders))
    SEs<-sqrt(diag(vcov(f)))
    names(SEs)=paste0("SE",seq(1:n.founders))
    return( c( logP, lp.perm, betas, SEs ))
  }, haplotype.dosage, y.transformed, mm.transformation, perms, mc.cores=mc.cores )
  
  logP = do.call( "rbind", logP )
  if ( nperm > 0 ) {
    logP.perm = logP[,seq(2,nperm+1)]
  } else {
    logP.perm = NULL
  }
  logP.true = data.frame( haplo.map, logP=logP[,1])
  
  results = list( logP=logP.true, logP.perm=logP.perm, beta=logP[, grepl("beta", colnames(logP))], se=logP[,grepl("SE", colnames(logP))])
  return(results)
}

# STITCH IMPORT FUNCTIONS

import.stitch.haplotype.dosages.wrapper <- function( dir) {
  files = list.files(path=dir, pattern=".hap.txt$",full.names=TRUE)
  lapply( files, import.stitch.haplotype.dosages )
}

import.stitch.haplotype.dosages <- function( file.name) {
  
  cat("reading ", file.name, "\n")
  df = read.delim( file.name )

  ex = strsplit(basename(file.name), split="\\.")[[1]]
  b = paste(ex[1:max(1,length(ex)-1)], collapse=".")
  out.file = paste0(dir,"/", b, ".RData")
  out.file2 = paste0(dir, "/", b, ".map.txt")
  
  magic = names(df)[4:ncol(df)]
  locus = paste0(df$chr, "." , df$bp )
  
  loci = unique(locus)
  nloci = length(loci)

  founders = unique(df$K)
  nfounders = length(founders)
  cat("samples=", length(magic), "loci=", nloci, "founders=", nfounders, "out.file=", out.file, "\n")
  a = array(0,dim=c(length(loci),length(magic),length(founders)))
  dimnames(a) = list( loci=loci, magic=magic, founders=founders)
  for( f in 1:nfounders) {
    dd = as.matrix(df[df$K==f,4:ncol(df)])
    a[,,f] = dd
  }
  splat = do.call("rbind", strsplit(loci,".",fixed=TRUE))
  locus.map = data.frame(chr=splat[,1],bp=as.numeric(splat[,2]), locus=loci,file.name=rep(out.file,length(loci)))

  chrs = unique(locus.map$chr)
  hap.dosage.data = list( file=file.name, map=locus.map, chrs=chrs, dosages=a)

  save(hap.dosage.data, file=out.file)
  write.table(locus.map,file=out.file2,sep="\t", quote=FALSE, row=FALSE)
  return(hap.dosage.data)
}

thin.dosages<-function(dosages, threshold=0.5 ) {
  nloci = dim(d)[1]
  d = abs(dosages[1:(nloci-1),,]-dosages[2:nloci,,])
  dd = apply(d,3,sum)
  ddd = apply(dd,2,max)
  w = which(ddd>threshold)
  ww = unique(sort(w,w+1))
  cat("thinning", length(ww)/length(ddd), "\n" )
  return(dosages[ww,,])
}

rename.dosage.arrays <- function( dir,lookup.file=NULL, old.name="Subject.Name", new.name="SUBJECT.NAME"  ) {
  if ( !is.null(lookup.file) & !is.null(old.name) & !is.null(new.name) ) {
    lookup =read.delim(lookup.file)
  } else {
    return(NULL)
  }
  files = list.files(path=dir,pattern=".map.txt")
  maps = lapply(files, read.delim)
  maps = do.call("rbind", maps)
  r = order(maps$chr, maps$bp)
  maps = maps[r,]
  chrs = unique(maps$chr)
  
  for( chr in chrs ) {
   file.name = paste0(dir,"/", chr, ".hap.dosage.RData")
   load(file.name)
      
    if (!is.null(lookup)) {
      browser()
      hap.dosages = hap.dosage.data$dosages
      old.ids = dimnames(hap.dosages)[[2]]
      ids = intersect(old.ids,lookup[[old.name]])
      hap.dosages = hap.dosages[,ids,]
      lookup = lookup[match(ids,lookup[[old.name]],nomatch=0),]
      new.ids = lookup[[new.name]]
      dn = dimnames(hap.dosages)
      dn[[2]] = new.ids
      dimnames(hap.dosages) = dn
      hap.dosage.data$dosages=hap.dosages
      save(hap.dosage.data, file=file.name)
      cat(file.name, length(old.ids), length(new.ids), "\n")
    }
  }
}

combine.dosage.arrays <- function( dir ) {
  files = list.files(path=dir,pattern=".map.txt")
  maps = lapply(paste0(dir,"/",files), read.delim)
  maps = do.call("rbind", maps)
  r = order(maps$chr, maps$bp)
  maps = maps[r,]
  chrs = unique(maps$chr)
  
  for( chr in chrs ) {
    maps.chr = maps[maps==chr,]
    ff = unique(maps.chr$file.name)
    hap.dosages=NULL
    for( f in ff) {
      load(f)
      cat(chr, f, dim(hap.dosage.data$dosages), "\n")
      if ( is.null(hap.dosages) ) {
        hap.dosages=hap.dosage.data$dosages
      } else {
        hap.dosages=abind(hap.dosages,hap.dosage.data$dosages,along=1)
      }
    }
    
    file.name = paste0(dir,"/", chr, ".hap.dosage.RData")
    hap.dosage.data=list(dosages=hap.dosages,chr=chr, map=maps.chr,file=file.name)
    
    save(hap.dosage.data, file=file.name)
    cat(file.name, dim(hap.dosages), "\n")
  }
}

