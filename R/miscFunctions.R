##' Calculate a multivariate sample matrix using LD
##' 
##' As an alternative to case/control permutation, sampling multivariate normal 
##' distribution is computationally quicker as well as requiring only sample
##' population genotypes.
##' 
##' @param gt a snpMatrix object
##' @param n The number of permutations to generate from sampling multivariate normal 
##'	default(1000)
##' @return A matrix of chi-square values suitable for wilcoxon testing
##' @author Olly Burren
##' @export
##' @seealso \code{\link{wilcoxon}}
##' @keywords htest
##' @examples
##'
##' data(ld.example)
##' rd<-mvs.perm(ceph.1mb,1000)

library(snpStats)
library(corpcor)
library(mvtnorm)
library(wgsea)


#given a genotype object creates a sigma matrix for use in multivariate sampling

mvs.sigma.ld<-function(gt){
	if(!is(gt,"SnpMatrix"))
		stop("Parameter is not a snpMatrix object!")
	#compute ld - note that for large amount of SNPs this might be computationally
	#expensive
	if((ncol(gt)-1)<=0)
		stop("Cannot calculate sigma for less than 2 snps")	
	ld <- ld(gt, stats=c("R.squared"), depth=ncol(gt)-1,symmetric=T)
	## set missing values to 0
	ld[which(is.na(ld))]<-0;
	## here we attempt different values for diag in the
	## attempt to get a positive definite matrix
	
	## set diag  equal to 1 to have a shot at being positive definite matrix
	diag(ld)<-1
	if(!is.positive.definite(ld,,method="chol")){
		#this recurses through various values of diag if we exceed 1 then
		#we compute the closest matrix that is positive definite.
		ld<-attempt.pos.def(ld)
	}
	ld
}

tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
				     warning = w.handler),
	 warning = W)
}


mvs.perm<-function(sigma,n=1000){
	if(!is.matrix(sigma))
		stop("sigma parameter is not a matrix")		
	if(!is.positive.definite(sigma,,method="chol"))
		stop("sigma is not positive definite")
	
	## in original paper method="chol" was not defined so I assume used eigen default
	## this is slower than the choleski decomp ! Perhaps we should contact the author ?
	rd<-rmvnorm(n,mean=rep(0,ncol(sigma)),sigma=sigma,method="chol")
	# convert to chi.squared and transpose
	#rd<-rd^2
	t(rd)
}

calc.mvs.wstar<-function(snps.in,sigma,n=1000){
				perms<-mvs.perm(sigma,n)
				Wstar<-wilcoxon(perms,snps.in=snps.in)
}


attempt.pos.def<-function(mat,diag.val=1.0001){
  print(paste("diag.val",diag.val))
	if(!is(mat,"Matrix"))
		stop("mat is not a Matrix!")
	if(diag.val >= 1.1){
	  print("Matrix is not positive definite. Finding closest approximation..")
		diag(mat)<-1
		return(as(make.positive.definite(mat),"Matrix"))
	}
	diag(mat)<-diag.val
	if(is.positive.definite(mat,,method="chol")==FALSE){
	  new.diag<-signif(1+((diag.val-trunc(diag.val))*10))
		mat<-attempt.pos.def(mat,new.diag)
	}else{
		return(mat)
	}
}

#helper function that takes a snpMatrix object and prunes SNPs based on r2.threshold
#note that if r2.threshold is a vector of length > 1 then returns a list

prune.r.squared<-function(gt,r2.threshold=0.9,as.snpMatrix=FALSE){
	if(!is(gt,"SnpMatrix"))
		stop("gt parameter is not a snpMatrix object!")
	r2 <- forceSymmetric(ld(gt,depth=ncol(gt)-1,stats="R.squared",symmetric=TRUE))
	r2[which(is.na(r2))]<-0;
	if(as.snpMatrix){
		results<-lapply(r2.threshold,function(x) gt[,prune.snps(r2,x)])
	}else{
		results<-lapply(r2.threshold,function(x) which(colnames(gt) %in% prune.snps(r2,x)))
	}
	names(results)<-paste("r2",r2.threshold,sep=".")
	if(length(results)==1)
	  return(unlist(results[1],use.names=FALSE))
	results
}

prune.snps<-function(r2,thr){
	D<-as.dist(1-r2)
	hc <- hclust(D, method="complete")
	clusters <- cutree(hc, h=1-thr)
	names(clusters[!duplicated(clusters)])
}


## function that takes a list of ensembl id's and creates GRanges objects the parameter tss.extension if
## set extends these regions to incorporate region 5' and 3' to transcriptional start site of gene (note that
## as it stands this is actually the min(start) of the gene so large 5' UTRs will mess things up.

library(biomaRt)
library(GenomicRanges)

convertEnsIDToGranges<-function(glist,ensembl,tss.extension=NULL){
  ## check paramaters
  ## attempt to guess species
  identifier<-table(gsub('[0-9]+$','',glist))
  if(length(identifier)>1)
  	stop("Trying to mix ensembl identifiers")
  external_id<-ifelse(names(identifier)=="ENSG",'hgnc_symbol','external_gene_id')
	gene.details<-getBM(
		filters= c("ensembl_gene_id"),
		attributes= c('ensembl_gene_id','chromosome_name','start_position','end_position','strand',external_id),
		values= list(ensembl_gene_id=glist),
		mart= ensembl)
	## Mitochondrial 'chromosome' is mostly named M
	
	mito.ndx<-which(gene.details$chromosome_name=="MT")
	if(length(mito.ndx)>1)
		gene.details[mito.ndx,]$chromosome_name="M"
		
	names(gene.details)<-gsub(external_id,'external_id',names(gene.details))
	## create GRange object
	gr<-with(gene.details,GRanges(
		seqnames=Rle(paste("chr",chromosome_name,sep ="")),
		ranges=IRanges(start=start_position,end=end_position),
		strand=strand,
		names=ensembl_gene_id,
		external_id=external_id))
	if(!is.null(tss.extension)){
		tmp.gr<-gr
		tss.neg<-which(strand(tmp.gr)=="-")
		tss.pos<-which(strand(tmp.gr)=="+")
		end(tmp.gr[tss.pos])<-ifelse(start(tmp.gr[tss.pos])+tss.extension>end(tmp.gr[tss.pos]),
			start(tmp.gr[tss.pos])+tss.extension,
			end(tmp.gr[tss.pos])
		)
		
		start(tmp.gr[tss.pos])<-start(tmp.gr[tss.pos])-tss.extension
		start(tmp.gr[tss.neg])<-ifelse(end(tmp.gr[tss.neg])-tss.extension < start(tmp.gr[tss.neg]),
			end(tmp.gr[tss.neg])-tss.extension,
			start(tmp.gr[tss.neg])
		)
		end(tmp.gr[tss.neg])<-end(tmp.gr[tss.neg])+tss.extension
		gr<-tmp.gr
	}
	if(length(gr) != length(unique(glist)))
		print("WARNING: Input and output gene lists have different lengths")
	gr
}



createLogicalGRanges<-function(geneset.list,ensembl,tss.extension=NULL){
	if(length(geneset.list)<=1)
		print("Warning you have 1 or less gene sets consider using convertEnsIDtoGRanges instead")
	#first step is to unlist to get the full set
	glist<-unique(unlist(geneset.list))
	gr<-convertEnsIDToGranges(glist,ensembl,tss.extension=tss.extension)
	##assigning to data meta is a pig as cannot use [[
	tmp.df<-as.data.frame(mcols(gr))
	for(i in names(geneset.list)){
		tmp.df[[paste("SET",i,sep=".")]]<-gr$names %in% geneset.list[[i]]
	}
	mcols(gr)<-tmp.df
	gr
}

snps2gene<-function(snp.gr,region.gr,add.region.id=FALSE){
	if(!is(snp.gr,"GRanges") || !is(region.gr,"GRanges"))
		stop("Both parameters are required to be GRanges objects")
	#first step is to assign genenames to each snp
	names.gr<-names(mcols(region.gr))
	tmp.mcols<-mcols(snp.gr)
	for(r in names.gr[grep("^SET\\.",names.gr)]){
		print(paste("Proceesing ",r))
		index<-which(mcols(region.gr)[[r]]==T)
		t.gr<-subsetByOverlaps(snp.gr,region.gr[index,])
		tmp.mcols[[r]]<-tmp.mcols$name %in% t.gr$name
	}
	mcols(snp.gr)<-tmp.mcols
	if(!add.region.id){
		
		return(snp.gr)
	}
	
	## add gene id's 	
	
	## if region.gr$name is not defined we make it up
	
	if(length(region.gr$names)==0){
		print("In here")
		region.gr$names<-paste("region",1:length(region.gr),sep=".")
	}
	ol<-as.data.frame(findOverlaps(snp.gr,region.gr))
	ol$region_id<-region.gr[ol$subjectHits,]$names
	#what to do if a SNP overlaps two sets of genes
	#here we just ignore it and take a first come first served.
	#for some applications this could be problematic. Take the 
	#example of multivariate sampling if a snp overlaps a test and control
	#region then we need to assign to the test region preferentially TODO !!
	nr<-ol[!duplicated(ol$queryHits),]
	##only care about snps that overlap our test regions
	snp.gr<-snp.gr[sort(nr$queryHits),]
	tmp.mcols<-mcols(snp.gr)
	tmp.mcols$region_id<-character(length=length(snp.gr))
	tmp.mcols$region_id<-as.character(nr$region_id)
	mcols(snp.gr)<-tmp.mcols
	snp.gr
}

##
## Code below uses LD to generate sigma's that can be 
## used to generate mvs'
	
generateSigmaGene<-function(snps.gr,gt.filepath){
	if(!is(snps.gr,"GRanges"))
		stop("snp.gr parameter must be a GRanges object")
	snps.grl<-split(snps.gr,seqnames(snps.gr))
	## remove chromosomes with no genes
	snps.grl<-snps.grl[which(sapply(snps.grl,length)>0)]
	## remove non standard chromosomes
	snps.grl<-snps.grl[grep("[0-9X]+",names(snps.grl))]
	## create a list of chr to gt files
	gt.flist<-assignGTToFile(names(snps.grl),gt.filepath)
	names(snps.grl)<-names(gt.flist)
	results<-list()
	for(chr in names(snps.grl)){
		print(paste("Processing ",chr))
		fname<-gt.flist[[chr]]
		objs<-load(fname)
		assign('gt',get(objs[grep('control',objs)]))
		#tmp.grl<-split(snps.grl[[chr]],snps.grl[[chr]]$gene_id)
		tmp.grl<-split(snps.grl[[chr]],snps.grl[[chr]]$region_id)
		gene.results<-lapply(tmp.grl,function(x){
			## snp out just the snps we are interested in
			sindex<-which(colnames(gt) %in% x$name)
			## catch error where there are no
			## snps in the gt file
			if(length(sindex)==0)
				return(NA)
			if(length(sindex)==1){
				snp.name<-colnames(gt[,sindex])
				return(Matrix(1,dimnames = list(snp.name,snp.name)))
			}
			mvs.sigma.ld(gt[,sindex])
		})
		results[[chr]]<-gene.results
	}
	results<-unlist(results)
	names(results)<-gsub("^[^\\.]+\\.","",names(results))
	results
}



assignGTToFile<-function(chrs,path){
	print(chrs)
	if(!file.exists(path))
                stop(paste("Cannot find path parameter dir",path))
	chrs<-gsub("^chr","",chrs)
	chrs<-gsub("X","23",chrs)
	files<-list.files(path=path)
	results<-lapply(chrs,function(x){
		f<-files[grep(paste('\\-',x,'.RData$',sep=""),files)]
		f<-paste(path,f,sep="/")
		if(length(f)==0)
			stop(paste("Cannot find genotype file for chr",x,'in path',path))
		f
	})
	names(results)<-chrs
	results
}

compute.mvs.perms_fast<-function(sig.gene,n.perms){
	perms<-do.call("rbind",lapply(seq_along(sig.gene),function(i){
	  mvs<-matrix()
	  x<-sig.gene[[i]]
	  if(!is(x,"Matrix"))
	  	return()
	  if(length(x)==1){
	  	mvs<-t(as.matrix(exp(-rexp(n.perms))))
	  }else{
	  	withCallingHandlers(
       mvs<-mvs.perm(as.matrix(x),n.perms),
	  			warning=function(w) {
	  				print(paste("WARNING: mvs problem gene:",names(sig.gene)[i],"proably not positive definite trying again"))
	  				x<-attempt.pos.def(x)
	  				mvs.perm(as.matrix(x),n.perms)
	  				#attempt to make positive definite
	  				invokeRestart("muffleWarning")
	  			}
			)
    }
    rownames(mvs)<-rownames(x)
    mvs
  }))
  ## is this needed ?
  perms
}

compute.mvs.perms_slow<-function(sig.gene,n.perms,snps.gr){
	perms<-do.call("rbind",lapply(seq_along(sig.gene),function(i){
	  mvs<-matrix()
	  x<-sig.gene[[i]]
	  if(!is(x,"Matrix"))
	  	return()
	  if(length(x)==1){
	  	mvs<-t(as.matrix(exp(-rexp(n.perms))))
	  }else{
	  	withCallingHandlers(
       mvs<-mvs.perm(as.matrix(x),n.perms),
	  			warning=function(w) {
	  				print(paste("WARNING: mvs problem gene:",names(sig.gene)[i],"proably not positive definite trying again"))
	  				x<-attempt.pos.def(x)
	  				mvs.perm(as.matrix(x),n.perms)
	  				#attempt to make positive definite
	  				invokeRestart("muffleWarning")
	  			}
			)
    }
    rownames(mvs)<-rownames(x)
    ## filter here and return only those SNPs that we are interested in ?
    index<-which(rownames(mvs) %in% snps.gr$name)
    mvs[index,]
  }))
  ## is this needed ?
  perms
}


compute.Wstar<-function(sig.gene,snps.gr,test.colname,ctrl.colname,n.perms){
	snps.df<-as.data.frame(mcols(snps.gr))
	snps.df<-snps.df[snps.df[[test.colname]] | snps.df[[ctrl.colname]],]
	perms<-do.call("rbind",lapply(seq_along(sig.gene),function(i){
		mvs<-matrix()
		x<-sig.gene[[i]]
		if(!is(x,"Matrix"))
			return()
		if(length(x)==1){
			mvs<-t(as.matrix(exp(-rexp(n.perms))))
		}else{
			withCallingHandlers(mvs<-mvs.perm(as.matrix(x),n.perms),
				warning=function(w) {
					print(paste("WARNING: mvs problem gene:",names(sig.gene)[i],"proably not positive definite trying again"))
					x<-attempt.pos.def(x)
					mvs.perm(as.matrix(x),n.perms)
					#attempt to make positive definite
					invokeRestart("muffleWarning")
				}
			)
		}
		rownames(mvs)<-rownames(x)
		mvs
	}))
	## need to filter perms so that only include SNPs 
	## we need
	all.snps<-which(rownames(perms) %in% snps.df$name)
	perms<-perms[all.snps,]	
	snps.in<-which(rownames(perms) %in% snps.df[snps.df[[test.colname]],]$name)
	wilcoxon(perms,snps.in=snps.in)
}



runZ<-function(sigma,snps.gr,test.set,ctrl.set,n.perms){
        #length of list sigma should be the same as snps.gr
        results<-list()
        iterator<-1:length(sigma)
        #filter snps so that only include test and control set snps
        ## NOTE returns a data.frame rather than granges object
        snps.df<-lapply(snps.gr,function(x){
                t<-as.data.frame(mcols(x))
                t[t[[test.set]] | t[[ctrl.set]],]
        })
        for(i in iterator){
                print(paste("Processing ",i))
                t.Z.in<-computeZ.input(sigma[[i]],snps.df[[i]],test.set,ctrl.set,n.perms)
        #       print(list(my.W=t.Z.in$W,my.Wstar=t.Z.in$Wstar))
                results[[i]]<-t.Z.in
        }
        W<-lapply(iterator,function(i) results[[i]]$W)
        Wstar<-lapply(iterator,function(i) results[[i]]$Wstar)
        n.in<-sapply(iterator,function(i) length(results[[i]]$snps.in))
        n.out<-sapply(iterator,function(i) length(results[[i]]$snps.out))

        #print(list(W=W,Wstar=Wstar,n.in=n.in,n.out=n.out))
        Z.value(W=W,Wstar=Wstar,n.in=n.in,n.out=n.out)
}

## this routine takes a special index file and 
## a list of snps and then returns a list of 
## corresponding sigma objects. Note
## that the fast is achieved by prefiltering
## sigma - this makes the downstream permutaion
## generation step significantly faster.

retrieve_sigma_fast<-function(index.gr,snps.gr){
	if(!is(snps.gr,"GRanges"))
		stop("snp.gr parameter must be a GRanges object")
	snps.grl<-split(snps.gr,seqnames(snps.gr))
	## remove chromosomes with no genes
	snps.grl<-snps.grl[which(sapply(snps.grl,length)>0)]
	## remove non standard chromosomes
	snps.grl<-snps.grl[grep("[0-9X]+",names(snps.grl))]
	## narrow down to files
	results<-list()
	for(chr in names(snps.grl)){
		snp.list<-snps.gr[seqnames(snps.gr)==chr,]$name
		ol.files<-subsetByOverlaps(index.gr,snps.gr)$file
		if(length(ol.files)==0){
			print("No sigmas found for this chr found in index")
			return();
		}
		retval<-lapply(ol.files,function(x){
  		assign('sigma',get(load(x)))
  		## we cannot do this as messes up variance estimates 
  		## see retrieve_sigma_no_filter
  		index<-which(rownames(sigma) %in% snp.list)
  		if(length(index)==0 ){
  			#ol.snps<-subsetByOverlaps(snps.gr,index.gr)$name
  			print("WARNING: Expecting some SNPs but none found, list of SNPs missing:")
  			#print(ol.snps)
  			return(NA)
  		}
  		if(length(index)==1){
  			snp.name<-colnames(sigma)[index]
  			return(Matrix(1,dimnames = list(snp.name,snp.name)))
  		}
  			sigma[index,index]
  	})
  	##end for
  	results[[chr]]<-unlist(retval)
  }
  results<-unlist(results)
  results
}

## this version of the routine is much slower 
## but more accurate. Envisage that as a first
## pass run fast version to see if worth running
## slow version

retrieve_sigma_slow<-function(index.gr,snps.gr){
	if(!is(snps.gr,"GRanges"))
		stop("snp.gr parameter must be a GRanges object")
	snps.grl<-split(snps.gr,seqnames(snps.gr))
	## remove chromosomes with no genes
	snps.grl<-snps.grl[which(sapply(snps.grl,length)>0)]
	## remove non standard chromosomes
	snps.grl<-snps.grl[grep("[0-9X]+",names(snps.grl))]
	## narrow down to files
	results<-list()
	for(chr in names(snps.grl)){
		snp.list<-snps.gr[seqnames(snps.gr)==chr,]$name
		ol.files<-subsetByOverlaps(index.gr,snps.gr)$file
		if(length(ol.files)==0){
			print("No sigmas found for this chr found in index")
			return();
		}
		retval<-lapply(ol.files,function(x){
  		assign('sigma',get(load(x)))
  		## we cannot do this as messes up variance estimates 
  		## see retrieve_sigma_no_filter
  		s.length<-length(rownames(sigma))
  		if(s.length==0){
  			#ol.snps<-subsetByOverlaps(snps.gr,index.gr)$name
  			print("WARNING: Expecting some SNPs but none found, list of SNPs missing:")
  			#print(ol.snps)
  			return(NA)
  		}
  		if(s.length==1){
  			snp.name<-colnames(sigma)
  			return(Matrix(1,dimnames = list(snp.name,snp.name)))
  		}
  			sigma
  	})
  	##end for
  	results[[chr]]<-unlist(retval)
  }
  results<-unlist(results)
  results
}


## with prefilter turned on we pull back each of the SNPs that is in our dataset
## and compute the mvs on this. If turned off we compute perms for all SNPs
## in underlying sigma files and then filter afterwards which is much slower.
## I think that the former is the best way forward as we extend to effectively sample
## LD within the recomb block that we have data for and also get the speed benefits
## and permute based on that irrespective of gene boundaries
##

## envisages that you have generated sigma for a 1kgenome population
## see code in this dir for that, and then you have indexed using 
## create_sig_index - the index file can be saved for retrieval

gen_perms_1kg<-function(index.gr,gene.gr,slist,n.perms,prefilter=TRUE){

  ol.files<-subsetByOverlaps(index.gr,gene.gr)$file
  
  ## as things are setup we have to process each block of sigmas independantly 
  ## and then join at the end
  
  	retval <- lapply(ol.files,function(x){
  		assign('t',get(load(x)))
  		sigma.gr<-subsetByOverlaps(t,gene.gr)
  		sigma.list<-lapply(seq_along(sigma.gr),function(i){
  	 		s<-sigma.gr[i]$sigma[[1]]
  	 		index<-which(rownames(s) %in% slist)
  	 		#print(index)
  	 		if(length(index)==0)
  	 			return(NA)
  	 		if(length(index)==1){
         print("Only one snp")
         snp.name<-colnames(s)[index]
         return(Matrix(1,dimnames = list(snp.name,snp.name)))
        }
        ## We need to capture the sigma for all variation in this block
        ## to get acurate sampling of LD 
        ## this will make things slower but I think it is neccessary.
        if(prefilter)
        	s<-s[index,index]
        if(!is.positive.definite(s,,method="chol")){
        	#this recurses through various values of diag if we exceed 1 then
        	#we compute the closest matrix that is positive definite.
        	s<-attempt.pos.def(s)
        }
        ## cut out the snps we want
        if(!prefilter){
        	s[index,index]
        }else{
        	s
        }
      })
      compute.mvs.perms_fast(sigma.list,n.perms)
    })
    do.call("rbind",retval)
}

create_sigma_dir_index<-function(in.dir){
	files<-list.files(path=ind.dir,pattern="*.RData$",full.name=TRUE)
	library(GenomicRanges)
	## parse filenames to make a GRanges object
	tmp<-gsub("\\.RData","",basename(files))

	chr=gsub("^([^_]+).*","\\1",tmp)
	start=as.numeric(gsub(".*_([^\\-]+)\\-.*","\\1",tmp))
	end=as.numeric(gsub(".*_[^\\-]+\\-(.*)","\\1",tmp))

	GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),file=files)
}

	
