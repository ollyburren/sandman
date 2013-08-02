#!/usr/bin/Rscript

library(GenomicRanges)

## DON'T CHANGE THESE UNLIKELY TO HAVE MORE THAN 10,000 PERMS
no.perms.per.file<-100
max.perms<-100*100

args<-commandArgs(TRUE)
if(length(args) < 4){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

## example parameters
#perm.dir<-'/stats/oliver/PERMUTATIONS/wtccc'
#snp.file<-'/stats/olly/SANDMAN/FU_QC_STATIC_GW/wtccc/scratch/sigma_snps/chr22.RData'
#get a snp file argument
#n.perms<-10000
#out.dir<-'/stats/oliver/TMP'


assign('snp.gr',get(load(snp.file)))

## function that retrieves and filters perms for a given snp
## object as efficiently as possible.
	
retrieve_perms<-function(snps.gr,perm.dir){
	if(!is(snps.gr,"GRanges"))
		stop("snp.gr parameter must be a GRanges object")
	snps.grl<-split(snps.gr,seqnames(snps.gr))
	## remove chromosomes with no genes
	snps.grl<-snps.grl[which(sapply(snps.grl,length)>0)]
	## remove non standard chromosomes
	snps.grl<-snps.grl[grep("[0-9X]+",names(snps.grl))]
	retval<-lapply(snps.grl,function(x){
		snp.list<-x$name
		chr<-unique(as.character(seqnames(x)))
		permfiles<-list.files(
			path=perm.dir,
			full.name=TRUE,
			pattern=paste(chr,".*RData$",sep=""))
		permfiles<-permfiles[1:(n.perms/no.perms.per.file)]
		##loop over each file retrieving only snps in snps.gr
		p.list<-lapply(seq_along(permfiles),function(i){
			p.file<-basename(permfiles[i])
			assign('p',get(load(permfiles[i])))
			index<-which(rownames(p) %in% snp.list)
			p[index,]
		})
		## we want batches of 1000 (100 x 'by' of 10)
		p.vect<-1:length(p.list)
		c.batches<-split(p.vect,
			cut(p.vect,seq(from=0,to=length(p.vect)+1,by=10)))
		internal<-lapply(seq_along(c.batches),function(g){
			perms<-do.call("cbind",p.list[c.batches[[g]]])
			perms
			#fname<-paste(out.dir,'/',chr,'.perms.',g,'.RData')
			#list(file=fname
		})
		names(internal)<-paste(chr,'.perms.',seq_along(c.batches),'.RData',sep="")
		internal
	})
	names(retval)<-names(snps.grl)
	retval
}



if(n.perms>max.perms)
	stop(paste("Maximum perms are",max.perms))

	if((n.perms %% no.perms.per.file)!=0 | n.perms < no.perms.per.file)
	stop(paste("n.perms must be a factor of and greater than",no.perms.per.file))	
	
	## for permutations we have 100 x 100 (each file contains 1000 permutations)
	## to get 1000 perms we would only need 10 files etc.
	



assign('snp.gr',get(load(snp.file)))	
res<-retrieve_perms(snp.gr,perm.dir)

## next write to a file/files

for(c in names(res)){
	o.list<-res[[c]]
	lapply(seq_along(o.list),function(i){
		perms<-o.list[[i]]
		save(perms,file=paste(out.dir,'/',names(o.list[i]),sep=""))
	})
}

print("Success")
												
	
