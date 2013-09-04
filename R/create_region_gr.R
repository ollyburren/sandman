#!/usr/bin/RScript

args<-commandArgs(TRUE)
if(length(args) < 2){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#region.file='~/GIT_REPOS/sandman/test_region.bed'
#out.dir='~/GIT_REPOS/sandman/'
library(GenomicRanges)

## we check all the args on the perl side as it's easier
## expect a pseudo bed format file where 'name' is the type  
## we allow an optional field of description which allows 
## a name to be given to a region which can be useful
## for debugging (if not supplied then software assigns region.x)

## chr1	123456	456789	test1	<my.region>
                                                            
regions<-read.table(file=region.file,header=F,sep="\t")
if(length(names(regions))==4){
	names(regions)<-c('seqname','start','stop','type')
}else{
	names(regions)<-c('seqname','start','stop','type','names')
}
	
	

regions.list<-split(regions,regions$type)


m.colnames<-paste("SET.",names(regions.list),sep="")


regions.grl<-GRangesList(lapply(seq_along(regions.list),function(i){
	x<-regions.list[[i]]
	gr<-with(x,
		GRanges(seqnames=Rle(seqname),
			ranges=IRanges(start=start,end=stop)
		)
	)
	ncols<-length(m.colnames)
	mcols<-as.data.frame(matrix(logical(length=ncols*length(gr)),ncol=ncols,nrow=length(gr)))
	names(mcols)<-m.colnames
	mcols[[m.colnames[i]]]<-T
	if(length(x$names)>0)
		mcols$names<-x$names
	mcols(gr)<-mcols
	print(gr)
	gr
}))

regions.gr<-unlist(regions.grl)
## not sure that chunksize makes that much sense in this context
## split by chromosome instead
grl<-split(regions.gr,seqnames(regions.gr))

for(i in seq_along(grl)){
	region.gr<-grl[[i]]
	file.out<-paste(out.dir,'/chunk.',i,'.RData',sep="")
	save(region.gr,file=file.out)
}

print("Success")
