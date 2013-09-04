#!/usr/bin/RScript

args<-commandArgs(TRUE)
if(length(args) < 6){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#gene.file='/stats/olly/SANDMAN/FU_RESOURCES/fu_gene_set.tab'
#misc.functions='/home/oliver/GIT_REPOS/sandman/R/miscFunctions.R'
#mart_host='feb2012.archive.ensembl.org'
#tss.extension=200000
#chunksize=100
#out.dir='/stats/oliver/TMP/'


source(misc.functions) #should source libs we need

## we check all the args on the perl side as it's easier



##ENS...	ctrl
genes<-read.table(file=gene.file,header=F,sep="\t")
names(genes)<-c('ensid','type')

gene.list<-split(genes$ensid,genes$type)


## Most performant way is perhaps to get all the SNPs overlapping our genes 
## of interest ? 

mart <- useMart('ENSEMBL_MART_ENSEMBL',host=mart_host)
ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)

gr<-createLogicalGRanges(gene.list,ensembl.gene,tss.extension=tss.extension)

#red.gr<-reduce(gr)
#chunksize<-100
grl<-split(gr,cut(seq_along(gr),seq(0,length(gr)+chunksize,by=chunksize)))


for(i in seq_along(grl)){
	region.gr<-grl[[i]]
	file.out<-paste(out.dir,'/chunk.',i,'.RData',sep="")
	save(region.gr,file=file.out)
}

print("Success")
