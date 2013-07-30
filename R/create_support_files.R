#!/usr/bin/RScript
source("/home/oliver/MVS/R/miscFunctions.R",echo=TRUE) #should source libs we need

test=0

## INPUTS
## MISC FUNCTIONS PATH - EVENTUALLY REPLACE
## TABLE OF SNP AND P.VALS
## TABLE OF GENES AND GENE_SETS
## OUTPUT FILE
## TABLE OF EXCLUDED REGIONS (BED FORMAT)

## OUTPUTS NONE SAVES AN RDATA OBJECT TO OUTPUT FILE
if(!test){
  args<-commandArgs(TRUE)
  if(length(args) < 5){
    cat("Error incorrect number of args","\n",sep="")
    q()
  }else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
  }
}else{
  snp.file='/dunwich/scratch/olly/SANDMAN/test_p.vals.tab'
  gene.file='/dunwich/scratch/olly/SANDMAN/test_gene_set.tab'
  excl.file='/dunwich/scratch/olly/SANDMAN/mhc.tab'
  out.file='/stats/olly/SANDMAN/test/wtccc/scratch/snps/snp.gr.RData'
  tabix.bin='~/src/tabix/tabix-0.2.6/tabix'
  ## eventually we will put this on the internet
  tabix.snp.catalogue.file='/stats/oliver/TABIX/dbSNP135.bed.gz'
  chunksize=100
  mart_host<-'feb2012.archive.ensembl.org'
  tss.extension<-200000
}

## we check all the args on the perl side as it's easier

##rs1234 0.01
snps<-read.table(file=snp.file,header=F,sep="\t")
names(snps)<-c('name','pval')

##ENS...	ctrl
genes<-read.table(file=gene.file,header=T,sep="\t")
names(genes)<-c('ensid','type')

gene.list<-split(genes$ensid,genes$type)

excl<-read.table(file=excl.file,header=F,sep="\t")
names(excl)<-c('chr','start','end','name')

## Most performant way is perhaps to get all the SNPs overlapping our genes 
## of interest ? 

mart <- useMart('ENSEMBL_MART_ENSEMBL',host=mart_host)
ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)

gr<-createLogicalGRanges(gene.list,ensembl.gene,tss.extension=tss.extension)

red.gr<-reduce(gr)
#chunksize<-100
grl<-split(red.gr,cut(seq_along(red.gr),seq(0,length(red.gr)+chunksize,by=chunksize)))

#for some reason this is quicker than using the -B option
snp.loc<-do.call("rbind",
	lapply(grl,function(x){
		df<-as.data.frame(x)[,1:3]
		#print(grep("chr[^0-9XY]+$",df[,1]))
		df<-df[grep("^chr[0-9XY]+$",df[,1]),]
		tabix.param<-paste(paste(df[,1],':',df[,2],'-',df[,3],sep=""),collapse=" ")
		tabix.cmd<-paste(tabix.bin,tabix.snp.catalogue.file)
		print(paste(tabix.cmd,tabix.param))
		t.snps<-read.table(pipe(paste(tabix.cmd,tabix.param)))
		#next we filter based on the snps input
		sindex<-which(t.snps$V4 %in% snps$name)
		t.snps[sindex,]
	})
)

## NOTE due to the nature of snp lookup (i.e. by name) if id's change between
## dbsnp builds they will not be taken forward for analysis.

names(snp.loc)<-c('chr','start','end','name')
snps.gr<-with(snp.loc,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),name=name))
#remove duplicate snps
snps.gr<-snps.gr[!duplicated(snps.gr$name),]

## add exclusion flags

excl.gr<-with(excl,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),name=name))

tmp.mcols<-as.data.frame(mcols(snps.gr))
excluded<-lapply(split(excl.gr,excl.gr$name),function(x){
	slist<-as.character(subsetByOverlaps(snps.gr,x)$name)
	print(paste("Processing",unique(as.character(x$name))))
	snps.gr$name %in% slist
})

for(n in names(excluded)){
	tmp.mcols[[n]]<-excluded[[n]]
}
mcols(snps.gr)<-tmp.mcols
	
## next add p.vals
snps.gr$id<-1:length(snps.gr)
new.mcols<-merge(mcols(snps.gr),snps,by.x="name",by.y="name",all.x=T)
mcols(snps.gr)<-new.mcols[order(new.mcols$id),]

snps.gr$name<-as.character(snps.gr$name)

## Here order is important ideally assign control sets before
## test sets perhaps add a new field to gene set file.

snps.gr<-snps2gene(snps.gr,gr,add.gene=TRUE)

#next we save the file to our outfile

save(snps.gr,file=out.file)

