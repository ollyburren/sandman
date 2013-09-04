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
source(misc.functions) #should source libs we need

region.file
snp.file
excl.file
out.dir
## NOTE: For time being this is set as TRUE in future 
## if no non-autosomal regions/genes are present could 
## set to false so snps2gene does not have extra overhead
## of adding names to each region/gene
add.name<-TRUE

snps<-read.table(file=snp.file,header=F,sep="\t")
names(snps)<-c('name','pval')

excl<-read.table(file=excl.file,header=F,sep="\t")
names(excl)<-c('chr','start','end','name')

assign('regions.gr',get(load(region.file)))
df<-as.data.frame(regions.gr)[,1:3]
df<-df[grep("^chr[0-9XY]+$",df[,1]),]
tabix.param<-paste(paste(df[,1],':',df[,2],'-',df[,3],sep=""),collapse=" ")
tabix.cmd<-paste(tabix.bin,tabix.snp.catalogue.file)
print(paste(tabix.cmd,tabix.param))
t.snps<-read.table(pipe(paste(tabix.cmd,tabix.param)))
#next we filter based on the snps input
sindex<-which(t.snps$V4 %in% snps$name)
snp.loc<-t.snps[sindex,]

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


snps.gr<-snps2gene(snps.gr,regions.gr,add.region.id=add.name)

#next we save the file to our outfile

out.file<-gsub(".RData",".gr.RData",basename(region.file))
out.file<-paste(out.dir,out.file,sep="/")

save(snps.gr,file=out.file)

print("Success")

