#!/usr/bin/Rscript
#datafile - RData object containing snpStats object
#output_dir - Directory to where output should be written
#must be set to something (i.e. 0 for no perms)
#nperms - number of permutations
#runnumber
#need to do something clever to sort out wtccc and t1dgc prefixes
args<-commandArgs(TRUE)
if(length(args) < 3){
    cat("Error incorrect number of args","\n",sep="")
    q()
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
   }
}
#datafile<-"/stats/sandman/GENOTYPES/wtccc/wtccc-22.RData"
#output_dir<-"/stats/sandman/PERMUTATIONS/wtccc/"
#runnumber=1
#nperms=2
library(snpStats);
library(wgsea);
objs<-load(datafile)
chr<-gsub("[^\\-]+\\-([0-9]+)\\.RData","\\1",basename(datafile))
dataset<-gsub("([^\\-]+)\\-[0-9]+\\.RData","\\1",basename(datafile))
case<-get(objs[grep("case",objs)])
control<-get(objs[grep("control",objs)])
snpnames<-colnames(case)
perms<-pairtest(case,control,n.perm=nperms)
rownames(perms)<-snpnames
fname<-paste(output_dir,'/',dataset,'/chr',paste(chr,nperms,runnumber,'RData',sep="."),sep="")
save(perms,file=fname)
print("\n")
print(paste("Success",fname))
