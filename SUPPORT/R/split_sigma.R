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

library(GenomicRanges)
#in.file<-''
#out.dir<-''

assign('t',get(load(in.file)))
fnames<-lapply(seq_along(t),function(y){
		gr<-t[y,]
		sigma<-gr$sigma[[1]]
		chr<-as.character(seqnames(gr))
		start<-start(gr)
		end<-end(gr)
		fname=paste(out.dir,chr,'_',start,'-',end,'.RData',sep="")
		print(fname)
		save(sigma,file=fname)
		fname
})
print(unlist(fnames))
print("Success")

