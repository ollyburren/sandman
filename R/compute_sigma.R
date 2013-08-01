#!/usr/bin/RScript

args<-commandArgs(TRUE)
if(length(args) < 3){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
source(misc.functions)

assign('snps.gr',get(load(snp.file)))
sigma<-generateSigmaGene(snps.gr,gt.dir)
save(sigma,file=out.file)
