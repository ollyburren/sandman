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
assign('index.gr',get(load(sigma.index.file)))
sigma<-retrieve_sigma(index.gr,snps.gr)
save(sigma,file=out.file)
