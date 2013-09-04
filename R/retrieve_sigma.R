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
## testing 02/09/2013
if(full.sigma==0){
	sigma<-retrieve_sigma_fast(index.gr,snps.gr)
}else{
	sigma<-retrieve_sigma_slow(index.gr,snps.gr)
}
save(sigma,file=out.file)
print("Success")
