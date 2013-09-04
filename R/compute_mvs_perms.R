## helper script to compute perms for a chromosome.

args<-commandArgs(TRUE)
if(length(args) < 4){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

source(misc.functions)


assign('sigma',get(load(sigma.file)))
assign('snps.gr',get(load(snp.file)))
## testing 02/09/2013
if(full.sigma==0){
	perms<-compute.mvs.perms_fast(sigma,n.perms)
}else{
	perms<-compute.mvs.perms_slow(sigma,n.perms,snps.gr)
}
save(perms,file=out.file)
print("Success")
