## helper script to compute perms for a chromosome.

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


assign('sigma',get(load(sigma.file)))

perms<-compute.mvs.perms(sigma,n.perms)

save(perms,file=out.file)
