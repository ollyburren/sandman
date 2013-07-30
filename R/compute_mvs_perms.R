## helper script to compute perms for a chromosome.

source("/home/oliver/MVS/R/miscFunctions.R",echo=TRUE)
test=0
if(!test){
  args<-commandArgs(TRUE)
  if(length(args) < 3){
    cat("Error incorrect number of args","\n",sep="")
    q()
  }else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
  }
}else{
	n.perms = 10
	sigma.file = '/dunwich/scratch/olly/SANDMAN/test_sigma/chr13.sigma.RData'
	snp.file = '/dunwich/scratch/olly/SANDMAN/test_space/chr13.RData'
	out.file = '//dunwich/scratch/olly/SANDMAN/test_space/perms/chr13.10.RData'
}

assign('snps.gr',get(load(snp.file)))
assign('sigma',get(load(sigma.file)))

snps.df<-as.data.frame(mcols(snps.gr))

perms<-compute.mvs.perms(sigma,snps.df,n.perms)

save(perms,file=out.file)
