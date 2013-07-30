#!/usr/bin/RScript

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
  snp.file='/stats/olly/SANDMAN/test/wtccc/scratch/sigma_snps/chr3.RData'
  gt.dir='/stats/oliver/GENOTYPES/wtccc'
  out.file='/stats/olly/SANDMAN/test/wtccc/sigma/chr3.sigma.RData'
}

assign('snps.gr',get(load(snp.file)))
sigma<-generateSigmaGene(snps.gr,gt.dir)

save(sigma,file=out.file)
