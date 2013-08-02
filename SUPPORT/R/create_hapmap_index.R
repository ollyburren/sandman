ind.dir<-'/stats/oliver/1KGenome/SIGMA'
out.file<-'/stats/oliver/1KGenome/SIGMA_INDEX/ceu.index.genome.RData'
source("/home/oliver/GIT_REPOS/sandman/R/miscFunctions.R")

index.gr<-create_sigma_dir_index(in.dir)
save(index.gr,file=out.file);
