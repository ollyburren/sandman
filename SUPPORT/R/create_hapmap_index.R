ind.dir<-'/stats/sandman/1KGenome/SIGMA'
out.file<-'/stats/sandman/1KGenome/SIGMA_INDEX/ceu.index.genome.RData'
source("/home/sandman/GIT_REPOS/sandman/R/miscFunctions.R")

index.gr<-create_sigma_dir_index(in.dir)
save(index.gr,file=out.file);
