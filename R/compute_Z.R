source("/home/oliver/MVS/R/miscFunctions.R",echo=TRUE)
test=0
if(!test){
  args<-commandArgs(TRUE)
  if(length(args) < 2){
    cat("Error incorrect number of args","\n",sep="")
    q()
  }else{
    for(i in 1:length(args)){
    	eval(parse(text=args[[i]]))
    }
  }
}else{
 #snp.file
 #test.colname
}


library(wgsea)
library(GenomicRanges)
assign('snps.gr',get(load('$snpfile')))
## here we remove anysnps that are missing after doing perms
if(file.exists('$missfile')){
  missing<-scan('$missfile',as.character())
  miss.index<-which(snps.gr\$name %in% missing)
  if(length(miss.index)>1){
    save(snps.gr,file='$snpfile.BACK')
    snps.gr<-snps.gr[-miss.index,]
    save(snps.gr,file='$snpfile');
  } ## else it's already been done !
}
#filter snps based on subset given
snps.gr<-snps.gr[$subset,]
snps.in<-which(snps.gr\$$test_set==TRUE)
W<-wilcoxon(snps.gr\$p.val,snps.in=snps.in)
save(W,file='$outfile');

assign('snps.gr',get(load(snp.file)))
#filter snps based on subset given
snps.gr<-snps.gr[$subset,]
all.snps<-which(rownames(perms) %in% snps.gr$name)
perms<-perms[all.snps,]
snps.in<-which(rownames(perms) %in% snps.df[snps.df[[$test_set]],]$name)
Wstar<-wilcoxon(perms,snps.in=snps.in)
save(Wstar,file='$outfile')
