#!/usr/bin/RScript

## This is a slave script that should be called by cal_sigma_1kg.pl
## It calculates sigma for multivariate sampling based on LD computed 
## from 1KGenomes 
## Please note that as it stands this means that chrX etc and Y are no
## included.

## Note the hardcoded variables - it's unlikely anyone but developer
## would run this.


library(GenomicRanges)
library(snpStats)
library(VariantAnnotation)

args<-commandArgs(TRUE)
if(length(args)<1){
	cat("Error incorrect number of args","\n",sep="")
	q()
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
	}
}

## need arg region.file  

## note that outfile is hardcoded to be a dir SIGMA

source("~/MVS/R/miscFunctions.R")

thou_gen_data_dir<-'/stats/oliver/1KGenome/VCF/CEU/CEU.'

call.rate<-0.99
z.HWE.co<-25

#region.file<-''

out.file<-gsub("RECOMB_0.1_REGIONS","SIGMA",region.file)

assign('regions.gr',get(load(region.file)))

seqlevels(regions.gr)<-gsub("chr","",seqlevels(regions.gr))

chr<-gsub("^(chr[^_]+).*","\\1",basename(region.file))

vcfname<-paste(thou_gen_data_dir,chr,'.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz',sep="")
vcfp<-ScanVcfParam(which=regions.gr)
thou<-readVcf(vcfname,"hg19",vcfp)
#this code to remove largescale structural deletions which break the snpMatrix conversion code
structural <- logical(length(alt(thou)))
structural[grep("<", unlist(alt(thou)), fixed=TRUE)] <- TRUE
#remove large features as these are indicative of largescale structural issues
structural[width(thou)>100]<-TRUE
thou <- thou[!structural, ]
if(class(alt(thou)) != "DNAStringSetList"){
  alt(thou)<-VariantAnnotation:::.toDNAStringSetList(unlist(alt(thou),use.names=F))
}
rd<-rowData(thou)
### do this next bit on an interval basis
regions.gr$sigma<-lapply(names(ranges(regions.gr)),function(x){
  t.vcf<-thou[rd$paramRangeID == x]
	calls<-geno(t.vcf)$GT
	a0<-ref(t.vcf)
	a1<-alt(t.vcf)
	
	## supress warnings 
	region.snpm<-suppressWarnings(genotypeToSnpMatrix(calls,a0,a1))
	## do some QC 
	## this means that there are no snps in this region
	if(nrow(region.snpm$map)==0 || sum(!region.snpm$map$ignore)==0)
		return(NA)
	sum<-col.summary(region.snpm$genotypes)
	ok.index<-which(with(sum, Call.rate >= call.rate & z.HWE^2 < z.HWE.co ))
	## now work out sigma
	print(paste(x,length(ok.index)))
	if(length(ok.index)==0)
		return(NA)
	if(length(ok.index)==1){
		snp.name<-rownames(sum[ok.index,])
		return(Matrix(1,dimnames = list(snp.name,snp.name)))
	}
	mvs.sigma.ld(region.snpm$genotypes[,ok.index])
})

seqlevels(regions.gr)<-paste("chr",seqlevels(regions.gr),sep="")
save(regions.gr,file=out.file)
print("Success")

