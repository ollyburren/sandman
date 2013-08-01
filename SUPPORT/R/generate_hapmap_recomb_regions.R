## create genome hapmap recomb regions for genome - this can then be fed into cal_sigma_1kg.pl
## for use in sandman

library(RMySQL)
library(GenomicRanges)

seqnames<-paste("chr",1:22,sep="")
recomb.thresh=0.1
out.dir='/stats/oliver/1KGenome/RECOMB_0.1_REGIONS/'
uname<-''
dbname<-''
hostname<-''

dbh<-dbConnect(MySQL(),user=uname,dbname=dbname,host=hostname);

calc.intervals.by.recomb<-function(r.df,thresh){
	#must be in positional order
	r.df<-r.df[order(r.df$position),]
	gmap<-r.df$genetic_map_position
	r.df$interval<-trunc(gmap/thresh)+1
	r.df.split<-split(r.df,r.df$interval)
	#here we adjust intervals slightly to use max recomb rate as within an interval this makes sense?
	interval.position<-sapply(r.df.split,function(x) x[which.max(x$recombination_rate),]$position)
	idx<-which(r.df$position %in% interval.position)
	data.frame(start=r.df[head(idx,n=length(idx)-1),]$position+1,end=r.df[idx[-1],]$position)	
}


gr.recomb<-function(dbh,seqname){
	recomb.df<-dbReadTable(dbh,paste(seqname,'_hapmap_recombination_rates',sep=""))
	regions.df<-calc.intervals.by.recomb(recomb.df,recomb.thresh)
	interval.gr<-with(regions.df,
		GRanges(seqnames=Rle(seqname),
						IRanges(start=start,end=end,name=1:nrow(regions.df)))
						)
	split.on<-trunc(as.numeric(names(ranges(interval.gr)))/10)
	sapply(split(interval.gr,split.on),function(x){
		fs<-min(start(x))
		fe<-max(end(x))
		ofile<-paste(out.dir,seqname,'_',fs,'-',fe,'.RData',sep="")
		region.gr<-x
		save(region.gr,file=ofile)
		region.gr
		}
	)
}

for(i in paste("chr",1:22,sep="")){
	gr.recomb(dbh,i)
}


