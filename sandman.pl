#!/usr/bin/perl

use strict;

=head1 NAME
	sandman

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHOR
 Olly Burren

=head1 BUGS

=cut

## LIBS CORE

use File::Find;
use Data::Dumper;
use Cwd;
use Config::IniFiles;

## LIBS CUSTOM
use FindBin;
use lib "$FindBin::Bin/lib";
use Sandman;
use Sandman::Step;
use Sandman::Step::Job;
use Sandman::GRIDDriverFactory;

## GLOBALS

use vars(qw/
	%DIRS $GRID_SCRATCH $SNP_CATALOGUE $MART_HOST
	%RSCRIPTS %MAND_PARAM %DEFAULT_PARAM $SANDMAN_ROOT 
  $BASE_DIR %SKIP_STEP $DRIVER $RLIBS $DEBUG
	/);

###########
## SETUP ##
###########

%RSCRIPTS=(
	lookup_snps => 'lookup_snps.R',
	create_gene_files => 'create_gene_gr.R',
	create_region_files => 'create_region_gr.R',
	create_support_files=>'create_support_files.R',
	compute_sigma=>'compute_sigma.R',
	retrieve_sigma=>'retrieve_sigma.R',
	compute_mvs_perms=>'compute_mvs_perms.R',
	retrieve_perms=>'retrieve_perms.R',
  misc_functions=>'miscFunctions.R' ## these need a tidy
);

## Define parameters that have to be set and their types 
## kind of out dated now we use ini files

%MAND_PARAM=(
	file=>{
		GLOBAL=>['geneset_file'],
		DATASET=>['pval_file']
	},
	dir=>{
		GLOBAL=>['base_dir'],
	},
	array=>{
		GLOBAL=>['analysis']
	});

$SANDMAN_ROOT = $ENV{SANDMAN};
if(! $SANDMAN_ROOT){
	warn("SANDMAN env not set using current working directory");
	$SANDMAN_ROOT = getcwd;
}


## these things get filled in for you unless you overwrite them in cfg.

%DEFAULT_PARAM=(
	GLOBAL=>{
		exclude_region_file=>"$SANDMAN_ROOT/default/exclude.bed",
		tss_extension=>200000,
		filter=>'!snps.gr$mhc',
		perm_number=>10000
	});
	

## Relative dir paths for support and temp files

%DIRS=(
	SCRATCH=>'scratch', 
	GR=>'scratch/gen_ranges', 
	TSNP=>'scratch/snps',
	SSNP=>'scratch/split_snps',
	TPERMS=>'scratch/perms',
	SIGMA=>'/scratch/sigma', 
	WSTAR=>'wstar',
	LOG=>'log',
	PERMS=>'perms',
	RESULTS=>'results');

##SET SANDMAN BOOTSTRAP VARIABLES

my $bs_sandman_file = "$SANDMAN_ROOT/.sandman.cnf";

if(!-e $bs_sandman_file){
   $bs_sandman_file = "./",basename($bs_sandman_file);
   die ("Cannot find .sandman.cnf file this is required to set installation variables\n") unless -e $bs_sandman_file;
}

## read in bootstrap config information
my %GLOBALS;
my $sm_cnf = tie %GLOBALS,'Config::IniFiles',(-file=>$bs_sandman_file);

$DEBUG = $GLOBALS{SM}{debug};

$DRIVER = Sandman::GRIDDriverFactory->instantiate(
	$GLOBALS{SM}{grid_driver},
	inifile => $GLOBALS{SM}{grid_config_path});

$RLIBS =  $GLOBALS{SM}{r_libs} || $ENV{R_LIBS};

$SNP_CATALOGUE=$GLOBALS{SM}{default_snp_catalogue}; # probably make this web enabled

## TODO think about rewriting so that biomart is not always need (i.e. keep gene locs in
## bed file ?)
$MART_HOST=$GLOBALS{SM}{default_mart_host}; ## can be overidden by analysis cnf

##SETUP R SCRIPT PATH
foreach my $k(keys %RSCRIPTS){
	my $script_path = "$SANDMAN_ROOT/R/$RSCRIPTS{$k}";
	if(-e $script_path){
		$RSCRIPTS{$k}= $script_path;
	}else{
		debug("Cannot find $script_path check SANDMAN installation",$DEBUG);                 
		exit(1);
	}
}

##

##################
## END OF SETUP ##
##################

########
##MAIN##
########

## replace with formal get options logic such as GetOpt:Long
my $inifile = shift;

my %ds_params;
my $cfg=parse_args($inifile);	 

$BASE_DIR=$cfg->val('GLOBAL','base_dir').'/';

if(-e "$BASE_DIR/results.tab"){
	if(yesno("Already completed this analysis: show results and exit ?")){
		print_results($BASE_DIR);
		debug("## SUCCESS SEE $BASE_DIR",1);
		exit(1);
	}
}

foreach my $ds ($cfg->GroupMembers('DATASET')){
	debug("## ANALYSING $ds",1);
	$ds_params{$ds}=analyse_dataset($cfg,$ds);
}

## final computation of Z scores combines results from both datasets and is separate

my %zparams;
foreach my $ds (keys %ds_params){
	my %a = %{$ds_params{$ds}};
	foreach my $al(keys %a){
		my %p = %{$a{$al}};
		next unless $p{wfile}; ## need this or else cannot compute !
		push @{$zparams{$al}{wstarfile}},$p{wstarfile};
		push @{$zparams{$al}{wfile}},$p{wfile};
	}
}


debug("## COMBINING STUDIES AND COMPUTING Z SCORES",1); 

my $zstep = Sandman::Step->new(
	logdir=>$BASE_DIR."/$DIRS{LOG}/compute_Z",
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$RLIBS"
	);
foreach my $analysis(keys %zparams){
	my $resoutfile = "$BASE_DIR/$analysis.Z.RData";
	compute_Z($zparams{$analysis}{wfile},
		$zparams{$analysis}{wstarfile},
		$resoutfile,$zstep
	);
}
do_step($zstep);
sleep(5); ## let FS catch up
print_results($BASE_DIR);

## if we get here then it's useful to have a copy of any conf files used
## so that we can reproduce analysis at a later date. 

`cp $inifile $BASE_DIR`;
`cp $bs_sandman_file $BASE_DIR`;

debug("## SUCCESS SEE $BASE_DIR",1);

exit(1);

###############
##END OF MAIN##
###############
################
##SUBROUTINES  #
################

sub parse_args{
	my $cnf = shift;
	if(! -e $cnf){
		die ("ERROR: cannot find analysis conf file: $cnf\n");
	}
	my $cfg = new Config::IniFiles( -file => $cnf);
	## set defaults
	foreach my $sect(keys %DEFAULT_PARAM){
		foreach my $key(keys %{$DEFAULT_PARAM{$sect}}){
			$cfg->setval($sect,$key,$DEFAULT_PARAM{$sect}{$key}) unless $DEFAULT_PARAM{$sect}{$key};
		}
	}
	#some checking
	foreach my $k(keys %MAND_PARAM){
		foreach my $sect(keys %{$MAND_PARAM{$k}}){
			my $gfound=0;
			foreach my $g(grep {/^$sect/} $cfg->GroupMembers('DATASET')){
				$gfound=1;
				foreach my $f(@{${MAND_PARAM}{$k}{$sect}}){
					if($k eq 'file'){
						die "Cannot find $k [$g] param $f\n" unless(-e $cfg->val($g,$f));
					}elsif($k eq 'dir'){
						die "Cannot find $k [$g] param $f\n" unless(-d $cfg->val($g,$f));
					}
				}
			}
		}
	}
	## overide defaults if found
	if(my $mh = $cfg->val('GLOBAL','mart_host')){
		$MART_HOST=$mh;
	}	
	return $cfg;
}


=head2 analyse_dataset

 NAME: analyse_dataset
 ARGS: Config::IniFiles object, SCALAR[STRING] 
 FUNCTION: Manages the analysis of a dataset
 RETURNS : A hashref, wfile - R data object to w score, 
 	wstarfile - R data object to permuted w score

=cut

sub analyse_dataset{
	my ($cfg,$dataset)=@_;
	if(my $snp_cat = $cfg->val($dataset,'snp_catalog')){	
  	  if(-e $snp_cat){
	    $SNP_CATALOGUE = $snp_cat
	  }else{
	    debug("Cannot find $snp_cat for this datasource: using $SNP_CATALOGUE",$DEBUG);
	  }
	}
	my %SKIP_STEP=();
	my %FORCE_STEP=(
		support=>0,
		sigma=>0);
		
	(my $DS_NAME = $dataset) =~ s/DATASET\s+//;
	
	##for ease add these to a hash
	my %targs = (
		base_dir=> $cfg->val('GLOBAL','base_dir').'/',
		pval_file => $cfg->val($dataset,'pval_file'),
		datasource_name => $DS_NAME,
		geneset_file => $cfg->val('GLOBAL','geneset_file'),
		regionset_file => $cfg->val('GLOBAL','regionset_file'),
		exclude_region_file => $cfg->val('GLOBAL','exclude_region_file'),
		tss_extension => $cfg->val('GLOBAL','tss_extension'),
		gt_dir=> $cfg->val($dataset,'gt_dir'),
		filter=>$cfg->val('GLOBAL','filter'),
		perm_number=>$cfg->val('GLOBAL','perm_number'),
		sigma_index=>$cfg->val('GLOBAL','sigma_index'),
		perm_dir=>$cfg->val($dataset,'precomp_perm_dir'),
		full_sigma=>$cfg->val('GLOBAL','full_sigma')
		
	);
	
	## set up analysis datastructure                      
	my %analysis;
	foreach my $al($cfg->val('GLOBAL','analysis')){
		my ($tname,$test,$ref)=split(/\s+/,$al);
		$analysis{$tname}=[$test,$ref];
	}
	## STEP SETUP
	
	## if the user has defined perm_dir then skip sigma
	
	$SKIP_STEP{sigma}=1 if $targs{perm_dir};

	my $dirs = make_support_dir($targs{base_dir}.$targs{datasource_name},0);
	
	## STEP CREATE SNP CATALOG
	
	## next make snp support files if required
	
	my $snp_file = $dirs->{'TSNP'}.'/snp.gr.RData';
	my $run_support = 1;
	if(-e $snp_file && !$FORCE_STEP{support}){
		$run_support = 0 unless yesno("support files found do you wish to overwrite?");
	}
	
	
	## 1a create support file if required based on genes
	if($run_support){
		debug("## CREATING SUPPORT FILES",1);
		run_support($dirs,$snp_file,\%targs);
		$FORCE_STEP{sigma}=1;
	}

	## STEP COMPUTE SIGMA'S these are used for the multivariate sampling.
	                                   
	my $rarray;  
	my $sigdotfile = $dirs->{SIGMA}."/.sigma";
	my $run_sigma = 1;
	if(-e $sigdotfile && !$FORCE_STEP{sigma} && !$SKIP_STEP{sigma}){
		$run_sigma = 0 unless yesno("sigma files found do you wish to recompute?");
	}	
	$run_sigma = 0 if $SKIP_STEP{sigma};
	if($run_sigma){
		debug("## COMPUTING SIGMA",1);
		run_sigma($dirs,$snp_file,$sigdotfile,\%targs);
		## automatically recompute perms
		$FORCE_STEP{perms}=1;
	}
	
	
	## STEP COMPUTE PERMS
	
	my $permsdotfile = $dirs->{PERMS}."/.perms";
	my $run_perms = 1;
	if(-e $permsdotfile && !$FORCE_STEP{perms}){
		$run_perms = 0 unless yesno("perm files found do you wish to recompute?");
	}
	if($run_perms){
		debug("## COMPUTING PERMS",1);
		run_perms($dirs,$snp_file,$sigdotfile,$permsdotfile,\%targs);
		## automatically recompute wstar
		$FORCE_STEP{wstar}=1;
		debug("## CHECKING PERMS",1);
		die "Error checking perms " unless
			run_check_perms($dirs,$snp_file,$permsdotfile)
	}   
	
	
	##STEP COMPUTE W FOR PERMS
	debug("## COMPUTING W*",1);
	my $wparams = run_wstar(\%analysis,$dirs,$snp_file,$FORCE_STEP{wstar},\%targs);
	debug("## CONSOLIDATING W* AND COMPUTING W",1);
	return run_consolidate_wstar_and_w($wparams,$dirs,$snp_file,\%targs); 
}


sub run_check_perms{
	my ($dirs,$snp_file,$permsdotfile) = @_;
	my $check_permsdotfile = $dirs->{PERMS}."/.checked";
	if(!-e $check_permsdotfile){
		check_perms($permsdotfile,$snp_file,$dirs->{LOG});
		if(-e $dirs->{PERMS}."/.error"){ ## notsure that this works
			debug("There has been an error generating permutations and they are malformed\n",$DEBUG);
			return 0;
		}elsif(-e $dirs->{PERMS}."/.missing"){
			open(MISS, $dirs->{PERMS}."/.missing") || die "cannot open ".$dirs->{PERMS}."/.missing\n";
			while (<MISS>) {}
			my $miss_count = $.;
			debug("There are $miss_count snps missing from perm files!",$DEBUG);
			## TODO we should attempt here to find these SNPs somehow !
			if(!yesno("Do you wish to continue?")){
				return 0;
			}
		}
	}
	return 1;
}
	

=head2 consolidate_wstar

 NAME: consolidate_wstar
 ARGS: SCALAR[wdir] - path to a dir of permuted wilcoxon statistics.
 	SCALAR[outfile] - path to a location for output file.
 	SCALAR[log_dir] - path to dir of log files
 FUNCTION: Consolidates a set of Wilcoxon statistics split over multiple 
 	files into one file.
 RETURNS : A jobid if program called in queue context.

=cut

sub consolidate_wstar{    
	my ($wdir,$outfile,$step,$log_dir) = @_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $R=<<REND;
wstar.dir<-'$wdir'
out.file<-'$outfile'
flist<-list.files(path=wstar.dir,pattern="^Wstar\\\\.[0-9]+\\\\.RData\$",full.names=TRUE)
Wstar<-unlist(lapply(flist,function(x) get(load(x))))
save(Wstar,file=out.file)
REND
  print $fh $R;
  close($fh);
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  my $job = Sandman::Step::Job->new(command=>$cmd);
  if($step->add_job($job)){
  	return 1;
  }else{
  	return 0;
  } 
}

=head2 compute_Z

 NAME: compute_Z
 ARGS: ARRAYREF[wfiles] - list of paths to wilcoxon statistics for each dataset,
 	ARRAYREF[wstarfiles] - list of paths to permuted wilcoxon statistics for each dataset,
 	SCALAR[outfile] - path to a location for output file.
 	SCALAR[log_dir] - path to dir of log files
 FUNCTION: Computes Z score for a set of datasets given computed wilcoxon and permuted wilcoxon scores
 RETURNS : A jobid if program called in queue context.

=cut

sub compute_Z{
	# need Wstar, W (i.e.pvals) and snpnames in and snpnames out for each dataset
	my ($wfiles,$wstarfiles,$outfile,$step) = @_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $wfname = join ',', map { qq/'$_'/ }@$wfiles;
	my $wstarfname = join ',', map { qq/'$_'/ }@$wstarfiles;
	my $R=<<REND;
library(wgsea)
## note these must be in the same order WRT dataset
wfname<-c($wfname)
wstarfname<-c($wstarfname)
wl<-lapply(wfname,function(x) get(load(x)))
Wstar<-lapply(wstarfname,function(x) get(load(x)))
W<-lapply(wl,function(x) x\$wilcoxon.stat)
n.in<-unlist(lapply(wl,function(x) x\$n.in))
n.out<-unlist(lapply(wl,function(x) x\$n.out))
result<-Z.value(W=W,Wstar=Wstar,n.in=n.in,n.out=n.out)
save(result,file="$outfile")
REND
  print $fh $R;
  close($fh);
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  my $job = Sandman::Step::Job->new(command=>$cmd);
	return($step->add_job($job));
  #return(dispatch_Rscript($cmd,"$logdir/compute_Z"));
  #`$cmd`;                                   
  #unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};    
}
	
=head2 make_support_dir

 NAME: make_support_dir
 ARGS: SCALAR[dir] - dir path to create,
 	SCALAR[force] - If top level dir already exists whether to remove it and recreate (data lost) w/o 
 		interactive prompt.
 FUNCTION: Safely create a set of dirs.
 RETURNS : HASHREF - set of dirs created.

=cut	


sub make_support_dir{
	my ($dir,$force)=@_;
	if(-d $dir && !$force){
		debug("Directory exists and force flag not set.",$DEBUG);
		if(!yesno("Directory $dir already exists are you sure you want to proceed ?")){
			exit(1);
		}
	}elsif(-d $dir && $force){
		debug("Directory exists and force flag set. Deleting directory ..",$DEBUG);
		## perhaps check interactively as this could be incredibly dangerous ?
		remove_tree($dir,{keep_root => 1,result=> \my $del_dirs});
		debug(join("\n",@$del_dirs),$DEBUG);
	}else{
		mkpath_wrapper($dir);
	}
	my %reshash;
	foreach my $k(keys %DIRS){
		if(mkpath_wrapper("$dir/".$DIRS{$k})){
			$reshash{$k}="$dir/".$DIRS{$k};
		}
	}
	return(\%reshash);
}

=head2 create_gene_files

 NAME: create_gene_files
 ARGS:  	SCALAR[genefile] - path to input gene file see docs for format
 	SCALAR[outdir] - path to where output should be written
 	SCALAR[tabix_chunksize] - Integer size (# of genes) for chunking SNP retrieval
 	SCALAR[mart_host] - biomart host to use see defaults.
 	SCALAR[tss_extension] - Integer for extension for region surrounding  5' TSS see docs.
 	SCALAR[log_dir] - path to dir of log files
 FUNCTION: Safely create a set of dirs.
 RETURNS : HASHREF - set of dirs created.

=cut	

sub create_gene_files{
	my ($genefile,$outdir,$tabix_chunksize,
			$mart_host,$tss_extension,$log_dir)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/create_gene_files",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{create_gene_files}  misc.functions=\\'$RSCRIPTS{misc_functions}\\' ";
	$cmd.= "gene.file=\\'$genefile\\' ";
	$cmd.= "out.dir=\\'$outdir\\' ";
	$cmd.= "chunksize=$tabix_chunksize ";
	$cmd.= "mart_host=\\'$mart_host\\' tss.extension=$tss_extension";
	my $job = Sandman::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	return do_step($step);
}

sub create_region_files{
	my ($regionfile,$outdir,$log_dir)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/create_region_files",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{create_region_files} ";
	$cmd.= "region.file=\\'$regionfile\\' ";
	$cmd.= "out.dir=\\'$outdir\\' ";
	my $job = Sandman::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	return do_step($step);
}

sub do_step{
	my $step = shift;
	if($step->execute()){
		$step->wait_on_complete();
		return 1;
	}else{
		return 0;
	}
}

=head2 make_support_files

 NAME: make_support_files [DEPRECATED]
 ARGS: SCALAR[snpfile] - path to GRanges formatted snpfile,
 	SCALAR[genefile] - path to input gene file see docs for format
 	SCALAR[outfile] - path to where output should be written
 	SCALAR[exclfile] - path to list of regions to exclude see docs for format
 	SCALAR[tabix] - path to tabix binary see docs
 	SCALAR[snpcat] - path to tabix indexed snp catalogue file
 	SCALAR[tabix_chunksize] - Integer size (# of genes) for chunking SNP retrieval
 	SCALAR[mart_host] - biomart host to use see defaults.
 	SCALAR[tss_extension] - Integer for extension for region surrounding  5' TSS see docs.
 	SCALAR[log_dir] - path to dir of log files
 FUNCTION: Given a set of gene lists (see docs for format), routine uses biomart
 	to retrieve genomic intervals, these can be extended using the tss_extension parameter.
 	A catalogue of SNPs is then interrogated using tabix to obtain a set of SNPs that overlap
 	these intervals, information is then added as to which set each SNP belongs. This forms the
 	basis for future analysis.
 RETURNS : SCALAR job_id if run in the context of queuing software.
 
 NOTE: Currently this routine is deprecated as functionality is replaced by create_gene_files 
 	and create_region_files

=cut	

sub make_support_files{
	my ($snpfile,$genefile,$outfile,
			$exclfile,$tabix,$snpcat,
			$tabix_chunksize,$mart_host,$tss_extension,$log_dir)=@_;
	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{create_support_files} snp.file=\\'$snpfile\\' misc.functions=\\'$RSCRIPTS{misc_functions}\\' ";
	$cmd.= "gene.file=\\'$genefile\\' excl.file=\\'$exclfile\\' ";
	$cmd.= "out.file=\\'$outfile\\' tabix.bin=\\'$tabix\\' ";
	$cmd.= "tabix.snp.catalogue.file=\\'$snpcat\\' chunksize=$tabix_chunksize ";
	$cmd.= "mart_host=\\'$mart_host\\' tss.extension=$tss_extension";
	return(dispatch_Rscript($cmd,"$log_dir/create_support_files"));
}

sub run_lookup_snps{
	my ($snpfile,$regiondir,$outdir,
		$exclfile,$tabix,$snpcat,$log_dir)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/lookup_snps",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my @jids;
	find(sub{
			if(/chunk\.[0-9]+\.RData/){
				my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{lookup_snps} misc.functions=\\'$RSCRIPTS{misc_functions}\\' ";
				$cmd.= "region.file=\\'$File::Find::name\\' excl.file=\\'$exclfile\\' ";
				$cmd.= "tabix.bin=\\'$tabix\\' ";
				$cmd.= "tabix.snp.catalogue.file=\\'$snpcat\\' ";
				$cmd.= "snp.file=\\'$snpfile\\' out.dir=\\'$outdir\\'";
				my $job = Sandman::Step::Job->new(command=>$cmd);
				$step->add_job($job);
			}
	},$regiondir);
	return do_step($step);
}

sub check_perms{
  my ($dotfile,$snpfile,$log_dir)=@_;
  	my $step = Sandman::Step->new(
		logdir=>"$log_dir/check_perms",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
  my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
  my $R=<<REND;
dotfile<-'$dotfile'
perm.files<-scan(dotfile,as.character())
library(GenomicRanges)
assign('snps.gr',get(load('$snpfile')))
assign('perm',get(load(perm.files[1])))
sindex<-which(snps.gr\$name \%in\%  rownames(perm))
if(length(sindex) != length(snps.gr)){
  missing<-snps.gr[-sindex,]\$name
  missing.fname<-paste(dirname(perm.files[1]),"/.missing",sep="")      
  write(missing,file=missing.fname)
}
check.fname<-paste(dirname(perm.files[1]),"/.checked",sep="")
file.create(check.fname)
REND
  print $fh $R;
  close($fh);
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  my $job = Sandman::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	return do_step($step);
}

sub consolidate_snps{
	my ($snpdir,$snpfile,$log_dir)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/consolidate_snps",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $R=<<REND;
library(GenomicRanges)
out.file<-'$snpfile'
in.dir<-'$snpdir'
sfiles<-list.files(path=in.dir,pattern="chunk.*",full.names=T)
library(GenomicRanges)
snps.gr<-unlist(GRangesList(lapply(sfiles,function(x){
		get(load(x))
		})
))
snps.gr<-snps.gr[!duplicated(snps.gr\$name),]
save(snps.gr,file=out.file)
REND
  print $fh $R;
  close($fh);
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
	my $job = Sandman::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	return do_step($step);
}

       

## helper function that takes a support file and splits it by chromosome
## as then we can compute sigma's on the q much more quickly

sub prepare_snps_for_q{
  my ($outdir,$snpfile,$log_dir)=@_;
  my $step = Sandman::Step->new(
		logdir=>"$log_dir/prepare_snps_for_q",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
  my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
  #open(TMP,">$fname") || die "Cannot open tmp file for writing\n";
  my $R=<<REND;
library(GenomicRanges)
out.dir<-'$outdir/'
assign('master.snps.gr',get(load("$snpfile")))
grl<-split(master.snps.gr,seqnames(master.snps.gr))
sapply(seq_along(grl),function(i){
 fname<-paste(out.dir,names(grl)[i],".RData",sep="")
 snps.gr<-grl[[i]]
 save(snps.gr,file=fname)
})
REND
  print $fh $R;
  close($fh);
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  my $job = Sandman::Step::Job->new(command=>$cmd);
  $step->add_job($job);
	if(do_step($step)){
  	return 1;
  }else{
  	return 0;
  }
}

sub create_collapse_perms_script{
	my ($flist,$outfile)=@_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
  my $R=<<REND;
fl=c($flist)
perms<-do.call("rbind",lapply(fl,function(x){
	get(load(x))
}))
save(perms,file='$outfile')
REND
	#print $R;
	print $fh $R;
	close($fh);
	return $fname;
}

sub collapse_perms{
	my($indir,$outdir,$log_dir)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/collapse_perms",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my %hash;
	find(sub{
			if(/\.perms\.([0-9]+)\.RData$/){
				push @{$hash{$1}},$File::Find::name;
			}
	},$indir);
	my @jids;
	my @tscripts;
	my @outfiles;
	foreach my $p(keys %hash){
		my $flist = join ',', map { qq/'$_'/ }@{$hash{$p}};
		my $outfile = $outdir."/$p.allperms.RData";
		my $fname=create_collapse_perms_script($flist,$outfile);
		push @outfiles,$outfile;
		my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
		my $job = Sandman::Step::Job->new(command=>$cmd);
		$step->add_job($job);
		push @tscripts,$fname;
	}
	if(do_step($step)){
		unlink @tscripts;
		return \@outfiles;
	}else{
		die("Unable to collapse perms");
	}
}
	
 
## TODO unlink tmp files 

sub retrieve_sigmas{
  my ($snp_file,$sigma_dir,$sig_index,$scratch_dir,$chrXgt_dir,$full_sigma,$log_dir)=@_;
  $full_sigma = 1 if $full_sigma; ## incase they have used true or false
  my @return;
	#we blow away contents of scratch_dir each time
	remove_tree($scratch_dir,{keep_root => 1,result=> \my $del_dirs});
	debug(join("\n",@$del_dirs),$DEBUG);
	
	my $outcome = prepare_snps_for_q($scratch_dir,$snp_file,$log_dir);
	die("could not prepare snps for the q") unless $outcome;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/retrieve_sigma",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	find(sub{
			debug($_,$DEBUG);
      if(/\.RData$/){
      	my $out_file = $sigma_dir.'/'.basename($_,'.RData').".sigma.RData";
        my $jid;
        if(/chrX\.RData$/){
        	if(! -e $chrXgt_dir){
        		my $prompt=<<PROMPT;
There are regions/genes in the analysis on chrX,but 1kg_gt_dir is either not set or does not exist	
in experiment conf file. Type [y]es if you wish to continue with the analysis and ignore these regions.
PROMPT
						if(!yesno($prompt)){
							exit(1);
						}
					}
        	## this means we need to calculate some sigma's rather than retrieve
        	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_sigma} snp.file=\\'$File::Find::name\\' ";
        	$cmd.= "gt.dir=\\'$chrXgt_dir\\' out.file=\\'$out_file\\' misc.functions=\\'$RSCRIPTS{misc_functions}\\'";
        	my $job = Sandman::Step::Job->new(command=>$cmd);
        	$step->add_job($job);
        	#$jid = dispatch_Rscript($cmd,"$log_dir/calculate_sigma");
        }else{
        	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{retrieve_sigma} snp.file=\\'$File::Find::name\\' ";
        	$cmd.= "sigma.index.file=\\'$sig_index\\' out.file=\\'$out_file\\' misc.functions=\\'$RSCRIPTS{misc_functions}\\' ";
        	$cmd.= "full.sigma=$full_sigma";
        	my $job = Sandman::Step::Job->new(command=>$cmd);
        	$step->add_job($job);
        	#$jid = dispatch_Rscript($cmd,"$log_dir/retrieve_sigma");
        }
        push @return,{
        	snpfile=>$File::Find::name,
          sigmafile=>$out_file
        };
      }
	},$scratch_dir);
  if(do_step($step)){
  	return \@return;
  }else{
  	die("Unable to retrieve sigmas");
  }
}

  

sub create_sigmas{
  my ($snp_file,$sigma_dir,$gt_dir,$scratch_dir,$log_dir)=@_;
  #create a temp dir
    my $step = Sandman::Step->new(
		logdir=>"$log_dir/compute_sigma",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
 #we blow away contents of scratch_dir each time
	remove_tree($scratch_dir,{keep_root => 1,result=> \my $del_dirs});
	debug(join("\n",@$del_dirs),$DEBUG);
	my $jid = prepare_snps_for_q($scratch_dir,$snp_file,$log_dir);
	do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished([$jid]);
	my @return;
	find(sub{
 		 debug($_,$DEBUG);
 		 if(/\.RData$/){
 		 	 my $out_file = $sigma_dir.'/'.basename($_,'.RData').".sigma.RData";
 		 	 my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_sigma} snp.file=\\'$File::Find::name\\' ";
 		 	 $cmd.= "gt.dir=\\'$gt_dir\\' out.file=\\'$out_file\\' misc.functions=\\'$RSCRIPTS{misc_functions}\\'";
 		 	 my $job = Sandman::Step::Job->new(command=>$cmd);
 		 	 $step->add_job($job);	
 		 	 push @return,{
		  	snpfile=>$File::Find::name,
		    sigmafile=>$out_file
		  };
		}
	},$scratch_dir);
	if(do_step($step)){
		return \@return;
	}else{
		die("Unable to create sigmas");
	}
}

## compute W

sub print_results{
	my ($base_dir)=@_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $R=<<REND;
result.files<-list.files(path="$base_dir",pattern="*.RData",full.name=TRUE)
rlist<-lapply(result.files,function(x){
  get(load(x))
})
names(rlist)<-gsub(".Z.RData\$","",basename(result.files))
res<-sapply(rlist,function(x){
  list(Z.emp=x\$Z.empirical\$statistic,emp.p.val=x\$Z.empirical\$p.value,
  Z.theo=x\$Z.theoretical\$statistic,theo.p.val=x\$Z.theoretical\$p.value)
})
write.table(res,file='$base_dir/results.tab',quote=FALSE,sep="\\t")
print(res)
REND
	print $fh $R;
  close($fh);
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  print `$cmd`; 
  unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};
}

sub compute_w{
	my ($snpfile,$outfile,$test_set,$ctrl_set,$filter,$missfile,$step,$logdir)=@_;
	my $subset = "(snps.gr\$$test_set | snps.gr\$$ctrl_set) & $filter";
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $R=<<REND;
library(wgsea)
library(GenomicRanges)
assign('snps.gr',get(load('$snpfile')))
if(file.exists('$missfile')){
  missing<-scan('$missfile',as.character())
  miss.index<-which(snps.gr\$name %in% missing)
  snps.gr<-snps.gr[-miss.index,]
}
## here snpfile should have had those snps that are missing in 
## perms removed
#filter snps based on subset given
snps.gr<-snps.gr[$subset,]
snps.in<-which(snps.gr\$$test_set==TRUE)
wilcoxon<-wilcoxon(snps.gr\$pval,snps.in=snps.in)
## we also need to save snp.in and snp.out values 
snps.out<-setdiff(which(snps.gr\$$ctrl_set==TRUE),snps.in)
W<-list(wilcoxon.stat=wilcoxon,n.in=length(snps.in),n.out=length(snps.out))
save(W,file='$outfile');
REND
  print $fh $R;
  close($fh);
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
	my $job = Sandman::Step::Job->new(command=>$cmd);
	if($step->add_job($job)){
		return 1;
	}else{
		return 0;
	}
} 


sub retrieve_perms{
	my($snp_file,$scratch_dir,$perm_in_dir,$perm_number,$perm_tmp,$log_dir)=@_;
		my $step = Sandman::Step->new(
		logdir=>"$log_dir/retrieve_perms",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
  remove_tree($scratch_dir,{keep_root => 1,result=> \my $del_dirs});
  debug(join("\n",@$del_dirs),$DEBUG);
  my $jid = prepare_snps_for_q($scratch_dir,$snp_file,$log_dir);
  do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished([$jid]);
  find(sub{
  		debug($_,$DEBUG);
  		if(/\.RData$/){
  			my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{retrieve_perms} snp.file=\\'$File::Find::name\\' ";
  			$cmd.="perm.dir=\\'$perm_in_dir\\' out.dir=\\'$perm_tmp\\' n.perms=$perm_number";
  			my $job = Sandman::Step::Job->new(command=>$cmd);
  			$step->add_job($job);
  		}
  },$scratch_dir);
  return do_step($step);
}

## for a given snp file and matching sigma file
## create a matrix of perms that we can use to
## compute Wstar
	
sub compute_perms{
	my ($plist,$permdir,$number_perms,$full_sigma,$log_dir,$chunksize)=@_;
	my $step = Sandman::Step->new(
		logdir=>"$log_dir/compute_perms",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	$full_sigma=1 if $full_sigma;
	$chunksize=1 unless $chunksize; # in case it's not set;
  foreach my $phash(@$plist){
  	my $snpfile = $phash->{snpfile};
  	my $sigmafile = $phash->{sigmafile};
  	for(my $x=1;$x<=$chunksize;$x++){
  		my $out_file = $permdir.'/'.basename($snpfile,'.RData').".perms.$x.RData";
    	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_mvs_perms} snp.file=\\'$snpfile\\' ";
    	$cmd.= "sigma.file=\\'$sigmafile\\' out.file=\\'$out_file\\' n.perms=$number_perms misc.functions=\\'$RSCRIPTS{misc_functions}\\' ";
    	$cmd.= "full.sigma=$full_sigma";
    	#debug($cmd);  
    	my $job = Sandman::Step::Job->new(command=>$cmd);
    	$step->add_job($job);
    }
  }
  return do_step($step);
}


sub create_compute_wstar_script{
	my ($test_set,$ctrl_set,$filter,$snpfile,$permfile,$missfile,$outfile)=@_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
  #open(TMP,">$fname") || die "Cannot open tmp
  #my $filter = '!snps.gr$MHC'; # test to exclude MHC
  my $subset = "(snps.gr\$$test_set | snps.gr\$$ctrl_set) & $filter";
  my $R=<<REND;
library(wgsea)
library(GenomicRanges)
assign('perms',get(load('$permfile')))
assign('snps.gr',get(load('$snpfile')))
## here we remove anysnps that are missing after doing perms
if(file.exists('$missfile')){
  missing<-scan('$missfile',as.character())
  miss.index<-which(snps.gr\$name %in% missing)
  snps.gr<-snps.gr[-miss.index,]
  ## this code below cause race conditions !!
  ## do not use !
  #if(length(miss.index)>1){
  #  save(snps.gr,file='$snpfile.BACK')
  #  snps.gr<-snps.gr[-miss.index,]
  #  save(snps.gr,file='$snpfile');
  #} ## else it's already been done !
}
#filter snps based on subset given
snps.gr<-snps.gr[$subset,]
all.snps<-which(rownames(perms) \%in\% snps.gr\$name)
perms<-perms[all.snps,]
snps.in<-which(rownames(perms) \%in\% snps.gr[snps.gr\$$test_set,]\$name)
Wstar<-wilcoxon(perms,snps.in=snps.in)
save(Wstar,file='$outfile')  
REND
	print $fh $R;
	close($fh);
	return $fname;
}

sub compute_wstar{
	my($test_set,$ctrl_set,$filter,$snpfile,$permdir,$wstardir,$missfile,$step,$logdir)=@_;
	my $outdir = "$wstardir/${test_set}_VS_$ctrl_set";
	mkpath_wrapper($outdir);
	my $retval = 0;
	my @jobs;
	find(sub{
			if(/([0-9]+)\.allperms\.RData$/){
				my $outfile = "$outdir/Wstar.$1.RData" ;       
				my $fname = create_compute_wstar_script($test_set,$ctrl_set,$filter,$snpfile,$File::Find::name,$missfile,$outfile);
				my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
				my $job = Sandman::Step::Job->new(command=>$cmd);
				$step->add_job($job);
				$retval=1;
			}
	},$permdir);
	return $retval;
}

## ANALYSIS STEPS

sub run_support{
	my($dirs,$snp_file,$args)=@_;
	if(-e $args->{'geneset_file'}){
		debug("FOUND GENESET creating GenomicRanges",$DEBUG);
		my $outcome = create_gene_files(
			$args->{'geneset_file'},
			$dirs->{GR},
			$GLOBALS{SM}{tabix_chunksize},
			$MART_HOST,$args->{tss_extension},
			$dirs->{LOG});
		die("Could not create GenomicRange") unless $outcome;
	}elsif(-e $args->{'regionset_file'}){ ## 1b or based on regions
		debug("FOUND REGIONSET creating GenomicRanges",$DEBUG);
			my $outcome = create_region_files(
				$args->{'regionset_file'},
				$dirs->{GR},
				$dirs->{LOG});
			die("Could not create GenomicRange") unless $outcome;
	}
	
	## using defined catalog assign positional locations for snps
	## in submitted pvalue file.
	
	my $outcome = run_lookup_snps(
		$args->{'pval_file'},
		$dirs->{GR},
		$dirs->{TSNP},
		$args->{'exclude_region_file'},
		$GLOBALS{SM}{tabix_bin},$SNP_CATALOGUE,
		$dirs->{LOG});
	die("Could not lookup SNPs") unless $outcome;
	$outcome=0;
	$outcome = consolidate_snps($dirs->{TSNP},$snp_file,$dirs->{LOG});
	die("Could not consolidate snps") unless $outcome;
}

sub run_sigma{
	my($dirs,$snp_file,$sigdotfile,$args)=@_;
	my $rarray;
	if($args->{sigma_index}){
	##retrieve sigmas
	$rarray = retrieve_sigmas(
		$snp_file,
		$dirs->{SIGMA},
		$args->{sigma_index},                                                   
		$dirs->{SSNP},
		$args->{gt_dir},
		$args->{full_sigma},
		$dirs->{LOG});
	}else{
		if(!-d $args->{gt_dir}){
			die("Cannot find genotype dir ".$args->{gt_dir}.": If you haven't got genotypes consider using precomputed sigmas (sigma_index parameter)\n.")
		}
		$rarray = create_sigmas(
			$snp_file,
			$dirs->{SIGMA},
			$args->{gt_dir},
			$dirs->{SSNP},
			$dirs->{LOG});
    }
		#write the snpfile to sigma mapping to a file
		#doubles as evidence that we have done things and also 
		#means we can run subsequent commands w/o having to run this
		open(SIG,">$sigdotfile") || die "Cannot open sigma dot file: $sigdotfile\n";
		foreach my $h(@$rarray){
			print SIG join("\t",$h->{snpfile},$h->{sigmafile})."\n";
		}
		close(SIG);
		if(!-z $sigdotfile){
			return 1;
		}else{
			return 0;
		}
}

sub run_perms{
		my($dirs,$snp_file,$sigdotfile,$permsdotfile,$args)=@_;
		if(-d $dirs->{PERMS}){
			## clear up first
			remove_tree(($dirs->{PERMS},$dirs->{TPERMS}),{keep_root => 1,result=> \my $del_dirs});
			debug(join("\n",@$del_dirs),$DEBUG);
		}
		## read in sigma .file
		my @jids;
		if($args->{perm_dir}){
			## use precomputed perms i.e from permuting genotypes in case/control
			my $outcome = retrieve_perms($snp_file,
				$dirs->{SSNP},
				$args->{perm_dir},
				$args->{perm_number},
				$dirs->{TPERMS},
				$dirs->{LOG});
			die("Unable to retrieve perms \n") unless $outcome;
		}else{ ## generate perms 
			open(SIG, $sigdotfile) || die "Cannot find sigma dot file $sigdotfile\n";
			my @plist;
			while(<SIG>){
				chomp;
				my($sfile,$sigfile) = split("\t",$_);
				my %hash = (snpfile=>$sfile,
										sigmafile=>$sigfile);
				push @plist,\%hash;
			} 
			close(SIG);
			## clear up
		
	
			my $csize = ceil($args->{perm_number}/$GLOBALS{SM}{default_perm_per_chunk});
			## $csize should be 1 if $DEFAULT_PERM_PER_CHUNK > perm_number
			my $outcome = compute_perms(\@plist,$dirs->{TPERMS},
					$GLOBALS{SM}{default_perm_per_chunk},$args->{full_sigma},$dirs->{LOG},$csize);
			die "Unable to compute perms" unless $outcome;
		} 
		my $outfiles = collapse_perms($dirs->{TPERMS},$dirs->{PERMS},$dirs->{LOG});
		open(SIG,">$permsdotfile") || die "Cannot open perms dot file: $permsdotfile\n";
		foreach my $f (@$outfiles){
			print SIG "$f\n";
		}
		close(SIG);
			if(!-z $permsdotfile){
			return 1;
		}else{
			return 0;
		}
}

sub run_wstar{
	my ($analysis,$dirs,$snp_file,$force,$args)=@_;
	my %wparams;
	my $step = Sandman::Step->new(
		logdir=>$dirs->{LOG}."/compute_wstar",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	my $miss_file = $dirs->{PERMS}."/.missing";
	my $skip_count=0;
	foreach my $al(keys %$analysis){
		$wparams{$al}{test_set}=$analysis->{$al}->[0];
		$wparams{$al}{ref_set}=$analysis->{$al}->[1];
		my $wstardir = $dirs->{WSTAR}."/".$wparams{$al}{test_set}.'_VS_'.$wparams{$al}{ref_set};
		$wparams{$al}{wstardir} = $wstardir;
		$wparams{$al}{missingfile}=$miss_file;
		#die($wstardir);
		my $run_wstar = 1;
		if((-d $wstardir || -e basename($wstardir).".Wstar.RData") && !$force){
			if(!yesno("wstar: $wstardir dir found do you wish to recompute?")){
				$skip_count++;
				next;
			}
		}
		if($run_wstar){
				## clear up
			if(-d $wstardir){
				## clear up first
				remove_tree($wstardir,{keep_root => 1,result=> \my $del_dirs});
				debug(join("\n",@$del_dirs),$DEBUG);
			}
	
		  my $outcome = compute_wstar(
		  	$wparams{$al}{test_set},$wparams{$al}{ref_set},
		  	$args->{filter},$snp_file,$dirs->{PERMS},
		  	$dirs->{WSTAR},$miss_file,$step,$dirs->{LOG});
		  die ("Unable to compute Wstar for ".$wparams{$al}{test_set}." VS "
		  	.$wparams{$al}{ref_set}) unless $outcome;
		}
	}
	## for the case where the user skips all analysis
	## else we end up dieing as do_step has nothing to do !
	return \%wparams if scalar(keys %wparams) == $skip_count;
	if(do_step($step)){
		return \%wparams;
	}else{
		die "Unable to compute wstar ";
	}
}

sub run_consolidate_wstar_and_w{
	my ($wparams,$dirs,$snp_file,$args)=@_;
	my ($wparams,$dirs) = @_;
	my $step = Sandman::Step->new(
		logdir=>$dirs->{LOG}."/consolidate_wstar_and_w",
		driver=>$DRIVER,
		env_settings=>"R_LIBS=$RLIBS"
	);
	foreach my $al(keys %$wparams){
		my $d = $wparams->{$al}->{wstardir};
		next unless $d; ## this means wstar was not computed on this run so below does not change
		my $wstarfile = $dirs->{WSTAR}.'/'.basename($d).".Wstar.RData";
		$wparams->{$al}->{wstarfile} = $wstarfile;
		debug("Consolidating $d",$DEBUG);
		my $outcome = consolidate_wstar($d,$wstarfile,$step,$dirs->{LOG});
		die("Could not consolidate wstar for $d") unless $outcome;
		my $wfile = $dirs->{WSTAR}.'/'.basename($d).".W.RData";
		debug("Computing W",$DEBUG);
		my $outcome = compute_w($snp_file,$wfile,$wparams->{$al}->{test_set},
			$wparams->{$al}->{ref_set},$args->{filter},
			$wparams->{$al}->{missingfile},$step,$dirs->{LOG});
		 die ("Unable to compute W for ".$wparams->{$al}->{test_set}." VS "
		  	.$wparams->{$al}->{ref_set}) unless $outcome;
		$wparams->{$al}->{wfile} = $wfile;
	}
	if(do_step($step)){
		return $wparams;
	}else{
		die "Unable to consolidate wstar or compute w ";
	}
}

