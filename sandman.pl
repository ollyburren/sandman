#!/usr/bin/perl

use strict;

=head1 NAME
	sandman

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHOR

=head1 BUGS

=cut

## LIBS
use File::Path qw/make_path remove_tree/;
use File::Basename;
use File::Find;
use File::Temp qw/tempfile/;
use POSIX qw/strftime ceil/;
use Term::ReadKey;
use Data::Dumper;
use Cwd;
use Config::IniFiles;
#use Storable qw( store retrieve );

## GLOBALS

use vars(qw/
	%DIRS $GRID_SCRATCH $SNP_CATALOGUE $MART_HOST
	%RSCRIPTS %MAND_PARAM %DEFAULT_PARAM $SANDMAN_ROOT
	/);

%RSCRIPTS=(
	create_support_files=>'create_support_files.R',
	compute_sigma=>'compute_sigma.R',
	compute_mvs_perms=>'compute_mvs_perms.R',
  misc_functions=>'miscFunctions.R' ## these need a tidy
);

%MAND_PARAM=(
	file=>{
		GLOBAL=>['geneset_file'],
		DATASET=>['pval_file']
	},
	dir=>{
		GLOBAL=>['base_dir'],
		DATASET=>['gt_dir'] ## currently MANDATORY as otherwise cannot compute SIGMA
	},
	array=>{
		GLOBAL=>['analysis']
	});

$SANDMAN_ROOT = $ENV{SANDMAN};
if(! $SANDMAN_ROOT){
	debug("SANDMAN env not set using current working directory");
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
	
	

%DIRS=(
	SCRATCH=>'scratch',
	TSNP=>'scratch/snps',
	SSNP=>'scratch/sigma_snps',
	TPERMS=>'scratch/perms',
	SIGMA=>'sigma', ## change this to use scratch/sigma as more logical
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

my %GLOBALS;
my $sm_cnf = tie %GLOBALS,'Config::IniFiles',(-file=>$bs_sandman_file); 





$SNP_CATALOGUE=$GLOBALS{SM}{default_snp_catalogue}; # probably make this web enabled

## TODO think about rewriting so that biomart is not always need (i.e. keep gene locs in
## bed file ?
$MART_HOST=$GLOBALS{SM}{default_mart_host}; ## can be overidden by analysis cnf



##SETUP R SCRIPT PATH
foreach my $k(keys %RSCRIPTS){
	my $script_path = "$SANDMAN_ROOT/R/$RSCRIPTS{$k}";
	if(-e $script_path){
		$RSCRIPTS{$k}= $script_path;
	}else{
		debug("Cannot find $script_path check SANDMAN installation");                 
		exit(1);
	}
}


## ARGUMENT PROCESSING
## replace with formal get options logic
my $inifile = shift;

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
	#die(Dumper($cfg->GroupMembers('DATASET')));
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

my %ds_params;
my $cfg=parse_args($inifile);	 
#if(!-e './params.store'){
	
	foreach my $ds ($cfg->GroupMembers('DATASET')){
		print $ds."\n";
		$ds_params{$ds}=analyse_dataset($cfg,$ds);
	}
	#store(\%ds_params, './params.store');
#}else{
#	my $hf = retrieve('./params.store');
#	%ds_params = %$hf;
#}

## for testing 



## final computation of Z scores combines results from both datasets and is separate

my %zparams;
foreach my $ds (keys %ds_params){
	my %a = %{$ds_params{$ds}};
	foreach my $al(keys %a){
		#debug("Computing Z score for $ds:$a");
		my %p = %{$a{$al}};
		next unless $p{wfile}; ## need this or else cannot compute !
		push @{$zparams{$al}{wstarfile}},$p{wstarfile};
		push @{$zparams{$al}{wfile}},$p{wfile};
	}
}

#die(Dumper(%zparams));

my $results_dir = $cfg->val('GLOBAL','base_dir');

foreach my $analysis(keys %zparams){
		my $resoutfile = $cfg->val('GLOBAL','base_dir')."/$analysis.Z.RData";
		compute_Z($zparams{$analysis}{wfile},
			$zparams{$analysis}{wstarfile},
			$resoutfile,$results_dir."/$DIRS{LOG}"
			);
}



exit(1);


sub analyse_dataset{
	my ($cfg,$dataset)=@_;
	if(my $snp_cat = $cfg->val($dataset,'snp_catalog')){	
  	  if(-e $snp_cat){
	    $SNP_CATALOGUE = $snp_cat
	  }else{
	    debug("Cannot find $snp_cat for this datasource: using $SNP_CATALOGUE");
	  }
	}
	my %FORCE_STEP=(
		support=>0,
		sigma=>0);
		
	(my $DS_NAME = $dataset) =~ s/DATASET\s+//;
	
	##for ease !!
	my %targs = (
	base_dir=> $cfg->val('GLOBAL','base_dir').'/',
	pval_file => $cfg->val($dataset,'pval_file'),
	datasource_name => $DS_NAME,
	geneset_file => $cfg->val('GLOBAL','geneset_file'),     
	exclude_region_file => $cfg->val('GLOBAL','exclude_region_file'),
	tss_extension => $cfg->val('GLOBAL','tss_extension'),
	gt_dir=> $cfg->val($dataset,'gt_dir'),
	#test_set=>$cfg->val('GLOBAL','test_set'),
	#ref_set=>$cfg->val('GLOBAL','ref_set'),
	filter=>$cfg->val('GLOBAL','filter'),
	perm_number=>$cfg->val('GLOBAL','perm_number'),
	);
	
	## set up analysis
	my %analysis;
	foreach my $al($cfg->val('GLOBAL','analysis')){
		my ($tname,$test,$ref)=split(/\s+/,$al);
		$analysis{$tname}=[$test,$ref];
	}
	## STEP 1: SETUP
	
	
	
	my $dirs = make_support_dir($targs{base_dir}.$targs{datasource_name},0);
	
	##next make snp support files
	
	my $snp_file = $dirs->{'TSNP'}.'/snp.gr.RData';
	my $run_support = 1;
	if(-e $snp_file && !$FORCE_STEP{support}){
		$run_support = 0 unless yesno("support files found do you wish to overwrite?");
	}
	
	if($run_support){ 
		debug("Generating support files $snp_file");
		my @jids = make_support_files(
			$targs{'pval_file'},
			$targs{'geneset_file'},
			$snp_file,
			$targs{'exclude_region_file'},
			$GLOBALS{SM}{tabix_bin},$SNP_CATALOGUE,$GLOBALS{SM}{tabix_chunksize},
			$MART_HOST,$targs{tss_extension},
			$dirs->{LOG});
		do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@jids);
		## automatically recompute sigmas
		$FORCE_STEP{sigma}=1;
	}
	## STEP 2: COMPUTE SIGMA'S
	
	
	my $rhash;
	my $sigdotfile = $dirs->{SIGMA}."/.sigma";
	my $run_sigma = 1;
	if(-e $sigdotfile && !$FORCE_STEP{sigma}){
		$run_sigma = 0 unless yesno("sigma files found do you wish to recompute?");
	}
	
	if($run_sigma){
		$rhash = create_sigmas(
		$snp_file,
		$dirs->{SIGMA},
		$targs{gt_dir},
		$dirs->{SSNP},
		$dirs->{LOG});
		my @jids = keys(%$rhash);
		return if !@jids;
		do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@jids);
		#write the snpfile to sigma mapping to a file
		#doubles as evidence that we have done things and also 
		#means we can run subsequent commands w/o having to run this
		open(SIG,">$sigdotfile") || die "Cannot open sigma dot file: $sigdotfile\n";
		foreach my $jid (keys(%$rhash)){
			my %h=%{$rhash->{$jid}};
			print SIG join("\t",$h{snpfile},$h{sigmafile})."\n";
		}
		close(SIG);
		## automatically recompute perms
		$FORCE_STEP{perms}=1;
	}
	
	## STEP 3a: COMPUTE PERMS
	
	## read in sigma .file
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
	
	my $permsdotfile = $dirs->{PERMS}."/.perms";
	my $run_perms = 1;
	if(-e $permsdotfile && !$FORCE_STEP{perms}){
		$run_perms = 0 unless yesno("perm files found do you wish to recompute?");
	}
	if($run_perms){
		## clear up
		if(-d $dirs->{PERMS}){
			## clear up first
			remove_tree(($dirs->{PERMS},$dirs->{TPERMS}),{keep_root => 1,result=> \my $del_dirs});
			debug(join("\n",@$del_dirs));
		}
	
		my $csize = ceil($targs{perm_number}/$GLOBALS{SM}{default_perm_per_chunk});
		## $csize should be 1 if $DEFAULT_PERM_PER_CHUNK > perm_number
		my $rhash = compute_perms(\@plist,$dirs->{TPERMS},
				$GLOBALS{SM}{default_perm_per_chunk},$dirs->{LOG},$csize);
		my @jids = keys(%$rhash);
		return if !@jids;
		do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@jids);
		my $outfiles = collapse_perms($dirs->{TPERMS},$dirs->{PERMS},$dirs->{LOG});
		open(SIG,">$permsdotfile") || die "Cannot open perms dot file: $permsdotfile\n";
		foreach my $f (@$outfiles){
			print SIG "$f\n";
		}
		close(SIG);
		## automatically recompute wstar
		$FORCE_STEP{wstar}=1;
	}
	
	
	## STEP 3b: CHECK PERMS
	
	#check the perms created and create missing snps list if the snpfile and perm rownames lists
	#are not the same
	my $check_permsdotfile = $dirs->{PERMS}."/.checked";
	#die( $check_permsdotfile);
	if(!-e $check_permsdotfile){
		check_perms($permsdotfile,$snp_file,$dirs->{LOG});
		if(-e $dirs->{PERMS}."/.error"){
			die("There has been an error generating permutations and they are malformed\n");
		}elsif(-e $dirs->{PERMS}."/.missing"){
			debug("There are snps missing from perm files!");
			if(!yesno("Do you wish to continue?")){
			exit(1);
			}else{
				open(MISS, $dirs->{PERMS}."/.missing") || die "cannot open ".$dirs->{PERMS}."/.missing\n";
				while(<MISS>){
					chomp;
					debug($_);
				}
			}
		}
	}
	
	
	##STEP 4: COMPUTE WSTARS FOR PERMS
	my @wstar_jids;
	my %wparams;
	foreach my $al(keys %analysis){
		$wparams{$al}{test_set}=$analysis{$al}->[0];
		$wparams{$al}{ref_set}=$analysis{$al}->[1];
		my $wstardir = $dirs->{WSTAR}."/".$wparams{$al}{test_set}.'_VS_'.$wparams{$al}{ref_set};
		#die($wstardir);
		my $run_wstar = 1;
		if((-d $wstardir || -e basename($wstardir).".Wstar.RData") && !$FORCE_STEP{wstar}){
			$run_wstar = 0 unless yesno("wstar: $wstardir dir found do you wish to recompute?");
		}
		if($run_wstar){
				## clear up
			if(-d $wstardir){
				## clear up first
				remove_tree($wstardir,{keep_root => 1,result=> \my $del_dirs});
				debug(join("\n",@$del_dirs));
			}
			$wparams{$al}{wstardir} = $wstardir;
			#next;
			#push @wdirs,$wstardir;
		  my $miss_file = $dirs->{PERMS}."/.missing";
			$wparams{$al}{missingfile}=$miss_file;
		  my $jids = compute_wstar(
		  	$wparams{$al}{test_set},$wparams{$al}{ref_set},
		  	$targs{filter},$snp_file,$dirs->{PERMS},
		  	$dirs->{WSTAR},$miss_file,$dirs->{LOG});
		  push @wstar_jids,@$jids;
		  $wparams{$al}{jids} = $jids;
		}
	}
	do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@wstar_jids);
	foreach my $al(keys %wparams){
		my $d = $wparams{$al}{wstardir};
		next unless $d; ## this means wstar was not computed on this run so below does not change
		my $wstarfile = $dirs->{WSTAR}.'/'.basename($d).".Wstar.RData";
		$wparams{$al}{wstarfile} = $wstarfile;
		debug("Consolidating $d");
		## TODO ALTER SO ALL CONSOLIDATION HITS Q AT THE SAME TIME
		my $jid = consolidate_wstar($d,$wstarfile,$dirs->{LOG});
		do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished([$jid]);
		my $wfile = $dirs->{WSTAR}.'/'.basename($d).".W.RData";
		debug("Computing W");
		my $jid = compute_w($snp_file,$wfile,$wparams{$al}{test_set},
			$wparams{$al}{ref_set},$targs{filter},
			$wparams{$al}{missingfile},$dirs->{LOG});
		do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished([$jid]);
		$wparams{$al}{wfile} = $wfile;
	}
	return \%wparams;	
}

## SUBROUTINES  

sub consolidate_wstar{    
	my ($wdir,$outfile,$log_dir) = @_;
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
  dispatch_Rscript($cmd,"$log_dir/consolidate_wstar");
  #`$cmd`;                                   
  #unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};    
}

sub compute_Z{
	# need Wstar, W (i.e.pvals) and snpnames in and snpnames out for each dataset
	my ($wfiles,$wstarfiles,$outfile,$logdir) = @_;
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
  return(dispatch_Rscript($cmd,"$logdir/compute_Z"));
  #`$cmd`;                                   
  #unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};    
}
	
	


sub make_support_dir{
	my ($dir,$force)=@_;
	if(-d $dir && !$force){
		debug("Directory exists and force flag not set.");
		if(!yesno("Directory $dir already exists are you sure you want to proceed ?")){
			exit(1);
		}
	}elsif(-d $dir && $force){
		debug("Directory exists and force flag set. Deleting directory ..");
		## perhaps check interactively as this could be incredibly dangerous ?
		remove_tree($dir,{keep_root => 1,result=> \my $del_dirs});
		debug(join("\n",@$del_dirs));
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

## FUNCTION CREATES SNP SUPPORT FILE
## TODO MAKE THIS COMPATIBLE WITH RUNNING ON THE Q ?
## AS FOR LOTS OF GENES THIS BIT WILL BE SLOW.

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

sub check_perms{
  my ($dotfile,$snpfile,$log_dir)=@_;
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
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  return(dispatch_Rscript($cmd,"$log_dir/check_perms"));
  #`$cmd`;                                   
  #unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};    
}

       

## helper function that takes a support file and splits it by chromosome
## as then we can compute sigma's on the q much more quickly

sub prepare_snps_for_q{
  my ($outdir,$snpfile,$logdir)=@_; 
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
  #print $R;
  #print TMP $R;
  print $fh $R;
  close($fh);
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  #`$cmd`;
  return(dispatch_Rscript($cmd,"$logdir/prepare_snps_for_q"));
  #unlink($fname);
}

sub create_collapse_perms_script{
	my ($flist,$outfile)=@_;
	my ($fh,$fname) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
  #open(TMP,">$fname") || die "Cannot open tmp file for writing\n";
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
	my($indir,$outdir,$logdir)=@_;
	my %hash;
	find(sub{
			if(/\.perms\.([0-9]+)\.RData$/){
				push @{$hash{$1}},$File::Find::name;
			}
	},$indir);
	debug(Dumper(\%hash));
	my @jids;
	my @tscripts;
	my @outfiles;
	foreach my $p(keys %hash){
		my $flist = join ',', map { qq/'$_'/ }@{$hash{$p}};
		my $outfile = $outdir."/$p.allperms.RData";
		my $fname=create_collapse_perms_script($flist,$outfile);
		push @outfiles,$outfile;
		my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
		my $jid = dispatch_Rscript($cmd,"$logdir/collapse_perms");
		push @jids,$jid;
		push @tscripts,$fname;
	}
	do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@jids);
	unlink @tscripts;
	return \@outfiles;
}
	
	
 
## TODO unlink tmp files 

sub create_sigmas{
  my ($snp_file,$sigma_dir,$gt_dir,$scratch_dir,$log_dir)=@_;
  #create a temp dir
  my @jids;
  my %return;
  if($GLOBALS{SM}{use_q}){
    #we blow away contents of scratch_dir each time
    remove_tree($scratch_dir,{keep_root => 1,result=> \my $del_dirs});
    debug(join("\n",@$del_dirs));
    my $jid = prepare_snps_for_q($scratch_dir,$snp_file,$log_dir);
    do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished([$jid]);
    find(sub{
      debug($_);
      if(/\.RData$/){
      	my $out_file = $sigma_dir.'/'.basename($_,'.RData').".sigma.RData";
      	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_sigma} snp.file=\\'$File::Find::name\\' ";
      	$cmd.= "gt.dir=\\'$gt_dir\\' out.file=\\'$out_file\\' misc.functions=\\'$RSCRIPTS{misc_functions}\\'";
      	my $jid = dispatch_Rscript($cmd,"$log_dir/compute_sigma");
      	debug("Running job $jid");
      	## create a hash of script params that we want to pass downstream
      	$return{$jid}={
      		snpfile=>$File::Find::name,
      		sigmafile=>$out_file
      	};
      }
    },$scratch_dir);
  }else{
    my $out_file = $sigma_dir.'/all.sigma.RData';
    my $cmd = "$GLOBALS{SM}{rscript} compute_sigma.R snp.file=\\'$snp_file\\' ";
    $cmd.= "gt.dir=\\'$gt_dir\\' out.file=\\'$out_file\\'";
    dispatch_Rscript($cmd,"$log_dir/compute_sigma");
    $return{NOQ}={
    		snpfile=>$snp_file,
      	sigmafile=>$out_file
    };
  }
  return \%return;
}

## compute W

sub compute_w{
	my ($snpfile,$outfile,$test_set,$ctrl_set,$filter,$missfile,$logdir)=@_;
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
  #my $cmd = "$RSCRIPT --vanilla $fname > /dev/null 2>&1";
  my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
  #`$cmd`; 
  return(dispatch_Rscript($cmd,"$logdir/compute_w"));
  #unlink($fname) unless $GLOBALS{SM}{keep_tmp_files};  
}  

## for a given snp file and matching sigma file
## create a matrix of perms that we can use to
## compute Wstar
	
sub compute_perms{
	my ($plist,$permdir,$number_perms,$log_dir,$chunksize)=@_;
	$chunksize=1 unless $chunksize; # in case it's not set;
	my %return;
	my @jids;
  if($GLOBALS{SM}{use_q}){
    foreach my $phash(@$plist){
    	my $snpfile = $phash->{snpfile};
    	my $sigmafile = $phash->{sigmafile};
    	for(my $x=1;$x<=$chunksize;$x++){
    		my $out_file = $permdir.'/'.basename($snpfile,'.RData').".perms.$x.RData";
      	#my $cmd = "$RSCRIPT ${RSCRIPT_DIR}compute_mvs_perms.R snp.file=\\'$snpfile\\' ";
  	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_mvs_perms} snp.file=\\'$snpfile\\' ";
      	$cmd.= "sigma.file=\\'$sigmafile\\' out.file=\\'$out_file\\' n.perms=$number_perms misc.functions=\\'$RSCRIPTS{misc_functions}\\'";
      	debug($cmd);
      	my $jid = dispatch_Rscript($cmd,"$log_dir/compute_perms");
      	debug("Running job $jid");
      	## create a hash of script params that we want to pass downstream
      	$return{$jid}={
      		permfile=>$out_file
      	};
      }
    }
  }else{
  	## it might be that we just run the code above and let dispatcher handle where it goes
  	## in this way we just run a number of permutations 
  	$number_perms = $number_perms * $chunksize;
  	debug("Running $number_perms");
  	my $snpfile = @$plist[0]->{snpfile};
  	my $sigmafile = @$plist[0]->{sigmafile};
  	my $out_file = $permdir.'/1.perms.RData';
  	#my $cmd = "$RSCRIPT ${RSCRIPT_DIR}compute_mvs_perms.R snp.file=\\'$snpfile\\' ";
  	my $cmd = "$GLOBALS{SM}{rscript} $RSCRIPTS{compute_mvs_perms} snp.file=\\'$snpfile\\' ";
  	$cmd.= "sigma.file=\\'$sigmafile\\' out.file=\\'$out_file\\' n.perms=$number_perms";
  	dispatch_Rscript($cmd,"$log_dir/compute_sigma");
  	$return{NOQ}={
  		permfile=>$out_file
  	};
  }                                      
  return \%return;
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
	my($test_set,$ctrl_set,$filter,$snpfile,$permdir,$wstardir,$missfile,$logdir)=@_;
	my $outdir = "$wstardir/${test_set}_VS_$ctrl_set";
	mkpath_wrapper($outdir);
	my @jids;
	my @outfiles;
	find(sub{
			if(/([0-9]+)\.allperms\.RData$/){
				my $outfile = "$outdir/Wstar.$1.RData" ;       
				my $fname = create_compute_wstar_script($test_set,$ctrl_set,$filter,$snpfile,$File::Find::name,$missfile,$outfile);
				my $cmd = "$GLOBALS{SM}{rscript} --vanilla $fname";
				my $jid = dispatch_Rscript($cmd,"$logdir/compute_wstar");
				push @jids,$jid;
				push @outfiles,$outfile;
			}
	},$permdir);
	#do{sleep($GLOBALS{SM}{q_poll_time})} until check_task_finished(\@jids);
	## are the tmp files cleared up automatically ?
	return \@jids;
	#return \@outfiles;
}

sub dispatch_Rscript{
	my ($cmd,$ldir)=@_;
	if(! -d $ldir){
		if(!mkpath_wrapper($ldir)){
			debug("Could not create log dir $ldir");
			return 0;
		}
	}
	my ($fh,$tfile) = &get_tmp_file($GLOBALS{SM}{grid_scratch});
	my $jname = basename($tfile);
	my $logfile = "$ldir/$jname.log";
	if($GLOBALS{SM}{use_q}){
		#first get a tmp file - we will write sh to this
		#use the 
		my $qcmd = "$GLOBALS{SM}{qsub} -o $logfile -j y $tfile";
		#open(SCRIPT,">$tfile") || die "Cannot open temp shell script $tfile\n";
		print $fh '#!/bin/bash'."\n$cmd";
		close($fh);
		debug("$qcmd");
		my $out = `$qcmd`;
		my $retval=0;
		if($out =~/Your job ([0-9]+) \("$jname"\) has been submitted/){
			$retval = $1;
		}else{
			die($tfile."\t:  ".$out."\n");
		}
		unlink("$tfile") unless $GLOBALS{SM}{keep_tmp_files};
		return $retval;
	}else{
		## here we attempt to run interactively.
		## TODO need to capture output to a log file nicely
		debug("Running interactively: $cmd");
		`$cmd > $logfile 2>&1`
	}
}
		
	

sub debug{
	my $msg=shift;
	my $ts = strftime("[%m/%d/%Y %H:%M:%S]", localtime);
	print join("\t",$ts,$msg)."\n";
}


sub get_tmp_file{
	my $dir = shift;
	my $template = 'sandman_XXXXX';
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}

sub mkpath_wrapper{
	my @dirs=@_;
	make_path(@dirs,{verbose=>1,error => \my $err});
	if(@$err){
		for my $diag (@$err) {
			my ($file, $message) = %$diag;
			if ($file eq '') {
				debug("general error: $message");
			}else {
				debug("problem creating $file: $message");
			}
		}
		return 0;
	}
	return 1;
}

sub get_current_joblist{
	my $optional_params = shift;
	open(QJF,"qstat $optional_params |") || die "$@\n";
	my %jobs;
	while(<QJF>){
		chomp;
		next unless /^[0-9]+/;
		my @vals=split(/\s+/,$_);
		$jobs{$vals[0]}={
				name=>$vals[2],
				status=>$vals[4]
		}
	}
	return \%jobs;
}

sub check_task_finished{
	my ($joblist,$qstat_params)=@_;
	my $jobs = &get_current_joblist();
	my %outstanding_jobs;
	foreach my $j(@$joblist){
		if(my $details = $jobs->{$j}){
			if($details->{status} !~/E/){
				push @{$outstanding_jobs{$details->{status}}},$j;
			}else{
				## if we get any errors going to be difficult to recover.
				die "Job $j has failed - please check logs\n";
			}
		}
	}
	my @report;
	foreach my $k(keys %outstanding_jobs){
		push @report,scalar(@{$outstanding_jobs{$k}}).": $k";
	}
	return 1 if !@report;
	debug(join("\t",@report,"running on the queue"));
	return 0;
}

=head2 yesno

 NAME: yesno
 ARGS: STRING - represents question to ask
 FUNCTION: Writes a question to the terminal and waits for a y or n keystroke
 RETURNS : Boolean true or false

=cut

sub yesno{
	my($question)=shift;
	if(!$question){
		warn "[WARNING] No question: returning false\n";
		return 0;
	}
	ReadMode('cbreak');
	print "$question?\n";
	my $dyes = "Please enter y or n. Use Ctrl-C to quit.\n";
	while(1){
		my $char = ReadKey(0);
		if($char!~/[yn]/i){
			print $dyes;
		}else{
			ReadMode('normal');
			return $char =~/y/i?1:0;
		}
	}
}

