#!/usr/bin/perl

use strict;
use File::Basename;
use File::Find;
use POSIX;
###############
#CONFIGURATION#
###############

my $outdir = '/stats/oliver/PERMUTATIONS/';
my $logdir = '/stats/oliver/gen_perms/';
my $tmpdir = '/tmp/';
my $datadir = '/stats/oliver/GENOTYPES/';
my $queue = 'all.q';
my $qsubcmd="qsub -q $queue -v R_LIBS=/stats/oliver/R_packages/";
my $RSCRIPT = '/home/oliver/GIT_REPOS/sandman/SUPPORT/R/run_perm_test.R';
my $filePattern='\.RData$';
#my $filePattern='t1dgc-22.RData$';
my $nperms=100;
my $total_runs=100;

sub wanted{
	if(/$filePattern/){
		my $count=$total_runs;
		while($count>0){
			my $filestub = basename($_,'.RData');
			debug($filestub);
			my $cmd = "$qsubcmd -o $logdir$filestub.$count.log -j y $tmpdir$filestub.$count.sh";
my $s=<<SHELL;
#!/bin/bash
$RSCRIPT datafile=\\'$File::Find::name\\' output_dir=\\'$outdir\\' nperms=$nperms runnumber=$count
SHELL
			open(SCRIPT,">$tmpdir$filestub.$count.sh") || die "Cannot open temporary shell script $tmpdir$filestub.$count.sh\n";
			print SCRIPT $s;
			close(SCRIPT);
			debug("$cmd");
			`$cmd`;
			unlink("$File::Find::name.sh");
			$count--;
		}
		
	}
}
find(\&wanted,  ($datadir));


sub debug{
	my $msg=shift;
	my $ts = POSIX::strftime("[%m/%d/%Y %H:%M:%S]", localtime);
	print join("\t",$ts,$msg)."\n";
}

=cut
