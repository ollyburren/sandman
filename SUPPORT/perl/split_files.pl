#!/usr/bin/perl

use strict;
use File::Basename;
use File::Path qw/make_path/;
use File::Find;
use POSIX;
use File::Temp qw/tempfile/;

	

###############
#CONFIGURATION#
###############
my $TEST = 0; #IF set to true allows us to test the script by running only a few jobs
my $ROOT_DIR = '/sandman/1KGenome/split_sigma/';
my $DATADIR = '/sandman/1KGenome/SIGMA_TMP/';
my $FILEPATTERN='RData$';
$FILEPATTERN='chr1_.*\.RData$' if $TEST;
my $QUEUE = 'all.q';
my $QSUBCMD="qsub -q $QUEUE -v R_LIBS=/sandman/R_packages";
my $RSCRIPT='/usr/bin/Rscript /home/sandman/GIT_REPOS/sandman/SUPPORT/R/split_sigma.R';

my $BASE_DIR = $ROOT_DIR;
my $logdir = "${BASE_DIR}log/";
if(! -d $logdir){
	mkdir($logdir);
}else{
	`rm $logdir*.log`;
}

find(sub {
		if(/$FILEPATTERN/){
			my ($fh,$fname) = &get_tmp_file();
			my $logfile = $logdir.basename($fname).".log";
			my $cmd = "$QSUBCMD -o $logfile -j y $fname";
my $s=<<SHELL;
#!/bin/bash
$RSCRIPT in.file=\\'$File::Find::name\\' out.dir=\\'/sandman/1KGenome/tmp/\\'
SHELL
			print $fh $s;
			close($fh);
			debug("$cmd");
			`$cmd`;
			unlink($fname) unless $TEST;
			exit if $TEST;
		}
	}, ($DATADIR));


sub debug{
	my $msg=shift;
	my $ts = POSIX::strftime("[%m/%d/%Y %H:%M:%S]", localtime);
	print join("\t",$ts,$msg)."\n";
}

sub get_tmp_file{
	my $dir = shift;
	my $template = 'split_XXXXX';
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}
