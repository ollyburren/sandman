package Sandman::GRIDDriver::SGE;

use strict;
use base ('Sandman::GRIDDriver');
use Sandman qw/:DEFAULT/;

## must implement submit and check_finished methods 
## see parent class Sandmand::GRIDDriver for details

sub submit{
	my $self = shift;
	my ($cmd,$ldir,$env_settings)=@_;
	if(! $cmd){
		debug("No command given for submission !");
		return (undef,-1);
	}
	if(! -d $ldir){
		if(!mkpath_wrapper($ldir)){
			debug("Could not create log dir $ldir");
			return (undef,-1);
		}
	}
	my ($fh,$tfile) = get_tmp_file($self->cfg->{qscratch});
	my $jname = basename($tfile);
	my $logfile = "$ldir/$jname.log";
	my $qcmd = $self->cfg->{qsub_bin}.' -q '.$self->cfg->{qname};
	if($env_settings){
		$qcmd.= " -v $env_settings ";
	}
	$qcmd.= " -o $logfile -j y $tfile ";
	my $shell = $self->cfg->{shell_path} || $ENV{SHELL};
	print $fh "#!$shell\n$cmd";
	close($fh);
	debug("$qcmd");
	my $out = `$qcmd`;
	my $retval=0;
	if($out =~/Your job ([0-9]+) \("$jname"\) has been submitted/){
		$retval = $1;
	}else{
		debug("Job not submitted !");
		$retval = -1;
		#die($tfile."\t:  ".$out."\n");
	}
	unlink("$tfile") unless $self->cfg->{keep_tmp_files};
	return ($logfile,$retval);
}

sub check_finished{
	my $self = shift;
	my ($joblist,$qstat_params)=@_;
	my $jobs = $self->_get_current_joblist($qstat_params);
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

sub _get_current_joblist{
	my $self = shift;
	my $optional_params = shift;
	my $cmd = $self->cfg->{qstat_bin}." $optional_params |";
	open(QJF,$cmd) || die "$@\n";
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

1;
	
