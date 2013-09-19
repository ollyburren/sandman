package Sandman::GRIDDriver;

use strict;
use Sandman qw/:DEFAULT/;
use Config::IniFiles;

## all grid drivers must implement the following routines

## init - this loads in and sets up grid 

sub new{
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {@_};
	bless($self,$class);
	$self->init();
	return $self;
}

sub init {
	my $self = shift;
	if(! $self->inifile && ! -e $self->inifile){
		die "Cannot configure grid please check ini file set/exists(".$self->inifile.")\n";
	}
	my %g;
	tie %g,'Config::IniFiles',(-file=>$self->inifile);
	## we only need the section for our particular driver
	my $d = (split("::",ref($self)))[-1];
	if(my $cfg = $g{$d}){
		$self->cfg($cfg);
	}else{
		die "Cannot find section for ".__PACKAGE__." ([$d]) in inifile: ".$self->inifile."\n";
	}
}

sub inifile {
	my $self = shift;
	if(@_)	{$self->{inifile} = shift}
	return $self->{inifile}
}

sub cfg {
	my $self = shift;
	if(@_)	{$self->{cfg} = shift}
	return $self->{cfg}
}


## submit - this submits a job to the q

=head2 submit

  NAME: submit 
  ARGS:	command to submit (SCALAR), log file dir path (SCALAR), environment settings (SCALAR OPTIONAL)
  FUNCTION: Submits a job to a GRID computing system
  RETURNS: logfile path (SCALAR), jobid
  
  NOTE logfile should return undef and jobid should be set to -1 if submission fails. 
=cut

sub submit{
	die ("submit routine must be implemented by a subclass of ".__PACKAGE__);
	return 0;
}

=head2 check_finished

  NAME: check_finished
  ARGS:	list of job ids (ARRAYREF), params (SCALAR - OPTIONAL)
  FUNCTION: Checks the queue to see if a job or list of jobs has completed
  RETURNS: TRUE if any job is still on queue or running.
  
=cut

sub check_finished{
	die ("status routine must be implemented by a subclass of ".__PACKAGE__);
	return 0;
}

1;


