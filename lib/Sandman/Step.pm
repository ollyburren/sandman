package Sandman::Step;

use Sandman qw/debug/;
## step is a collection of jobs

use strict;

sub new{
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {@_};
	$self =  bless($self,$class);
	return $self;
}

sub poll_time {
	my $self = shift;
	if(@_)	{$self->{poll_time} = shift}
	return $self->{poll_time}
}

sub env_settings {
	my $self = shift;
	if(@_)	{$self->{env_settings} = shift}
	return $self->{env_settings}
}

sub driver {
	my $self = shift;
	if(@_)	{$self->{driver} = shift}
	return $self->{driver}
}       

sub logdir {
	my $self = shift;
	if(@_)	{$self->{logdir} = shift}
	return $self->{logdir}
}

sub submitted_jobs {
	my $self = shift;
	return $self->{submitted_jobs}
}

sub failed_jobs {
	my $self = shift;
	return $self->{failed_jobs}
}

sub add_job{
	my $self = shift;
	my $job = shift;
	unless($job->isa('Sandman::Step::Job')){
			debug "Job must be an object of type Sandman::Step::Job";
			return 0;
	}
	## we set job all job atts except command at this level
	## in this way it's the same for 
	## all jobs in a step
	$job->driver($self->driver);
	$job->logdir($self->logdir);
	$job->env_settings($self->env_settings);
	if(!$self->{_outstanding_jobs}){
		$self->{_outstanding_jobs}=[$job];
	}else{
		push @{$self->{_outstanding_jobs}},$job;
	}
	return 1;
}

sub execute{
	my $self = shift;
	my (%submitted_jobs,%failed_jobs);
	while(my $j = shift @{$self->{_outstanding_jobs}}){
	#foreach my $j(@{$self->_outstanding_jobs}){
		if($j->dispatch()){
			$submitted_jobs{$j->id}=$j
		}else{
			$failed_jobs{$j->id}=$j;
		}
	}
	$self->{submitted_jobs}=\%submitted_jobs;
	$self->{failed_jobs}=\%failed_jobs;
	return 1 if %submitted_jobs;
}

sub finished{
	my $self=shift;
	my @jids = keys %{$self->submitted_jobs};
	return 1 if !@jids;
	if($self->driver->can('check_finished')){
		return $self->driver->check_finished(\@jids);
	}else{
		debug("check_finished must be implemented by ".ref($self->driver));
		return 0;
	}
}

sub wait_on_complete{
	my $self = shift;
	do{sleep($self->poll_time || 5)} until $self->finished;
}

1;
	
	
