package Sandman::Step::Job;

use strict;

use Sandman qw/debug get_tmp_file mkpath_wrapper/;

sub new{
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {@_};
	$self =  bless($self,$class);
	return $self;
}

sub command {
	my $self = shift;
	if(@_)	{$self->{command} = shift}
	return $self->{command}
}

sub env_settings {
	my $self = shift;
	if(@_)	{$self->{env_settings} = shift}
	return $self->{env_settings}
}           

sub logdir {
	my $self = shift;
	if(@_)	{$self->{logdir} = shift}
	return $self->{logdir}
}           


sub driver {
	my $self = shift;
	if(@_)	{$self->{driver} = shift}
	return $self->{driver}
}


## getters only - user shouldn't touch these

sub id {
	my $self = shift;
	return $self->{id}
}

sub logfile {
	my $self = shift;
	return $self->{logfile}
}


## private methods 
## lazy get setter only used by this class

sub _set_id {
	my $self = shift;
	if(@_)	{$self->{id} = shift}
	return $self->{id}
}

sub _set_logfile {
	my $self = shift;
	if(@_)	{$self->{logfile} = shift}
	return $self->{logfile}
}
		
	
## interface

## dispatch
## dispatch a job to queue
## sets job id, if not then set id to -1
## calling software can then report this by
## polling a set of jobs

sub dispatch{
	my $self = shift;
	if($self->driver->can('submit')){
			my ($logfile,$jid) =  $self->driver->submit(
				$self->command,
				$self->logdir,
				$self->env_settings);
			$self->_set_id($jid);
			$self->_set_logfile($logfile);
	}else{
			debug("submit must be implemented by ".ref($self->driver));
			return 0;
	}
	return 1 if $self->id() > 0;
	return 0;
}

sub finished{
	my $self=shift;
	if(!$self->id){
		debug("No job id set. Are you sure that this job has been dispatched?");
		return 1;
	}
	if($self->driver->can('check_finished')){
		return $self->driver->check_finished([$self->id]);
	}else{
		debug("check_finished must be implemented by ".ref($self->driver));
		return 0;
	}
}

1;


