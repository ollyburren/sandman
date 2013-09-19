package Sandman::GRIDDriverFactory;
 
use strict; 

sub instantiate {
	my $invoker	= shift;
	my $requested_type = shift;
	my $location	= "Sandman/GRIDDriver/$requested_type.pm";
	my $class	= "Sandman::GRIDDriver::$requested_type";
	require $location;
	return $class->new(@_);
}

1;                                                         
