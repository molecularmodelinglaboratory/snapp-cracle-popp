#!/usr/bin/perl -w

package AIP::Director_nt;
require Exporter;

our @ISA      = qw( Exporter AIP::Core );
our @EXPORT   = qw(  ); # insert sub names here
use vars qw( );
our $version  = 1.00;

#####
#use diagnostics;
use strict;

use Carp;
use Switch;
use AIP::Core;
use AIP::Initialize_nt;
use AIP::Extract;
use AIP::AXLE;

##### Ensure DESTROY is called
#use sigtrap qw(handler the_rest_is_silence INT QUIT);
#sub the_rest_is_silence { die(); } # Calling this on a kill signal ensures DESTROY is called
sub the_rest_is_silence { exit; }
$SIG{'INT'} = \&the_rest_is_silence;
$SIG{'QUIT'} = $SIG{'INT'};

####
## CONSTRUCTOR & DESTRUCTOR

# NEW
sub new {
	my $class = shift;
	my $self = AIP::Core->new();
        
	bless($self,$class);
	$self->date();
    $self->_open_log();
	$self->_open_log();
	return $self;
}

# DESTROY
sub DESTROY {
	#my $self = shift;
	
	print 'Director done.'."\n";
	#Total run time: '.$self->time_passed()."\n";
	#die('done.'."\n");
	return 1;
}

####
## ERROR AND & WARNING

# _open_log
sub _open_log {
    my $self = shift;
    $self->{HANDLE}{LOG} = $self->SUPER::_open_log($self->{DIR}{MAIN});
    #open(WARNING,'>&'.$self->{HANDLE}{WARNING})
    #    or croak("Warning file not opened: $@");
}


#####
## STATUS & RECOVERY

# CHECK STATUS
sub check_status {
    my $self = shift;
    my $recovery = 0;
    
    if(!-d $self->{DIR}{BLAST}) {
		eval {
				mkdir($self->{DIR}{BLAST},0755)  || print $!;
		};
		if($@) {
			croak('Unable to create directory.'."\n".'$@'."\n");
		}
	}
	else {
		my @args = ( 'find '.$self->{DIR}{BLAST}.'*.blast -size -400c | xargs -r rm');
		eval{
			system(@args) == 0
				or carp "system @args failed: $?";
		};
		if($@) {
			carp('bad BLAST files not deleted. '.$@);
		}
	}
    
    if(-d $self->{DIR}{MAIN}) {
        
        # XML? easier search
        if(-d $self->{DIR}{XML}) {
			#&& -f $self->{DIR}{XML}.$self->{FILES}{XML}) {
            # open XML file
            # read 2nd line (skip first -- project data)
            # Find last step saved
            # return 1 if($self->_recovery($recovery));
        }
		else {
			eval {
				mkdir($self->{DIR}{XML},0755)  || print $!;
			};
			if($@) {
				croak('Unable to create directory.'."\n".'$@'."\n");
			}
		}
        
        # Directory? for pre-XML versions
        if(-d $self->{DIR}{DATA}) {
            
            # HOMOLOGY
            if(-f $self->{DIR}{DATA}.$self->{FILES}{H_HASH}) {
                $recovery++;
                if(-f $self->{DIR}{DATA}.$self->{FILES}{H_MATRIX}) {
                    $recovery++;
                    if(-f $self->{DIR}{DATA}.$self->{FILES}{H_LODS}) {
                        $recovery++;
                        if(-f $self->{DIR}{DATA}.$self->{FILES}{H_PCORR}) {
                            $recovery++;
                            #if(-f $self->{DIR}{DATA}.$self->{FILES}{H_IDENT}) {
                            #    $recovery++;
                            #    if(-f $self->{DIR}{DATA}.$self->{FILES}{H_PCORR2}) {
                            #        $recovery++;
                            #    }
                            #}
                        }
                    }
                }
            }
            return 1 if($self->_recovery($recovery));
            
            # FUNCTION
            $recovery = 10;
            if(-f $self->{DIR}{DATA}.$self->{FILES}{F_HASH}) {
                $recovery++;
                if(-f $self->{DIR}{DATA}.$self->{FILES}{F_MATRIX}) {
                    $recovery++;
                    if(-f $self->{DIR}{DATA}.$self->{FILES}{F_LODS}) {
                        $recovery++;
                        if(-f $self->{DIR}{DATA}.$self->{FILES}{F_PCORR}) {
                            $recovery++;
                            #if(-f $self->{DIR}{DATA}.$self->{FILES}{F_IDENT}) {
                            #    $recovery++;
                            #    if(-f $self->{DIR}{DATA}.$self->{FILES}{F_PCORR2}) {
                            #        $recovery++;
                            #    }
                            #}
                        }
                    }
                }
            }
            return 1 if($self->_recovery($recovery));
            
        }
		else {
			eval {
				mkdir($self->{DIR}{DATA},0755)  || print $!;
			};
			if($@) {
				croak('Unable to create directory.'."\n".'$@'."\n");
			}
		}
        
        # Graphs? for graph only runs
        #if(-d $self->{DIR}{GRAPHS}) {
        #    {$self->{HANDLE}{E}} 
        #}
		if(!-d $self->{DIR}{GRAPH}) {
			eval {
				mkdir($self->{DIR}{GRAPH},0755)  || print $!;
			};
			if($@) {
				croak('Unable to create directory.'."\n".'$@'."\n");
			}
		}
    }
	else {
		eval {
			mkdir($self->{DIR}{MAIN},0755);
			mkdir($self->{DIR}{DATA},0755);
			mkdir($self->{DIR}{XML},0755);
			mkdir($self->{DIR}{GRAPH},0755);
		};
		if($@) {
			croak('Unable to create directories.'."\n".'$@'."\n");
		}
	}
   
    return 0;
}

# _RECOVERY
sub _recovery {
    my $self = shift;
    my $recovery = shift;
    my $choice;
    my $blurb = 'AIP has found an earlier session for your target. Would you like to ';
    my $msg = '';
    
    # handle old sessions
    switch($recovery) {
        case [0]     { $blurb = 'No earlier homology searches found'; }
        case [10]    { $blurb = 'No earlier function searches found'; }
        case [1..6]  { $msg = 'finish your homology search'; }
        case [7..9]  { $msg = '...how did you get here'; }
        case [11..16]{ $msg = 'finish your function search'; }
        case [17..20]{ $msg = '...how did you get here'; }
    }
    
    # no old sessions, just go back
    if(!$msg) {
        print $blurb."\n";
        return 0;
    }
    
    # warn and check before continuing
    #carp($blurb.$msg.'?'."\n");
    #do{
        #print '(y/n): ';
        #chomp($choice = <STDIN>);
    #} while(!$choice =~ /[yn]/);
    
    # nope, stop
    #if($choice =~ /n/) {
    #    return 0;
    #}
    
    # go on
    else {
        $self->{FLAGS}{RECOVERY} = $recovery;
    }
    
    return 1;
}

####
# RUN

sub run {
    my $self = shift;
    my $init = AIP::Initialize_nt->new($self);
    my $axle = AIP::AXLE->new($self->{DIR});
    
	### Create the network
	# $init should already have all the variables needed to run
	#	-- inherited from Director
	# Will run consecutive, iterative BLAST queries to formalize
	# the network.
	# Will output the resultant hash in the DATA directory
	# Should retain the hash -- no need to save
	#	-- created in $init->new()
	print $self->generate_timestamp.'> Initializing network...'."\n";
	$init->initialize_network;
	$axle->log_init($self->{TARGET},$self->{PARAM},$self->{DATE});
	$axle->commit;
    
	### Normalize the network
	# $init should have everything it needs to run
	# Will iterate through the hash, collecting sequence average
	#	-- will  be biased based upon number of times the sequence shows up
	#	-- i.e. resolution error -- fuzzy near the edges
	# Normalizes from the hash -- not the matrix!
	# Should retain the hash, -- no need to save
	print $self->generate_timestamp.'> Normalizing the network...'."\n";
	eval {
		$init->normalize_network();
	};
	if($@) {
		carp($@."\n");
	}
	
	### Create the Matrix and Correlate
	# Again, should have access to everything
	# This will output the hash in full matrix format and
	# call the R package from Perl (using the module, not a cmd
	# line call) to perform the correlation and output the
	# correlated matrix.
	# Note: Initialize only correlates the FIRST time
	# Extract must perform the second correlation after
	# creating the identity table from the threshold.
	print $self->generate_timestamp.'> Correlating the data...'."\n";
	$init->create_matrix();
	my $original = $init->correlate();
	
	### Save the Parameters
	# Send parameters to AXle-P
	

	### Run Extract Threads
	# Extract at thresholds 95..[5]..70
	# Join after 4 threads are started
	# Finish last two
    print $self->generate_timestamp.'> Isolating by threshold and extracting data...'."\n";
	$self->{subnet} = {};
    $self->{thresh} = [
        [0.95,1.00],
        [0.90,1.00],
		[0.85,1.00],
        [0.80,1.00],
        [0.70,1.00]
       # [0.80,0.90],
       # [0.70,0.80]
    ];

	$self->threshold_thread($init,$original);

   
#    my $t = {};
#	my $m = {};
#	my $thread_count = 0;
#    my $ct = {
#		't95_100' => [0.95,1.00],
#        't90_100' => [0.90,1.00],
#        't80_100' => [0.80,1.00],
#        't70_100' => [0.70,1.00],
#        't80_90' => [0.80,0.90],
#        't70_80' => [0.70,0.80]
#        };
    #$threshold->[$thread_count] = $ct;
#    foreach my $label (keys %{$ct}) {
#		my $lower = $ct->{$label}[0];
#		my $upper = $ct->{$label}[1];
#		$t->{$label} = threads->create(\&$self->threshold_thread($init,$lower,$upper));
#        $thread_count++;
#			
#        if($thread_count >= $self->{PARAM}{PROCESSORS}) {
#            foreach my $i (keys %{$t}) {
#				$m->{$i} = $t->{$i}->join;
#            }
#            undef $t;
#            $thread_count = 0;
#	    }
#    }
    
	#map {
	#	my $g = $_;
	#	map {
	#		print $g.'::'.$_.'::'.$self->{subnet}{$g}{$_}."\n";
	#	} keys %{$self->{subnet}{$g}};
	#} keys %{$self->{subnet}};
	
	$axle->log_subnets($self->{subnet},$self->{TARGET}{ID});
   	$axle->log_key($self->{ABOUT});
	$axle->commit;
    $axle->finalize;
	
	if($self->{PARAM}{GRAPH}) {
		print $self->generate_timestamp.'> Graphing [background]...'."\n";
		while(my $gcmd = shift(@{$self->{to_graph}})) {
			system(@{$gcmd});
		}
	}
	
    return 1;
}

# EXTRACT_THREADS
sub threshold_thread {
	my $self = shift;
	my $init = shift;
	my $original = shift;
	#my $lower = shift;
	#my $upper = shift;


    # Create extraction method
    my $extract = AIP::Extract->new($self);
	$extract->set_original($original);
    
    while(my $range = shift(@{$self->{thresh}})) {
		
		$init->matrix_reset($original);
		
		my ($lower,$upper) = @{$range};
		my $ll = $lower * 100;
		my $ul = $upper * 100;
		my $label = 't'.$ll.'-'.$ul;
        $init->update_filenames_with_threshold($lower,$upper);
        
        # EMERGENCY: we may have variable conflicts!!
        # doesn't look like it, but...
        # Note: use 'threads', not 'Thread' -- this should
        # eliminate any potential conflicts...
        
		# Note: we have an error here
		# Likely our $init stopped working
		
        ### Identity table
		print $self->generate_timestamp.'> Applying threshold identity and correlation '.$ll.' and '.$ul.'...'."\n";
        next if(!$init->find_identity($lower,$upper));
        
        ### Correlation
        my $grid = $init->correlate();
    
        ### Subnet Identification
        # Potential problem: do we need the second threshold for this?
        # Ideally, no. These cliques should have been isolated from the
        # threshold application of the identity table.
        # We are still in the thread b/c every threshold still needs
        # each step done individually.
        print $self->generate_timestamp.'> Extracting...'."\n";
		my $do_network = 0;
		if($upper < 1) { # correct for the slice
			$lower = 0.85;
			$upper = 1.0;
		}
		$extract->matrix_to_hash($grid,$lower);
        $self->{subnet}{$label} = $extract->find_subnet($self->{TARGET}{ID});
        if($self->{PARAM}{GRAPH}) {
#			print 'adding graph '.$label."\n";
            my $gcmd = $extract->graphviz_subnet(
				$self->{subnet}{$label},
				$label,
				$self->{TARGET}{ID});
            push(@{$self->{to_graph}},$gcmd);
			
			if($lower == 0.85) {
				my $gcmd = $extract->graphviz_network($lower);
				push(@{$self->{to_graph}},$gcmd);
			}
        }
        
    }
    
	if($self->{PARAM}{GRAPH}) {
		my $gcmd = $extract->graphviz_network(0.5,1);
            push(@{$self->{to_graph}},$gcmd);
	}
	
    return 1;
}

#***************************
#** END FILE ***************
#***************************

1;

__END__


=head1 NAME

Director.pm

=head1 PACKAGE

AIP

=head1 AUTHOR

Stephen Bush, E<gt> sjbush@unc.edu E<lt>

=head1 COPYRIGHT & LICENSE

Copyright (c) 2007 by Stephen Bush

=head1 SYNOPSIS

Allows complete recovery from any step of a started AIP session.
Typically useful for finishing Collector data before beginning
extraction.

=head1 DESCRIPTION

=head1 METHODS

=head1 BUGS

=over 4

=head3 * No known bugs

=back

=head1 SEE ALSO


=cut
