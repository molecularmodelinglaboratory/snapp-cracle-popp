#!/usr/bin/perl -w

package AIP::Core;
require Exporter;

our @ISA      = qw( Exporter );
our @EXPORT   = qw(  ); # insert sub names here
use vars qw( );
our $version  = 1.00;

#####
use strict;
use Carp;

#####
use AIP::Seq;

##### Ensure DESTROY is called
use sigtrap qw(handler the_rest_is_silence INT QUIT);
sub the_rest_is_silence { die(); } # Calling this on a kill signal ensures DESTROY is called

####
## CONSTRUCTOR AND DESTRUCTOR

# NEW
sub new {
	my $class = shift;
	my $self = {
		'HEADING' => '[A.I.P.] v0.4.0.0 02.09',
		'ABOUT' => {},
		'FILES' => {
			# Logs
			'LOG'		=> 'log.aip',
			'XML' => '',
			
			# Homology
			'H_HASH' => 'h00_hash.out',
			'H_LODS' => 'h01_lods.out',
			'H_MATRIX' => 'h02_matrix.out',
			'H_PCORR' => 'h03_pcorr1.out',
			#'H_IDENT' => 'h04_ident',
			#'H_PCORR2' => 'h05_pcorr',
			
			# Function
			F_HASH => 'f_hash.out',
			F_MATRIX => 'f_matrix.out',
			F_LODS => 'f_lods.out',
			F_PCORR1 => 'f_pcorr1.out',
			F_IDENT => 'f_ident.out',
			F_PCORR2 => 'f_pcorr2.out' },
		
		FLAGS => {
			DEBUG 		=> 0,
			RECOVERY	=> 0,
			COUNT		=> 0 },
		PARAM => {
			DB => 'nr',
			NUMITERATIONS => 500,
			SAVEALL => 1,
			HYPOTHETICAL => 0,
			GRAPH => 0,
			THRESHOLD => 0.8,
			PROGRAM => 'blastp',
			EVALUE => 1e-4,
			IS_PROTEIN => 0,
			IS_STRUCTURE => 0,
			SCORE_TYPE => 'bitscore',
			PROCESSORS => 1,
            REFSEQ_ONLY => 0 },
		DATE => {
			STAMP	=> undef,
			YEAR	=> undef,
			MONTH	=> undef,
			DAY		=> undef,
			HOUR	=> undef,
			MINUTE	=> undef,
			SECOND	=> undef },
		DIR => {
			MAIN	=> undef,
			BLAST 	=> 'BLAST/',
			DATA	=> undef,
			GRAPHS	=> undef,
			XML		=> undef },
		TARGET => {
			ID	=> undef,
			SEQ	=> undef,
			TYPE => undef },
		HANDLE => {
			LOG => undef,
		}
	};
	
		
	bless($self,$class);
	$self->date();
	return $self;
}

# DESTROY
sub DESTROY {
	my $self = shift;
	
	print "\n";
	print "Destroying Core...\n";
	print "Total run time: ".$self->time_passed()."\n";
	die("done.\n");
}

####
## ERROR & WARNING

# _open_warn
sub _open_log {
	my $self = shift;
	my $dir	= shift || $self->{DIR}{MAIN};
	open(my $fh,'>'.$self->{FILES}{LOG})
		or croak('Error file not opened: '.$@);
	#$self->{HANDLE}{W} = fileno(WARNING);
	return $fh;
}

####
## DATE & TIME

#***************************************
#	Format:		YEAR-MONTH-DAY_HOUR-MINUTE
#	Example:	08-01-01_00-01
#	e.g.:		January 1, 2008, 12:01 a.m.

# DATE
sub date {
	my $self = shift;
	my $temp;
	
	##### Get time
	(	$self->{DATE}{SECOND},
		$self->{DATE}{MINUTE},
		$self->{DATE}{HOUR},
		$self->{DATE}{DAY},
		$self->{DATE}{MONTH},
		$self->{DATE}{YEAR},
		$temp,
		$temp,
		$temp	) = localtime(time);
	undef $temp;
	
	##### Modify time (where needed)
	$self->{DATE}{YEAR} = sprintf("%02d", 
		($self->{DATE}{YEAR} + 1900) % 100);			# getting last two digits of the year
	$self->{DATE}{MONTH} += 1;					# correcting MONTH for proper viewing 
	$self->{DATE}{STAMP} = sprintf("%02d-%02d-%02d_%02d-%02d-%02d",
		$self->{DATE}{YEAR},
		$self->{DATE}{MONTH},
		$self->{DATE}{DAY},
		$self->{DATE}{HOUR},
		$self->{DATE}{MINUTE},
        $self->{DATE}{SECOND});					# generating date stamp

	return $self->{DATE}{STAMP};
}

# TIME PASSED
sub time_passed {
	my $self = shift;
	my @newtime = localtime(time);
	my $temp;

	#print "\t$newtime[3]:$newtime[2]:$newtime[1]:$newtime[0]\n";
	#print "\t$self->{DATE}{DAY}:$self->{DATE}{HOUR}:$self->{DATE}{MINUTE}:$self->{DATE}{SECOND}\n";
	
	$newtime[3] -= $self->{DATE}{DAY};
	$newtime[2] -= $self->{DATE}{HOUR};
	$newtime[1] -= $self->{DATE}{MINUTE};
	$newtime[0] -= $self->{DATE}{SECOND};
	
	if($newtime[3] < 0) {
		$newtime[3] = 60 - abs($newtime[3]); }
	if($newtime[2] < 0) {
		$newtime[2] = 60 - abs($newtime[2]);
		$newtime[3]--; }
	if($newtime[1] < 0) {
		$newtime[1] = 60 - abs($newtime[1]);
		$newtime[2]--; }
	if($newtime[0] < 0) {
		$newtime[0] = 60 - abs($newtime[0]);
		$newtime[1]--; }
	
	return "$newtime[2] hours, $newtime[1] minutes, and $newtime[0] seconds";
}

# GENERATE TIMESTAMP
sub generate_timestamp {
    my @t = localtime(time);
    return sprintf("%02d_%02d_%02d_%02d",$t[3],$t[2],$t[1],$t[0]);
}

####
## GET / SET

sub params {
	my $self = shift;
	my $opt = shift;

	$opt->{id} ? $self->id($opt->{id}) : 0; # save id
	$opt->{db} ? $self->database($opt->{db}) : 0;
	$opt->{evalue} ? $self->evalue($opt->{evalue}) : 0;
	$opt->{refseqonly} ? $self->refseq_only($opt->{refseqonly}) : 0;
    $opt->{numiterations} ? $self->iterations($opt->{numiterations}) : 0;
	$opt->{saveall} ? $self->saveall($opt->{saveall}) : 0;
	$opt->{hypothetical} ? $self->hypotheticals($opt->{hypothetical}) : 0;
	$opt->{graph} ? $self->graph($opt->{graph}) : 0;
	$opt->{threshold} ? $self->threshold($opt->{threshold}) : 0;
	$opt->{processors} ? $self->processors($opt->{processors}) : 0;

	
	return 1;
}

# ID
sub id {
	my $self = shift;
	if (@_) {
		$self->{TARGET}{ID} = shift;
		$self->{PARAM}{IS_STRUCTURE} = isPDB($self->{TARGET}{ID});
		$self->_initiate;
	}
	return $self->{TARGET}{ID};
}

# _INITIATE
sub _initiate {
	my $self = shift;

	$self->{DIR}{MAIN} = './'.$self->{TARGET}{ID}.'/';
	$self->{DIR}{DATA} = $self->{DIR}{MAIN}.'DATA/';
	$self->{DIR}{GRAPH} = $self->{DIR}{MAIN}.'GRAPHS/';
	$self->{DIR}{XML} = $self->{DIR}{MAIN}.'XML/';
	$self->sequence if(!$self->{PARAM}{IS_STRUCTURE});

	return 1;
}

# SEQUENCE
sub sequence {
	my $self = shift;
	
    if (@_) { $self->{TARGET}{SEQ} = shift; }
	elsif (!$self->{TARGET}{SEQ} && $self->{TARGET}{ID} ) {
		eval {
            my $t = acc2seq($self->{TARGET}{ID});
			$self->{TARGET}{SEQ} = $t->[0];
			$self->{ABOUT}{$self->{TARGET}{ID}} = $t->[1];
            $self->check_protein($self->{TARGET}{SEQ});
        };
        if($@ && !$self->{TARGET}{SEQ}) {
            croak('No sequence returned for '.$self->{TARGET}{ID}."\n");
        }
	}
	return $self->{TARGET}{SEQ};
}

# GRAPH
sub graph {
	my $self = shift;
	if (@_) { $self->{PARAM}{GRAPH} = shift; }
	return $self->{PARAM}{GRAPH};
}

# HYPOTHETICALS
sub hypotheticals {
	my $self = shift;
	if (@_) { $self->{PARAM}{HYPOTHETICALS} = shift; }
	return $self->{PARAM}{HYPOTHETICALS};
}

# ITERATIONS
sub refseq_only {
	my $self = shift;
	if (@_) { $self->{PARAM}{REFSEQ_ONLY} = shift; }
	return $self->{PARAM}{REFSEQ_ONLY};
}

# ITERATIONS
sub iterations {
	my $self = shift;
	if (@_) { $self->{PARAM}{NUMITERATIONS} = shift; }
	return $self->{PARAM}{NUMITERATIONS};
}

# SAVE ALL
sub saveall {
	my $self = shift;
	if (@_) { $self->{PARAM}{SAVEALL} = shift; }
	return $self->{PARAM}{SAVEALL};
}

# PROCESSORS
sub processors {
	my $self = shift;
	if (@_) { $self->{PARAM}{PROCESSORS} = shift; }
	return $self->{PARAM}{PROCESSORS};
}

# DB
sub database {
	my $self = shift;
	if (@_) { $self->{PARAM}{DB} = shift; }
	return $self->{PARAM}{DB};
}

# THRESHOLD
sub threshold {
	my $self = shift;
	if (@_) { $self->{PARAM}{THRESHOLD} = shift; }
	if (@_) { $self->{PARAM}{UPPER_THRESHOLD} = shift; }
	return $self->{PARAM}{THRESHOLD};
}

# EVALUE
sub evalue {
	my $self = shift;
	if (@_) { $self->{PARAM}{EVALUE} = shift; }
	return $self->{PARAM}{EVALUE};
}

# PROTEIN
sub check_protein {
	my $self = shift;
	if (@_) {
        my $check = shift;
        $self->{PARAM}{IS_PROTEIN} = isProtein($check);
    }
	return $self->{PARAM}{IS_PROTEIN};
}

# STRUCTURE
sub check_structure {
	my $self = shift;
	if (@_) { $self->{PARAM}{IS_STRUCTURE} = shift; }
	return $self->{PARAM}{IS_STRUCTURE};
}

############################################################
##### EXTRA METHODS ########################################
############################################################

## OUTPUT
sub _output {
	my $self = shift;
	my $file = shift;
	my $output = shift;
	
	eval{
		open(my $fh, '>', $file);
		print $fh $$output;
		close($fh);
	};
	if($@) {
		carp('File error: '.$file.'. '.$@."\n");
	}
	
	return 1;
}

# OUTPUT HASH
sub _output_hash {
	my $self = shift;
	my $file = shift;
	my $hash = shift;
	my $list = @_ ? shift : 0;
	my $meta = @_ ? shift : 0;
	
	my $out = '';
	my $omit = {};
	
	# meta data, i.e. min, max
	if($meta) {
		foreach my $i (@{$meta}) {
			#print $i."\n";
			$out .= '#'.$i.' '.$hash->{$i}."\n";
			$omit->{$i} = 1;
		}
	}
	
	# list, i.e. target first
	if($list) {
		foreach my $i (@{$list}) {
			$out .= '#'.$i.' '.scalar(keys(%{$hash->{$i}}))."\n";
			$omit->{$i} = 1;
			foreach my $j (keys %{$hash->{$i}}) {
				#print $i.' :: '.$j.' :: '.$hash->{$i}{$j}."\n";
				$out .= sprintf('%s: %.5f'."\n",$j,$hash->{$i}{$j});
			}
		}
	}
	
	# the rest of it
	foreach my $i (keys %{$hash}) {
		next if(exists $omit->{$i});
		$out .= '#'.$i.' '.scalar(keys(%{$hash->{$i}}))."\n";
		foreach my $j (keys %{$hash->{$i}}) {
			$out .= sprintf('%s: %.5f'."\n",$j,$hash->{$i}{$j});
		}
	}
	
	chomp($out);
	$self->_output($file,\$out);
    return 0;
}

sub _output_matrix {
    my $self = shift;
    my $file = shift;
    my $grid = shift;
    
    my $output = '';
    
    map {
        map {
            $output .= $_."\t";
        } @{$_};
        chomp($output);
        $output .= "\n";
    } @{$grid};
    chomp($output);
    
    $self->_output($file,\$output);
    return 0;
}

sub _input_matrix {
    my $self = shift;
    my $file = shift;
    my $matrix = [];
    my $count = 0;
    
    open(my $fh,'<',$file) or carp('Bad matrix file!');
    
    while (my $line = <$fh>) {
        chomp($line);
        $matrix->[$count] = [];
        @{$matrix->[$count]} = split(/[,\s]+/,$line);
        $count++;
    }

    return $matrix;    
}

############################################################

=for later

***************************************
This function is not complete

sub debug {
	my $self = shift;
	confess "usage: thing->debug(level)"	unless @_ == 1;
	my $level = shift;
	if (ref($self))  {
		$self->{"_DEBUG"} = $level;         # just myself
	} else {
		$Debugging        = $level;         # whole class
	}
}

=cut


1;
__END__

=head1 NAME

Core.pm

=head1 PACKAGE

AIP

=head1 SYNOPSIS

Methods needed by all extensions of the AIP module

=head1 DESCRIPTION

=head1 METHODS

=cut

=head1 BUGS

=over 4

=item * No known bugs

=back

=head1 SEE ALSO

=head1 AUTHOR

Stephen J. Bush, E<gt>muppetjones at gmail dot com E<lt>

=head1 COPYRIGHT & LICENSE

Copyright (c) 2007 by Stephen Bush

=cut
