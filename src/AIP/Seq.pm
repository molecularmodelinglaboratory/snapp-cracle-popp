#!user/bin/perl -w
# Filename:	Seq.pm
# Author:	Stephen Bush
# Date:		01.28.07 [1.0 creation]
# Ver:		2.0

package AIP::Seq;
require Exporter;

our @ISA      = qw( Exporter );
our @EXPORT   = qw( acc2seq isProtein isPDB batch2seq ); # insert sub names here
use vars qw( );
our $version  = 1.00;

use strict;
use warnings;
use Carp;
use AIP::Core;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::GenBank;									# used to grab a GenBank sequence

#use sigtrap qw(handler the_rest_is_silence INT QUIT);
sub the_rest_is_silence { exit; }
$SIG{'INT'} = \&the_rest_is_silence;
$SIG{'QUIT'} = $SIG{'INT'};


############################################################
############################################################
sub new {
	my $class = shift;
	my $self = shift;
	$self->{HANDLE} = shift;
	
	bless($self,$class);
	return $self;
}

############################################################
sub _initiate_seq {
	my $self = shift;
	open(ERR,">&$self->{HANDLE}{E}")
		or croak("Error file not opened: $@");
	open(WARNING,">&$self->{HANDLE}{W}")
		or croak("Error file not opened: $@");
}

#################################################################
#The seq_file2seq() function opens a file with a DNA or protein
#sequence and returns the sequence as a string
#Input: 1 string sequence
#    1) filename of file containing sequence
#        a) any format accepted by Bio:SeqIO
#Returns: Bio::Seq object
sub isProtein( $ ) {
	my $seq = shift;
	if($seq !~ /[ilvfpywqhed]+/i) {
		return 0;
	} else { return 1; }
}

# PDB: <number><3 num or lett><maybe lowercase lett>
sub isPDB( $ ) {
	my $id = shift;
	if($id =~ /^\d\w{3}[a-z]?$/) {
		return 1; 
	} else { return 0; }
}


#################################################################
sub random_sequence($ $) {
	my $self = shift;
	$self->{SEED}{LENGTH} = shift;
	my $type = shift;
	
	my @blocks;
	$self->{SEED}{SEQ} = '';
	if($type eq 'protein') {
		@blocks = ("A","R","N","D","B","C","E","Q","Z","H","I","L","K","M","F","P","S","T","W","Y","V");
		for(my $i = 0; $i < $self->{SEED}{LENGTH}; $i++) {
			$self->{SEED}{SEQ} .= $blocks[rand(20)]; }
	}elsif($type eq 'nucleotide') {
		@blocks = ("C","G","T","A");
		for(my $i = 0; $i < $self->{SEED}{LENGTH}; $i++) {
			$self->{SEED}{SEQ} .= $blocks[rand(4)]; }
	}

	$self->{SEED}{ACC} = "random$type";
	return $self->{SEED}{SEQ};
}

######################################################
# Get a single sequence
#	-- Send in an accession number
#	-- get a reference to your sequence
######################################################
sub acc2seq( $ ) {
	#my $self = shift;
	
	my $count = 0; my $seq;
	my $acc = shift;								# get accession number
	my $tacc = [ $acc ];
	while($count < 5) {
		my $gb = Bio::DB::GenBank->new(				# create new GenBank object
			-retrievaltype => 'tempfile' ,
			-format => 'Fasta');
		my $seqio;

		## Getting an unknown error here:
		# Use of uninitialized value in concatenation (.) or string at
		# /usr/lib64/perl5/5.8.8/x86_64-linux-thread-multi/Scalar/Util.pm line 30.
		# thrown somewhere in get_Stream. Still returns sequences though...
		eval{
			$seqio = $gb->get_Stream_by_acc($tacc);
		};	# contact NCBI for data
		if($@ || !$seqio) {										# handle NCBI error
			carp('Sequence problem? '.$@."\n");
			sleep 2;
			$count++;
			next;
		}
		my $seqobj = $seqio->next_seq();			# return SeqIO object
		#if($seqobj->alphabet ne 'protein') {
		#	return -1;
		#}
		$seq = [$seqobj->seq(),$seqobj->desc()];						# return the string sequence
		
		##### Error handling for no sequence returned
		if($seq->[0] =~ /(Error)/) {
			sleep 2;								# pause and keep going
			$count++;								
			next;									# reiterate and reevaluate loop
		}
		$count = 100;									
    }

    $count == 100 ? return $seq : return 0;		# if good seq, return, else 0
}

#################################################################
#batch2seq($)
#***************************************
#I: array reference
#O: haip reference
#
#This function may randomly die with a "Floating point exception."
#because of this function:
#	$gb->get_Stream_by_acc($refbatch)
#NCBI does not like to be queried too often. <batch2seq> is more efficient
#than <acc2seq> because we're getting several sequences at once. While this
#dramatically reduced NCBI query time and thus failure rate, it is not
#foolproof.
#
#We use the function like so:
#
#@$array_ref = keys(%{$haip_ref});
#$new_haip_ref = $self->batch2seq($array_ref);
#
#e.g.
#	A => ( B => 1, C => 2, D => 3)
#	~ [B,C,D]
#	~ ( B => ANKDIL..., C => MQUI..., D => MSSR )
sub batch2seq {
	#my $self = shift;															# GET class data
	my $refbatch = shift;
	my $bad_seq = shift;
	my $isprot = @_ ? shift : 1;
	my $batch = {};

	#print "\t\t".'s1'."\n";
	
	my $type = $isprot ? 'p' : 'na';
	$DB::single=2; # insert at line 9!
	my $count = 0; my $seq;														# INITIALIZE count and seq vars
	my $gb = new Bio::DB::GenBank(												# INITIALIZE a GenBank object
		-retrievaltype => 'tempfile' ,
		-format => 'Fasta');
	my $seqio;																	# INITIALIZE a SeqIO object var
	my $seqobj;																	# INITIALIZE a Seq object var
	
	#print "\t\t".'2..'."\n";
	
	while($count < 5) {															# WHILE we're under 5 tries
		eval{
			#print "\t\t".'3..'."\n";
			$DB::single=2; # insert at line 9!
			$seqio = $gb->get_Stream_by_acc($refbatch);
			$DB::single=2; # insert at line 9!
			#print "\t\t".'4..'."\n";
		};
		if($@ || !$seqio) {														# IF we encountered an error
			if($count == 4) {
				carp('Enough tries. Going home without sequences.'."\n");
				return 1;
			}
			#print "\t\t\t".'a..'."\n";
			#if($@ =~ /the rest is silence/) { croak; }
			
			carp('Bad batch: '.$count.'. Trying again...'."\n");
			sleep 2;															#	-- PAUSE
			$count++;															#	-- INCREMENT our number of tries
			next;																#	-- SKIP to next loop increment	
		}
		else {
			#print "\t\t\t".'b..'."\n";
			$count = 6;
		}
	}

	#print "\t\t".'5..'."\n";
	
	my $reflist = join('|',@{$refbatch});
	my $regex = qr/((?>$reflist))/;
	while($seqobj = $seqio->next_seq()) {									# WHILE we still have new sequences
		$DB::single=2; # insert at line 9!
		
		#print "\t\t".'6..'."\n";
		
		if($seqobj->primary_id =~ /$regex/) {
			my $id = $1;
			
			# seq error?
			if($seqobj->seq() =~ /(ERROR)/) {
				
				#print "\t\t\t".'a..'."\n";
				
				$bad_seq->{$id} = 1;
			}
			
			# proper type?
			elsif($seqobj->alphabet =~ /$type/) {
				
				#print "\t\t\t".'b..'."\n";
				
				$batch->{$id} = [$seqobj->seq(),$seqobj->desc()];
			}
			
			# just plain bad?
			else {
				
				#print "\t\t\t".'c..'."\n";
				
				$bad_seq->{$id} = 1;
			}
		}
		else {}
		
		#print "\t\t".'7..'."\n";
		
		#if($seqobj->seq() =~ /(ERROR)/) {							#*IF it's not a qood sequence
		#	$bad_seq->{$refbatch->[$l]} = 1;			#*	-- ADD to bad sequences and LOG
		#}
		#
		#
		#for(my $l = 0; $l <= $#{$refbatch}; $l++) {							# LOOP through our batch list
		#	if($seqobj->primary_id =~ /$refbatch->[$l]/) {					# IF we found a seq|acc match
		#		if($seqobj->seq() =~ /(ERROR)/) {							#*IF it's not a qood sequence
		#			$bad_seq->{$refbatch->[$l]} = 1;			#*	-- ADD to bad sequences and LOG
		#			#pring WARNING "$refbatch->[$l] [".$seqobj->primary_id()."] gave a bad sequence.\n";
		#		}
		#		elsif($seqobj->alphabet =~ /$type/) {					#*ELSE it's good and a protein sequence
		#			$batch->{$refbatch->[$l]} = $seqobj->seq();				#*ADD the sequence to our batch and TOC
		#			#$self->{TABLE_OF_CONTENTS}{$refbatch->[$l]} = $seqobj->desc()."\t".$seqobj->primary_id();
		#		} else { # alphabet is not protein							#*ELSE it's not a protein
		#			$bad_seq->{$refbatch->[$l]} = 1;			#* 	-- ADD to bad sequences and LOG
		#			#pring WARNING "$refbatch->[$l] is not a protein sequence.\n";
		#		}
		#		last;														#*EXIT loop since we found our sequence
		#	} else {														# ELSE we did not find our seq|acc match
		#	}
		#}
	}
	
	#print "\t\t".'s8'."\n";
	
	return $batch;
}

sub _open_error {
	my $self = shift;
	open(ERR,">>",$self->{DIR}{THRESHOLD}.$self->{FILES}{ERROR})
		or croak("Error file not opened: $@");
}
sub _open_warn {
	my $self = shift;
	open(WARNING,">>",$self->{DIR}{THRESHOLD}.$self->{FILES}{WARNING})
		or croak("Error file not opened: $@");
}

1;
__END__

