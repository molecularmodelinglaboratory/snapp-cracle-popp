#!/usr/bin/perl -w

# Filename:	Initialize.pm
# Author:	Stephen Bush
# Date:		01.28.07 [v0.2]
# Modified: 03.03.09 [v0.4]

package AIP::Initialize_nt;
require Exporter;

our @ISA      = qw( Exporter AIP::Core );
our @EXPORT   = qw(  ); # insert sub names here
use vars qw( );
our $version  = 1.00;

use strict;
use Carp;
use AIP::Core;
use AIP::Seq;
use AIP::Stats;
use Bio::DB::GenBank;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;	# for blasting the sequence locally
use Switch;

##### Ensure DESTROY is called
#use sigtrap qw(handler the_rest_is_silence INT QUIT);
#sub the_rest_is_silence { die(); } # Calling this on a kill signal ensures DESTROY is called
sub the_rest_is_silence { exit; }
$SIG{'INT'} = \&the_rest_is_silence;
$SIG{'QUIT'} = $SIG{'INT'};

####
# CONSTRUCTOR AND DESTRUCTORS

# NEW
sub new {
	my $class = shift;
	#my $self = shift;   # self from Director
	my $sibling = shift;
    
    my $self = AIP::Core->new();
    $self->{DIR} = $sibling->{DIR};
    $self->{FILES} = $sibling->{FILES};
    $self->{PARAM} = $sibling->{PARAM};
    $self->{TARGET} = $sibling->{TARGET};
    $self->{DATE} = $sibling->{DATE};
    $self->{HANDLE} = $sibling->{HANDLE};
    $self->{FLAGS} = $sibling->{FLAGS};
	$self->{ABOUT} = $sibling->{ABOUT};
    
    $self->{HASH} = {}; #{ a => { b => 1 } };
    $self->{MATRIX} = [];
	
	$self->{CORRFILES} = {
		PRE => '',
		POST => ''
	};
	
	bless($self,$class);
	return $self;
}

# DESTROY
sub DESTROY {
	print 'Initialize done.'."\n";
	return 0;
}


#### CLASS METHODS

## RESET
sub matrix_reset {
	my $self = shift;
	$self->{MATRIX} = shift;
}

## initialize network
sub initialize_network {
	my $self = shift;
	
	if($self->{PARAM}{IS_STRUCTURE}) {
		return $self->initialize_PDUG;
	}
	else {
		return $self->initialize_homology_network;
	}
}

## INITIALIZE HOMOLOGY NETWORK
sub initialize_homology_network {
    my $self = shift;
    my $no_new_sequence_count = 0;	# fail-safe for redundant batch queries
    #my $leaves = {};
    #share($leaves);

	my $blast_temp = '';
	my $hash_temp = {};
	$self->{to_gbseq} = [];
	$self->{to_blast} = [];
	$self->{to_parse} = [];
	$self->{to_blast_sec} = [];
	$self->{seq_data} = {};
	$self->{bad_seq} = {};
	
	print $self->generate_timestamp.'> initializing target hub...'."\n";
	
	
    #### First Blast Query
    eval {
		$blast_temp = $self->sablast(\$self->{TARGET}{SEQ},\$self->{TARGET}{ID}); 
		chmod(0777,$blast_temp); # open permissions													# CHANGE the file permissions
	};
	if($@) {
		croak('Error in BLAST query :: '.$self->{TARGET}{ID}."\n".$@."\n");
	}
	# No BLAST file
	elsif(!-f $blast_temp) {
		croak('No BLAST file returned :: '.$self->{TARGET}{ID}."\n");
	}
	# Empty BLAST file
	elsif(!-s $blast_temp) {
		unlink($blast_temp);
		croak('Empty BLAST file returned :: '.$self->{TARGET}{ID}."\n");		
	}
	else {}
    
    #### Parse Blast Query
    eval {
		$hash_temp = $self->blastParse($blast_temp,$self->{TARGET}{ID});
    };
	if($@) {
		croak('Error in BLAST results :: '.$self->{TARGET}{ID}."\n".$@."\n");
	}
	elsif(!(scalar keys %{$hash_temp})) {
		croak('No BLAST results :: '.$self->{TARGET}{ID}."\n".
			'Empty query results?'."\n");
	}
	else {}
	
	## We have results
	# push the reference to the array keys onto the seq queue
	# we'll pull it off and use the entire set of sequences at once.
	$self->{HASH}{$self->{TARGET}{ID}} = $hash_temp;
	push(@{$self->{to_gbseq}},keys %{$hash_temp});
	
	## Get the sequence data
	# if we already have the sequence data for the first few, we can begin
	# blasting right away
	#my $seq_temp = shift(@{$self->{to_gbseq}});
	my $tseq;
	eval {
		$tseq = batch2seq(
			$self->{to_gbseq},
			$self->{bad_seq},
			$self->{PARAM}{IS_PROTEIN});
	};
	if($@ || !$tseq || $tseq == 1) {
		croak('No sequence data!'.$@."\n");
	}
	else {
		map {
			$self->{seq_data}{$_} = $tseq->{$_}[0];
			$self->{ABOUT}{$_} = $tseq->{$_}[1];
		} keys %{$tseq};
		$self->{to_gbseq} = [];
	}

	## Copy the sequences
	## Add the keys to to_blast queue
	#map { print 'hey'.$_."\n" } keys %{$seq_temp};
	#print '2...'.$self->{to_gbseq}->pending."\n";
	#print '3...'.$self->{to_blast}->pending."\n";
	#$self->{seq_data} = shared_clone($seq_temp);
	push(@{$self->{to_blast}},keys %{$self->{seq_data}});
	#print '4...'.$self->{to_blast}->pending."\n";
	#my $tcount = 0;
	#map { print $_.'::'.$seq_temp->{$_}."\n"; die if(++$tcount == 5); } keys %{$seq_temp};
	
	$self->blast_thread();

#Up to this point:
#1) BLASTed our target
#2) Parsed our target
#3) Iterativaly BLASTed subsequent results
#    -- using threads (max 4 at a time)
#    -- bottleneck -- performs 4 and waits, repeat
#4) Repeat for up to <iterations> times
#

    ##### Write data
	my $file = $self->{DIR}{DATA}.$self->{FILES}{H_HASH};
	my $list = [$self->{TARGET}{ID}];
	$self->_output_hash($file,$self->{HASH},$list);
    
    return $self->{HASH};															# RETURN our haip reference
}

sub blast_thread {

	my $self = shift;
	my $seq = {};
    my $work_to_do = $#{$self->{to_blast}};
    my $tried_seq = {};
    my $counter = 0;
	
    print $self->generate_timestamp.'> initializing past target...'."\n";
    
    # Thread 'cancellation' signal handler
    #$SIG{'KILL'} = sub { threads->exit(); };

	my $no_seq = 0;
	my $num_seqs = '';
	# Create the network
	while(
        (
            scalar @{$self->{to_gbseq}} ||
            scalar @{$self->{to_parse}} ||
            scalar @{$self->{to_blast}} 
        )# && 
        #$self->{FLAG}{COUNT} <= $self->{PARAM}{NUMITERATIONS})
	) {
		#print '1..'."\n";
		
		last if($counter >= $self->{PARAM}{NUMITERATIONS});
		
		## SEQUENCES
		if($#{$self->{to_gbseq}} > 100) {
			#$work_to_do--;
            #print $tseq."\n"; die;
            my $tseq;
			
			#$num_seqs .= $#{$self->{to_gbseq}}.'.';
			#print $self->generate_timestamp.'> '.$num_seqs."\n";
			
			#print '2'."\n";
            
			eval{
				#print $self->generate_timestamp.'> sequences...'.threads->tid."\n";
				$tseq = batch2seq(
					$self->{to_gbseq},
					$self->{bad_seq},
					$self->{PARAM}{IS_PROTEIN}
				);
			};
			if($@ || !$tseq || $tseq == 1) {
				
				#print "\t".'a..'."\n";
				
				carp('Failed to retrive sequence data! '."\n".$@."\n");
                
				map { print $_."\t"; } @{$self->{to_gbseq}};
				
                # trying again may be extraneous...batch tries up to 5 times
                my $tstr = "$tseq";
                if(!exists $tried_seq->{$tstr}) {
                    $tried_seq->{$tstr} = 1;
                    push(@{$self->{to_gbseq}},$tseq);
                }
			}
			else {
				#print "\t".'b..'."\n";
                map {
					$self->{seq_data}{$_} = $tseq->{$_}[0];
					$self->{ABOUT}{$_} = $tseq->{$_}[1];
				} keys %{$tseq};
                
                push(@{$self->{to_blast}},keys %{$tseq});
                #$work_to_do += scalar keys %{$tseq};
            }
			$self->{to_gbseq} = [];
		}
		
		## PARSE
		# $tpar is the file name to parse
		while(my $tpar = shift(@{$self->{to_parse}})) {
			#$work_to_do--;
            my $tph;
				#print '3'."\n";
            
			## Pull ID from filename
			# <id>_<evalue>.blast; NP_11111_0.0001.blast
			$tpar =~ /^(?:BLAST\/)?(.*)\_\d*?\.\d*?\.blast$/;
			my $tpar_id = $1;
			
			# Parse the file
			eval {
				#print $self->generate_timestamp.'> parsing...'.threads->tid."\n";
				$tph = $self->blastParse($tpar,$tpar_id);
			};
			if ($@) {
				#print "\t".'a..'."\n";
				carp('Failed to parse '.$tpar."\n".$@."\n");
			}
			else {
				#print "\t".'b..'."\n";
				my $to_get = {};
				map {
					delete($tph->{$_}) if(exists $self->{bad_seq}{$_});
					$to_get->{$_} = 1
						if(!exists $self->{HASH}{$_}
						 && !$self->blast_file_exists($_)
						 && !exists $self->{seq_data}{$_});
				} keys %{$tph};
				
				# Save the parse results
				if(!exists $self->{HASH}{$tpar_id}) {
					$self->{HASH}{$tpar_id} = {};
				}
				else {
					# we shouldn't be here...hopefully
					# this is adding the results of a BLAST query
					# should only be done when we don't have
					# the info...
				}
				$self->{HASH}{$tpar_id} = $tph;
				push(@{$self->{to_gbseq}},keys %{$to_get});
				#$work_to_do++;
			}
		}
		
		## BLAST
		if(my $tb = shift(@{$self->{to_blast}})) {
            #$work_to_do--;
				#print '4'."\n";
			
			# BLAST the sequence
			my $bf;
			eval{
				#print $self->generate_timestamp.'> blasting...'.threads->tid."\n";
				#print $counter.' :: '.$work_to_do.' :: '."\n";
				$DB::single=2; # insert at line 9!
				$bf = $self->sablast(\$self->{seq_data}{$tb},\$tb);
				chmod(0777,$bf);					# CHANGE permissions
			};
			# Add results to be parsed
			if($@ || !-f $bf || !-s $bf) {
				#print "\t".'a..'."\n";
				carp('No BLAST data returned: '.$tb."\n".$@."\n");
				#lock($self->{bad_seq});
				$self->{bad_seq}{$tb} = 1;
				$no_seq++;
				#print $ENV{PATH}."\n";
				sleep(5); # prevent ncbi from screaming
				#print $bf.'...bad'."\n";
				#{ lock($counter); $counter--; }
			}
			#elsif(!-s $bf) {
			#	unlink $bf;
			#	lock($self->{bad_seq});
			#	$self->{bad_seq}{$tb} = 1;
			#}
			else {
				#print "\t".'b..'."\n";
				push(@{$self->{to_parse}},$bf);
				delete($self->{seq_data}{$tb});
                #$work_to_do++;
				$counter++;
				$no_seq = 0;
			}
		}
	}
				#print '5'."\n";
	
    #print 'Leaving '.threads->tid.' :: '.
    #    $self->{to_gbseq}->pending.' '.
    #    $self->{to_parse}->pending.' '.
    #    $self->{to_blast}->pending.' '.
    #    $work_to_do."\n";
	
	return 1;
}

## INITIALIZE PDUG
sub initialize_PDUG {
	my $self = shift;
	
	my $dali_file = 'dali.detox';
	my $regex = qr/([0-9a-zA-Z]{4})([0-9a-zA-Z]?)\,([0-9a-zA-Z]{4})([0-9a-zA-Z]?)\,([0-9\.]*)\;/;
	
	open(my $fh,'<',$dali_file) or die('Couldn\'t open file'."\n");
	while(my $line = <$fh>) {
		$line =~ /$regex/;
		my $p1 = $2 ? uc($1.$2) : uc($1).'-';
		my $p2 = $4 ? uc($3.$4) : uc($4).'-';
        $self->{HASH}{$p1} = {} if(!exists $self->{HASH}{$p1});
        $self->{HASH}{$p2} = {} if(!exists $self->{HASH}{$p2});
        $self->{HASH}{$p1}{$p2} = $5;
        $self->{HASH}{$p2}{$p1} = $5;
	}
	
	return 0;
}

## NORMALIZE NETWORK
sub normalize_network() {
	my $self = shift;															# GET our class data
	
	# count hashes
	my $treeNumSeq = 0;															# SCALAR count number of total seq
	my $queryNumSeq = 0;														# SCALAR count number of query seq
	my %seqNumSeq = ();															# HAIP count number of seq appearances
	
	# average haipes
	my $tree_avg = 0;															# SCALAR average of total sequences
	my %query_avg = ();															# HAIP average per query
	my %seq_avg = ();															# HAIP average per sequence
	
	my $max = 0;																# SCALAR initiate a max value
	my $min = 0;																# SCALAR initiate a min value
	
	
	##### get and calculate averages
	#####	-- except the sequence averages
	foreach my $branch (keys(%{$self->{HASH}})) {								# LOOP through our haip -- 1st LVL
		$query_avg{$branch} = 0;												# RESET the query average
		$queryNumSeq = 0;														# RESET the query count
		foreach my $leaf (keys(%{$self->{HASH}{$branch}})) {					# LOOP through our haip -- 2nd LVL
			$tree_avg += $self->{HASH}{$branch}->{$leaf};						# ADD score to total average
			$query_avg{$branch} += $self->{HASH}{$branch}->{$leaf};			# ADD score to query average
			if(!exists($seq_avg{$leaf})) { 									# IF we don't have an average for this seq
				$seq_avg{$leaf} = 0; 											#	-- INITIATE seq average
				$seqNumSeq{$leaf} = 0; }										#	-- INITIATE seq count
			$seq_avg{$leaf} += $self->{HASH}{$branch}->{$leaf};				# ADD score to seq average
			$treeNumSeq++;														# INCREMENT total count
			$queryNumSeq++;														# INCREMENT query count
			$seqNumSeq{$leaf}++;												# INCREMENT sequence count
		}
		if($queryNumSeq == 0) {													# IF the query had no data (should have caught already)
			print "Bad Sequence '$branch'-> No BLAST results in hash!!!";	#	-- LOG in ERROR file
			delete $self->{HASH}{$branch};									#	-- DELETE the branch (useless)
			next;																#	-- SKIP ahead to the next branch
		}
		eval{																	# EVALUATE to make sure we haven't missed it
			$query_avg{$branch} = $query_avg{$branch} / $queryNumSeq;			# CALCULATE the query average
		}; if($@) {																# IF we had an error (i.e. divide by 0)
			print {$self->{HANDLE}{LOG}} "$branch has no leaves.\n";								#	-- LOG in ERROR file
			delete $self->{HASH}{$branch};									#	-- DELETE the branch
		}
	}
	$tree_avg = $tree_avg / $treeNumSeq;										# CALCULATE tree average
	foreach my $k (keys(%seq_avg)) {											# LOOP through our valid sequences
		$seq_avg{$k} = $seq_avg{$k} / $seqNumSeq{$k};							#	-- CALCULATE the sequence average
	}
	
	##### calculate logodds
	my $a = 0; my $b = 0;														# INITIATE two helper variables
	foreach my $i (keys(%{$self->{HASH}})) {								# LOOP through 1st tier
		foreach my $j (keys(%{$self->{HASH}{$i}})) {					# LOOP through 2nd tier
			if(!exists($seq_avg{$i})) { 									# IF this sequence wasn't returned (prob not in DB) 
				$seq_avg{$i} = $query_avg{$i};							#	-- USE it's query avg
				#pring WARN "$i was not found in any BLAST queries.\n"; }		#	-- LOG a warning
			}
			$a = $self->{HASH}{$i}{$j} / $query_avg{$i};			# CALCULATE ratio of query score to query average
			$b = $seq_avg{$i} / $tree_avg;									# CALCULATE ratio of seq avg to total avg
			$self->{HASH}{$i}{$j} = log($a/$b) / log(2);			# CALCULATE log odds score
			if($self->{HASH}{$i}{$j} > $max) {						# IF this score is the max as of yet
				$max = $self->{HASH}{$i}{$j};						#	-- SAVE it as such
			} elsif($self->{HASH}{$i}{$j} < $min) {					# IF this score is the min as of yet
				$min = $self->{HASH}{$i}{$j};						#	-- SAVE it as such
			}
		}
	}

	$self->{HASH}{MAX} = $max;												# ADD max value to tree
	$self->{HASH}{MIN} = $min;												# ADD min value to tree
	
	##### Write data
	my $file = $self->{DIR}{DATA}.$self->{FILES}{H_LODS};
	my $list = [$self->{TARGET}{ID}];
	my $meta = ['MAX','MIN'];
	$self->_output_hash($file,$self->{HASH},$list,$meta);													# CLOSE output file
	
	return $self->{HASH};														# RETURN ref to our tree
}

## 
sub create_matrix() {
	my $self 	= shift;														# GET our class data
	my @header  = (); # shouldn't matter...right?															# INITIATE array
	my $grid	= { HEADER => '', DATA => '' };									# INITIATE matrix output string
	
	my $max = $self->{HASH}{MAX};										# GET max value and MODIFY to emmulate perfect score
	my $min = $self->{HASH}{MIN} * 2;											# GET min value and MODIFY to emmulate no score
	delete $self->{HASH}{MAX};												# REMOVE max from tree
	delete $self->{HASH}{MIN};												# REMOVE min from tree
	
	##### get header line
	push @header, $self->{TARGET}{ID};											# PUSH seed onto header as first value
	foreach my $branch (keys(%{$self->{HASH}})) {								# LOOP through our tree
		##print $branch.'..pushed'."\n";
		if($branch ne $self->{TARGET}{ID}) {										# IF the branch is not our seed
			push @header, $branch;												# PUSH branch onto header
		}
	}

	my $maxs = sprintf('%.5f',$max);
	my $mins = sprintf('%.5f',$min);

	##### arrange data lines
	for(my $i = 0; $i <= $#header; $i++) {										# LOOP through 1st tier
		$grid->{HEADER} .= $header[$i]."\t";										# ADD current branch to output string header
		$self->{MATRIX}[$i] = [];
        for(my $j = 0; $j <= $#header; $j++) {									# LOOP through 2nd tier
			#if($i eq $j) {														# IF we're at a diagonal value (mirror)
			#	$grid->{DATA} .= $maxs;											#	-- ADD the max to output string data
            #    $self->{MATRIX}[$i][$j] = $max;
            #} elsif(exists($self->{HASH}{$header[$i]}{$header[$j]})) {		# IF we're at a valid data point
            if(exists($self->{HASH}{$header[$i]}{$header[$j]})) {		# IF we're at a valid data point
				$grid->{DATA} .= sprintf('%.5f',$self->{HASH}{$header[$i]}{$header[$j]});	# 	-- ADD value to the output string data
                $self->{MATRIX}[$i][$j] = $self->{HASH}{$header[$i]}{$header[$j]};
            }
			elsif($i eq $j) {
				$grid->{DATA} .= $maxs;											#	-- ADD the min to output string data
                $self->{MATRIX}[$i][$j] = $max;
			}else {															# ELSE we have no value
				$grid->{DATA} .= $mins;											#	-- ADD the min to output string data
                $self->{MATRIX}[$i][$j] = $min;
            }
			$grid->{DATA} .= "\t";												# ADD delimiting tab to output string
		} 
		##### NOTE: we remove these extra tabs and newlines to avoid
		##### extra rows or columns being added by accident when the file
		##### is read in and split. Otherwise, our file would have increasingly
		##### more newlines and tabs at the end of each one for every iteration
		##### we ran on the matrix		
		$grid->{DATA} =~ s/\s*$//; 												# REMOVE extra tabs and newlines from row
		$grid->{DATA} .= "\n";													# REMOVE a newline on
	} 
	$grid->{DATA} =~ s/\s*$//;													# REMOVE extra tabs and newlines from graph
	$grid->{HEADER} =~ s/\s*$//;												# REMOVE same from header
	unshift(@{$self->{MATRIX}},\@header);

	chomp($grid);
	
	open(my $mfh,'>'.$self->{DIR}{DATA}.$self->{FILES}{H_MATRIX});				# OPEN filehandle for the matrix file
	print $mfh $grid->{HEADER}."\n".$grid->{DATA}."\n";						# PRINT out our header & data string
	close($mfh);																# CLOSE the out file
	
	return $self->{FILES}{H_MATRIX};												# RETURN the filename
}

## CORRELATE
sub correlate() {
	my $self = shift;															# GET our class data

	## The problem
	# 1) Calling R from a cluster is a huge pain -- have to sleep and check for file
	#	results. Not very efficient
	# 2) Unfortunately, we cannot connect to the internet when submitted to the cluster
	# 3) Our correlation still has some errors
	# 4) Suck it up and use R again

	#my $header = shift(@{$self->{MATRIX}});

	#print 'c '.$self->{MATRIX}[0][0].' :: ';

	# Create R commands
	my $rcmd = 'infile=read.table(\''.$self->{DIR}{DATA}.$self->{FILES}{H_MATRIX}.
		'\', header=TRUE, sep=\'\t\');'.
		'a=cor(x=infile,use=\'pairwise.complete.obs\',method=\'pearson\');'.
		'write.table(a,file=\''.$self->{DIR}{DATA}.$self->{FILES}{H_PCORR}.
		'\', sep=\'\t\','.'row.names=FALSE,col.names=TRUE,quote=FALSE);';
	
	# Write R script
	my $rfile = $self->generate_timestamp().'.r';
	open(my $fh,'>',$rfile);
	print $fh $rcmd;
	close($fh);
	
	# Execute R script
	my @args = ('R','--vanilla --slave -f '.$rfile); # INITIATE an array with our command line options
	eval{
		system(@args) == 0															# CALL R from the command line
			or croak "system @args failed: $?";
	};
	if($@) {
		carp(' :: Failed call to R. '.$@);
	}
	
	#print 'c '.$self->{MATRIX}[0][0]."\n";
	
	# Get R script results
	$self->{MATRIX} = $self->_input_matrix($self->{DIR}{DATA}.$self->{FILES}{H_PCORR});
	

	
	#unshift(@{$self->{MATRIX}},$header);
	
#	print threads->tid.' '.$self->{MATRIX}[0][0]."\n";
	unlink($rfile);
	
	return $self->{MATRIX};
	
	#map {
	#	map {
	#		print $_.' ';
	#	} @{$_};
	#	print "\n";
	#} @{$self->{MATRIX}};

### OLD :: USE OUR CORRELATION
#	print  'a'.$#{$self->{MATRIX}}.' '.$#{$self->{MATRIX}->[0]}."\n";
#	my $header = shift(@{$self->{MATRIX}});
#    $self->{MATRIX} = pcoeff_matrix($self->{MATRIX});
#	$self->_output_matrix($self->{DIR}{DATA}.$self->{FILES}{H_PCORR1},$self->{MATRIX});
#
#	## - 2 :: one to count for header, one to start at 0
#	my $size = scalar @{$self->{MATRIX}} - 1;
#	my $i = 0;
#	my $graph = '';
#	my $head = '';
#	for my $i (0 .. $size) {
#		$head .= $header->[$i].' ';
#		foreach my $j (0 .. $size) {
#			#if($i != $j) {
#				$graph .= sprintf('%.4f ',$self->{MATRIX}[$j][$i]);
#			#}
#			#elsif($i < $size) {
#			#	$graph .= $header->[$i+1]."\n";
#			#}
#			#else {}
#		}
#		$graph .= "\n";
#	}
#	unshift(@{$self->{MATRIX}},$header);
#	print  'b'.$#{$self->{MATRIX}}.' '.$#{$self->{MATRIX}->[0]}."\n";
#
#	open(my $cfm,'>',$self->{DIR}{DATA}.$self->{FILES}{H_PCORR1});
#	print $cfm $head."\n".$graph;
#	close($cfm);
#
#	return $self->{MATRIX};
	
#Statistics::R is outdated
#    
#    &R::initR("--silent");
#    my $r_cmd = 'infile=read.table(\''.
#        $self->{DIR}{DATA}.$self->{FILES}{H_MATRIX}.
#        '\', header=TRUE, sep=\'\t\');'.
#        'a=cor(x=infile,use=\'pairwise.complete.obs\',method=\'pearson\');'.
#        'write.table(a,file=\''.
#        $self->{DIR}{DATA}.$self->{FILES}{H_PCORR1}.'\', sep=\'\t\','.
#        'row.names=FALSE,col.names=TRUE,quote=FALSE);';
#    &R::eval($r_cmd);
#
#	return $self->{DIR}{DATA}.$self->{FILES}{H_PCORR1};					# RETURN our correlation file name



#infile=read.table(\''.
#        $self->{DIR}{DATA}.$self->{FILES}{H_MATRIX}.
#        '\', header=TRUE, sep=\'\t\');'.
#        'a=cor(x=infile,use=\'pairwise.complete.obs\',method=\'pearson\');'.
#        'write.table(a,file=\''.
#        $self->{DIR}{DATA}.$self->{FILES}{H_PCORR1}.'\', sep=\'\t\','.
#        'row.names=FALSE,col.names=TRUE,quote=FALSE);';


}

#sub init_test {
#	my $self = shift;
#	print threads->tid().' '.$self->{MATRIX}.' '.$self->{MATRIX}[0][0]."\n";
#}

sub update_filenames_with_threshold {
	my $self = shift;
	my $lower = shift;
	my $upper = shift;
	
	$lower *= 100;
	$upper *= 100;
	
	$self->{FILES}{H_MATRIX} = 'h04_ident_t'.$lower.'-'.$upper.'.out';
	$self->{FILES}{H_PCORR} = 'h05_pcorr_t'.$lower.'-'.$upper.'.out';
	
	return 0;
}

## FIND_IDENTITY
sub find_identity {
	my $self = shift;
	my $lower = shift;
    my $upper = shift;
	my $grid = '';					# initiate an output string

	# Note: This function apparently used to calculate the total
	# number of times a node appeared as a positive identity.
	# I don't remember why it did this, so I took it out.
	# Refer to the old code to get it back.

	my $header = shift(@{$self->{MATRIX}});

	my $positives = 0;
	my $identity;
	my $dimension = $#{$header};
	$grid = join("\t",@{$header})."\n";
	for(my $i = 0; $i <= $#{$self->{MATRIX}}; $i++) {
		$identity = [];
		for(my $j = 0; $j <= $#{$self->{MATRIX}[$i]}; $j++) {
			if($self->{MATRIX}[$i][$j] > $lower &&
				$self->{MATRIX}[$i][$j] <= $upper) {
				$identity->[$j] = 1;
				$positives++;
			}
			else {
				$identity->[$j] = 0;
			}
		}
		$grid .= join("\t",@{$identity})."\n";
	}
	#print $#{$self->{MATRIX}}."\n";
	#print $self->{MATRIX}[0][0].' :: ';
	unshift(@{$self->{MATRIX}},$header);
	#print $self->{MATRIX}[0][0]."\n";
	#print $#{$self->{MATRIX}}."\n";
	chomp($grid);
#die;
	
	if(!$positives) { return 0; }
	
	##### output the truth table to a file
	$upper *= 100;
	$lower *= 100;
	my $idf = $self->{DIR}{DATA}.$self->{FILES}{H_MATRIX}.'_t'.$lower.'-'.$upper.'.out';
	open(my $fh,'>'.$self->{DIR}{DATA}.$self->{FILES}{H_MATRIX})						# open the file or die
		or croak("Cannot open file.\n\t$@");
	print $fh $grid;
	close($fh);													# close the file
	
	return $self->{MATRIX};												# return the outfile name
}

## SABLAST (STAND ALONE BLAST)
sub sablast( $ $ ) {
	my $self = shift;															# GET our class data
	my $seq = shift or croak("No sequence data.\n"); # GET sequence reference
	my $acc = shift or croak("No accession number.\n"); # GET accession reference
	#print 'HERE!!'.$$acc.' '.$$seq."\n";
	#my $blast_file = $self->{DIR}{BLAST}.$$acc.'_'.$self->{PARAM}{EVALUE}.'.blast';
	my $blast_file = $self->blast_file_name($$acc);
    #print 'hi.'.$blast_file."\n";
	#if(-f $blast_file && !-s $blast_file) {
	#	unlink($blast_file);
	#}
    
	if(!-f $blast_file || !-s $blast_file) {
		my @params = (																# INITIATE array for BLASTALL parameters
			'program' => $self->{PARAM}{PROGRAM},
			'database' => $self->{PARAM}{DB},
			'outfile' => $blast_file,
			'expect' => $self->{PARAM}{EVALUE},
			'a' => 1
		);
		print 'BLAST: '.$$acc."\n";
		my $factory = Bio::Tools::Run::StandAloneBlast->new(@params); # INITIATE SAB handler
		my $seqobj = Bio::Seq->new(-id => $$acc,-seq => $$seq); # INITIATE a Seq obj with our seed data
		$factory->blastall($seqobj); # RUN our BLAST query
	}
    #else { print 'already done '.$$acc."\n"; }
    #print {$self->{HANDLE}{LOG}} $$acc.'::'.$blast_file.'::'.$$seq."\n\n";
    if(!-s $blast_file) {
        unlink($blast_file);
        return 0;
    }
    
    
	return $blast_file; # RETURN the filename of our BLAST data
}

sub blast_file_name {
    my $self = shift;
    my $id = shift;
   	my $filename = $self->{DIR}{BLAST}.$id.'_'.$self->{PARAM}{EVALUE}.'.blast';

    return $filename;
}

## blast file name
sub blast_file_exists {
	my $self = shift;
	my $filename = shift;
	$filename = $self->{DIR}{BLAST}.$filename.'_'.$self->{PARAM}{EVALUE}.'.blast';
	
	if(-f $filename) {
		return 1;
	}
	else {
		return 0;
	}
}

## BLAST PARSE
sub blastParse($) {
    my $self = shift;															# GET our class data
    my $infile = shift;
	my $parid = shift;
    my $scores = {};															# INITIATE an array ref for scores
    my $count = 0;																# INITIATE a count for total no. of scores
    my $acc = 0; my $bs = 0; my $ev = 0; # INITIATE important BLAST variables
    my $prefix = $self->{PARAM}{IS_PROTEIN} ? 'NP' : 'NT';
    my $save = 0;
	my $resave = 0;
	my $line = '';
	my $out = '';
	
    #print 'hi :: '.$self->{PARAM}{IS_PROTEIN}.'::'.$self->{PARAM}{IS_STRUCTURE}.'::'.$prefix."\n";
    
    ##### read in the BLAST file data
    open(my $BF,$infile);															# OPEN the blast file
    my @blast_file = <$BF>;														# READ in entire file
    close($BF);
    # CLOSE the file handle
    if($blast_file[0] =~ /\>SAVED\:/) {
		$line = shift(@blast_file);
		$line =~ /^\>SAVED\:(.*)\s+$/;
		$self->{ABOUT}{$parid} = $1;
	}
	else {
		splice(@blast_file,0,20); # DELETE the unnecessary file data
		$resave = 1;
	}
	
    my $np = qr/^ref\|((?:$prefix)_\d{3,})(?:[\w\W[:punct:]]*)*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $gb = qr/^(?:\w+\|)?(?:\w+\d*_\w+)?([a-zA-Z]+\d{3,})\s*(?:[\w\W[:punct:]]*)*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $sp = qr/^[\w\d\_]*\s\(([a-zA-Z]+[\w\d]*)\).*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $junk = qr/^.*\(?(\w+\d{3,})\)?.*([\d\.]+)\s+(\d+(e-|\.)\d*)/;
	
	##### parse through the data
    foreach $line (@blast_file) { # LOOP through the file
        if($line =~ m/^>/ || $line =~ m/Matrix\:/) { # IF we're past the initial data
            @blast_file = {};													#	-- DELETE the remainder of the file
            last;																#	-- EXIT the loop
        }
        if($line =~ /(hypothetical)/ && !$self->{PARAM}{HYPOTHETICAL}) { # IF we don't want the hypothetical data
            next; }																#	-- SKIP it and reiterate the loop
            
        ### NOTE: these regular expressions will pick out any and every
        ### valid sequence accession number, the bit score, and the e-value
        ### for every line.
        ### The first expression picks almost every line out
        ### The second picks out lines with an accession right at the front
        ### The third should pick out any leftovers we might want
        ### The fourth seems only to pick out blank lines
        ### I have never seen the fifth needed
		## NOTE (03.31.09): changed 3rd to check for [a-zA-Z] -- \w picks up
		## numbers as well as letters.
        elsif($line =~ /$np/) {
            $acc = $1;															# SAVE accession number
            $bs = $2;															# SAVE raw score
            $ev = 3;
			$save = 1;
			#print 'a....'.$1.':'.$2.':'.$3."\n";
        } elsif(!$self->{PARAMS}{REFSEQ_ONLY}) {
            if($line =~ /$gb/) {
                $acc = $1;															# SAVE accession number
                $bs = $2;															# SAVE raw score
                $ev = $3;															# SAVE e-value	
				$save = 1;
                #print 'b....'.$1.':'.$2.':'.$3."\n";
            } elsif($line =~ /$sp/) {
                $acc = $1;															# SAVE accession number
                $bs = $2;															# SAVE raw score
                $ev = $3;															# SAVE evalue
				$save = 1;
                #print 'c....'.$1.':'.$2.':'.$3."\n";
            }
        } elsif($line =~ /$junk/) {
            #print {$self->{HANDLE}{LOG}} 'NOT USED: '.$line; 
        } elsif($line =~ /\s*/) {
            # LOG an error
        } else {
            carp('Not sure how we got here: '.$line); # LOG an error
        }
        
        # saving the good data or printing to error file
        if(!$acc || !$ev || !$bs) {
            next; # IF we're missing data SKIP
        } elsif($acc =~ m/XP_/i) { # IF it's an experimental sequence
            carp('NOT USED, experimental: '.$line);	#-- LOG a warning
            next; #-- SKIP
        } elsif(exists($self->{bad_seq}{$acc})) {
            next; # IF it's a known bad seq, SKIP it
        } elsif($self->{PARAM}{SCORE_TYPE} eq 'bitscore') { # IF we want to save the bitscore
            $scores->{$acc} = $bs; # -- ADD score to our array
            $count++; # -- INCREMENT our score count
        } elsif($self->{PARAM}{SCORE_TYPE} eq 'evalue') { # IF we want to save the evalue
            $scores->{$acc} = $ev; # -- ADD score to our array
            $count++; 															#	-- INCREMENT our score count
        } else {}
		
		if($resave && $save) {
			$out .= $line;
		}
    }
    
    #print $count.'::'.$scores."\n";
    if($resave && $out) {
		open(my $fh,'>',$infile);
		print $fh '>SAVED:'.$self->{ABOUT}{$parid}."\n".$out;
		close($fh);
	}
    
    $count > 0 ? return $scores : return 0; # IF we found scores, RETURN them, ELSE RETURN 0
}

1;

__END__

=head1 NAME

=head1 METHODS

=item NORMALIZE NETWORK

I: none
O: haip reference
Desc: Calculates the log odds score for a given haip

	log( ( query_score / query average )
		 -------------------------------
		 ( acc_average / total average ) )
	________________________________________
	
				log( 2 )
				
e.g. for [query A=>B] in:
	query A		query B		query C
	B:1			A:1			D:1
	c:2			C:2			B:2
	D:3			D:4			A:3

	log(   ( 1 / 2 )
		 ---------------
		 ( 1.5 / 2.111 )   )
	_________________________
	
		  log( 2 )


=cut


