#!/usr/bin/perl -w

# Filename:	Initialize.pm
# Author:	Stephen Bush
# Date:		01.28.07 [v0.2]
# Modified: 03.03.09 [v0.4]

package AIP::Initialize;
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
use threads;
use threads::shared;
use Thread::Queue;

use sigtrap qw(handler the_rest_is_silence INT QUIT);
sub the_rest_is_silence { die(); } # Calling this on a kill signal ensures DESTROY is called

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
    
    $self->{HASH} = &share({}); #{ a => { b => 1 } };
    $self->{MATRIX} = &share([]);
	
	$self->{CORRFILES} = {
		PRE => '',
		POST => ''
	};
	
	bless($self,$class);
	return $self;
}

# DESTROY
sub DESTROY {
	print 'Initialize '.threads->tid.' done.'."\n";
	return 0;
}


#### CLASS METHODS

## INITIALIZE HOMOLOGY NETWORK
sub initialize_homology_network {
    my $self = shift;
    my $no_new_sequence_count = 0;	# fail-safe for redundant batch queries
    #my $leaves = {};
    #share($leaves);

	my $blast_temp = '';
	my $hash_temp = {};
	$self->{to_gbseq} = Thread::Queue->new();
	$self->{to_blast} = Thread::Queue->new();
	$self->{to_parse} = Thread::Queue->new();
	$self->{to_blast_sec} = Thread::Queue->new();
	$self->{seq_data} = &share({});
	$self->{bad_seq} = &share({});
	
	print $self->generate_timestamp.'> first blast...'."\n";
	
	
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
    
	print $self->generate_timestamp.'> first parse...'."\n";
	
    #### Parse Blast Query
    eval {
		$hash_temp = $self->blastParse($blast_temp);
    };
	if($@) {
		croak('Error in BLAST results :: '.$self->{TARGET}{ID}."\n".$@."\n");
	}
	elsif(!(scalar keys %{$hash_temp})) {
		croak('No BLAST results :: '.$self->{TARGET}{ID}."\n".
			'Empty query results?'."\n");
	}
	
	print $self->generate_timestamp.'> first batch...'."\n";
	
	## We have results
	# push the reference to the array keys onto the seq queue
	# we'll pull it off and use the entire set of sequences at once.
	$self->{HASH}{$self->{TARGET}{ID}} = shared_clone($hash_temp);
	$self->{to_gbseq}->enqueue([keys %{$hash_temp}]);
	
	## Get the sequence data
	# if we already have the sequence data for the first few, we can begin
	# blasting right away
	my $seq_temp = $self->{to_gbseq}->dequeue();
	eval {
		$seq_temp = batch2seq(
			$seq_temp,
			$self->{bad_seq},
			$self->{TARGET}{TYPE});
	};
	if($@) {
		croak('No sequence data!'.$@."\n")
	}
	
	print $self->generate_timestamp.'> begin threading...'."\n";
	
	## Copy the sequences
	## Add the keys to to_blast queue
	#map { print 'hey'.$_."\n" } keys %{$seq_temp};
	#print '2...'.$self->{to_gbseq}->pending."\n";
	#print '3...'.$self->{to_blast}->pending."\n";
	$self->{seq_data} = shared_clone($seq_temp);
	$self->{to_blast}->enqueue(keys %{$seq_temp});
	#print '4...'.$self->{to_blast}->pending."\n";
	#my $tcount = 0;
	#map { print $_.'::'.$seq_temp->{$_}."\n"; die if(++$tcount == 5); } keys %{$seq_temp};

	## Start our threads!
	my $t1 = threads->create(
		sub{ $self->blast_thread() } );
	my $t2 = threads->create(
		sub{ $self->blast_thread() } );
	my $t3 = threads->create(
		sub{ $self->blast_thread() } );
	my $t4 = threads->create(
		sub{ $self->blast_thread() } );
	
    while(
        !$t1->is_joinable() ||
        !$t2->is_joinable() ||
        !$t3->is_joinable() ||
        !$t4->is_joinable()
    ) {
        sleep(2);
        #print '1 ready'."\n" if $t1->is_joinable();
        #print '2 ready'."\n" if $t2->is_joinable();
        #print '3 ready'."\n" if $t3->is_joinable();
        #print '4 ready'."\n" if $t4->is_joinable();
        #sleep(2);
    }
    
    #print 'joining threads...'."\n";
	$t1->join();
	$t2->join();
	$t3->join();
	$t4->join();
	#print 'threads joined...'."\n";

#Up to this point:
#1) BLASTed our target
#2) Parsed our target
#3) Iterativaly BLASTed subsequent results
#    -- using threads (max 4 at a time)
#    -- bottleneck -- performs 4 and waits, repeat
#4) Repeat for up to <iterations> times
#

    ##### Write data
    #$self->table_of_contents();												# CREATE table of contents
    my $temp_out = '#'.$self->{TARGET}{ID}."\t".keys(%{$self->{HASH}{$self->{TARGET}{ID}}})."\n";	# print seed acc no.
    foreach my $k (keys(%{$self->{HASH}{$self->{TARGET}{ID}}})) {		    # PRINT out our seed data first
        $temp_out .= $k.': '.$self->{HASH}{$self->{TARGET}{ID}}{$k}."\n";
    }
    
    foreach my $i (keys(%{$self->{HASH}})) {								# rinse and repeat for each other sequence
        if($i ne $self->{TARGET}{ID}) {
            $temp_out .= '#'.$i."\t".scalar(keys(%{$self->{HASH}{$self->{TARGET}{ID}}}))."\n";
            foreach my $j (keys(%{$self->{HASH}{$i}})) {
                $temp_out .= "$j: ".$self->{HASH}{$i}{$j}."\n";
            }
        }
    }
    
	$self->_output($self->{DIR}{DATA}.$self->{FILES}{H_HASH},\$temp_out);
	#open(my $hfh,'>'.$self->{DIR}{DATA}.$self->{FILES}{H_HASH});				# OPEN the output file
    #print $hfh $temp_out;
    #close($hfh);																	# CLOSE our output file
    
    return $self->{HASH};															# RETURN our haip reference
}

sub blast_thread {

	my $self = shift;
	my $seq = &share({});
    my $work_to_do :shared = $self->{to_blast}->pending;
    my $tried_seq = &share({});
    my $counter :shared = 0;
	
    #print 'Inside: '.threads->tid.' '.@{$self->{to_blast}}."\n";
    
    # Thread 'cancellation' signal handler
    $SIG{'KILL'} = sub { threads->exit(); };

	# Create the network
	while(
        (
            $self->{to_gbseq}->pending ||
            $self->{to_parse}->pending ||
            $self->{to_blast}->pending
        )# && 
        #$self->{FLAG}{COUNT} <= $self->{PARAM}{NUMITERATIONS})
	) {
		last if($counter >= $self->{PARAM}{NUMITERATIONS});
		
		## SEQUENCES
		threads->yield;
		if(my $tseq = $self->{to_gbseq}->dequeue_nb) {
			$work_to_do--;
            #print $tseq."\n"; die;
            
            eval{
				#print $self->generate_timestamp.'> sequences...'.threads->tid."\n";
				$tseq = batch2seq(
					$tseq,
					$self->{bad_seq},
					$self->{TARGET}{TYPE}
				);
			};
			if($@ || !$tseq) {
				carp('Failed to retrive sequence data!'."\n".$@."\n");
                
                # trying again may be extraneous...batch tries up to 5 times
                my $tstr = "$tseq";
                if(!exists $tried_seq->{$tstr}) {
                    $tried_seq->{$tstr} = 1;
                    $self->{to_gbseq}->enqueue($tseq);
                }
			}
			else {
                threads->yield;
                {
                    lock($self->{seq_data});
                    map { $self->{seq_data}{$_} = $tseq->{$_} } keys %{$tseq};
                }
                
                $self->{to_blast}->enqueue(keys %{$tseq});
                $work_to_do += scalar keys %{$tseq};
            }
		}
		
		## PARSE
		# $tpar is the file name to parse
		threads->yield;
		while(my $tpar = $self->{to_parse}->dequeue_nb) {
			$work_to_do--;
            my $tph;
            
			# Parse the file
			eval {
				#print $self->generate_timestamp.'> parsing...'.threads->tid."\n";
				$tph = $self->blastParse($tpar);                
			};
			if ($@) {
				carp('Failed to parse '.$tpar."\n".$@."\n");
			}
			
			map {
				delete($tph->{$_}) if(exists $self->{bad_seq}{$_});
				#delete($tph->{$_}) if(exists $self->{HASH}{$_});
			} keys %{$tph};
			
			## Pull ID from filename
			# <id>_<evalue>.blast; NP_11111_0.0001.blast
			$tpar =~ /^(?:BLAST\/)?(.*)\_\d*?\.\d*?\.blast$/;
			my $tpar_id = $1;
			
			# Save the parse results
			threads->yield;
			{
				lock($self->{HASH});
				if(!exists $self->{HASH}{$tpar_id}) {
					$self->{HASH}{$tpar_id} = {};
				}
				else {
					# we shouldn't be here...hopefully
					# this is adding the results of a BLAST query
					# should only be done when we don't have
					# the info...
				}
				$self->{HASH}{$tpar_id} = shared_clone($tph);
				$self->{to_gbseq}->enqueue([keys %{$tph}]);
                $work_to_do++;
			}
		}
		
		## BLAST
		threads->yield;
		if(my $tb = $self->{to_blast}->dequeue_nb) {
            $work_to_do--;
			
			# BLAST the sequence
			my $bf;
			eval{
				#print $self->generate_timestamp.'> blasting...'.threads->tid."\n";
				{lock($counter); $counter++;}
				print $counter.' :: '.threads->tid.' :: '."\n";
				$bf = $self->sablast(\$self->{seq_data}{$tb},\$tb);
				chmod(0777,$bf);					# CHANGE permissions
			};
			# Add results to be parsed
			if($@ || !-f $bf || !-s $bf) {
				carp('No BLAST data returned.'.$tb."\n".$@."\n");
				lock($self->{bad_seq});
				$self->{bad_seq}{$tb} = 1;
				#print $bf.'...bad'."\n";
				#{ lock($counter); $counter--; }
			}
			#elsif(!-s $bf) {
			#	unlink $bf;
			#	lock($self->{bad_seq});
			#	$self->{bad_seq}{$tb} = 1;
			#}
			else {
				threads->yield;
				lock($self->{seq_data});
				$self->{to_parse}->enqueue($bf);
				delete($self->{seq_data}{$tb});
                $work_to_do++;
			}
		}
	}
	
    #print 'Leaving '.threads->tid.' :: '.
    #    $self->{to_gbseq}->pending.' '.
    #    $self->{to_parse}->pending.' '.
    #    $self->{to_blast}->pending.' '.
    #    $work_to_do."\n";
	return 1;
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
	
	##### Write data
    my $temp_out = '#MAX '.$max."\n";											# PRINT the max value
    $temp_out .= '#MIN '.$min."\n";											# PRINT the min value
    $temp_out .= '#'.$self->{TARGET}{ID}."\t".keys(%{$self->{HASH}{$self->{TARGET}{ID}}})."\n";	# print seed acc no.
    foreach my $y (keys(%{$self->{HASH}{$self->{TARGET}{ID}}})) {			# print the seed results first
        $temp_out .= sprintf('%s: %.5f'."\n",$y,$self->{HASH}{$self->{TARGET}{ID}}{$y});
    }
    foreach my $g (keys(%{$self->{HASH}})) {							# rinse and repeat for each other sequence
        if($g ne $self->{TARGET}{ID}) {
            $temp_out .= '#'.$g."\t".scalar(keys(%{$self->{HASH}{$g}}))."\n";
            foreach my $h (keys(%{$self->{HASH}{$g}})) {
                $temp_out .= sprintf('%s: %.5f'."\n",$h,$self->{HASH}{$g}{$h});
            }
        }
    }
    open(my $lfs,'>'.$self->{DIR}{DATA}.$self->{FILES}{H_LODS});				# OPEN output file
    print $lfs $temp_out;
    close($lfs);																# CLOSE output file

	$self->{HASH}{MAX} = $max;												# ADD max value to tree
	$self->{HASH}{MIN} = $min;												# ADD min value to tree
	
	return $self->{HASH};														# RETURN ref to our tree
}

## 
sub create_matrix() {
	my $self 	= shift;														# GET our class data
	my @header	:shared = (); # shouldn't matter...right?															# INITIATE array
	my $grid	= { HEADER => '', DATA => '' };									# INITIATE matrix output string
	
	my $max = $self->{HASH}{MAX} * 1.5;										# GET max value and MODIFY to emmulate perfect score
	my $min = $self->{HASH}{MIN} * 5;											# GET min value and MODIFY to emmulate no score
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
		$self->{MATRIX}[$i] = &share([]);
        for(my $j = 0; $j <= $#header; $j++) {									# LOOP through 2nd tier
			if($i eq $j) {														# IF we're at a diagonal value (mirror)
				$grid->{DATA} .= $maxs;											#	-- ADD the max to output string data
                $self->{MATRIX}[$i][$j] = $max;
            } elsif(exists($self->{HASH}{$header[$i]}{$header[$j]})) {		# IF we're at a valid data point
				$grid->{DATA} .= sprintf('%.5f',$self->{HASH}{$header[$i]}{$header[$j]});	# 	-- ADD value to the output string data
                $self->{MATRIX}[$i][$j] = $self->{HASH}{$header[$i]}{$header[$j]};
            } else {															# ELSE we have no value
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
		carp(threads->tid.' :: Failed call to R. '.$@);
	}
	
	# Get R script results
	$self->{MATRIX} = $self->_input_matrix($self->{DIR}{DATA}.$self->{FILES}{H_PCORR});
	
	#unshift(@{$self->{MATRIX}},$header);
	
	print threads->tid.' '.$self->{MATRIX}[0][0]."\n";
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
	$grid = join("\t",@{$header})."\n";
	for(my $i = 0; $i <= $#{$self->{MATRIX}}; $i++) {
		for(my $j = 0; $j <= $#{$self->{MATRIX}[$i]}; $j++) {
			if($self->{MATRIX}[$i][$j] >= $lower &&
               $self->{MATRIX}[$i][$j] < $upper) {
				$self->{MATRIX}[$i][$j] = 1;
			}
			else {
				$self->{MATRIX}[$i][$j] = 0;
			}
		}
		$grid .= join("\t",@{$self->{MATRIX}[$i]})."\n";
	}
	unshift(@{$self->{MATRIX}},$header);
	chomp($grid);

	
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
	my $blast_file = $self->{DIR}{BLAST}.$$acc.'_'.$self->{PARAM}{EVALUE}.'.blast';
    
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

## BLAST PARSE
sub blastParse($) {
    my $self = shift;															# GET our class data
    my $infile = shift;														# GET filename
    my $scores = {};															# INITIATE an array ref for scores
    my $count = 0;																# INITIATE a count for total no. of scores
    my $acc = 0; my $bs = 0; my $ev = 0; # INITIATE important BLAST variables
    my $prefix = $self->{TARGET}{TYPE} eq 'p' ? 'NP' : 'NT';
    
    ##### read in the BLAST file data
    open(BF,$infile);															# OPEN the blast file
    my @blast_file = <BF>;														# READ in entire file
    close(BF);
    # CLOSE the file handle
    splice(@blast_file,0,20); # DELETE the unnecessary file data
    
    my $np = qr/^ref\|((?:$prefix)_\d{3,})(?:[\w\W[:punct:]]*)*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $gb = qr/^(?:\w+\|)?(?:\w+\d*_\w+)?([a-zA-Z]+\d{3,})\s*(?:[\w\W[:punct:]]*)*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $sp = qr/^[\w\d\_]*\s\(([a-zA-Z]+[\w\d]*)\).*\s+([\d\.]+)\s+(\d*(e-|\.)\d*)/;
	my $junk = qr/^.*\(?(\w+\d{3,})\)?.*([\d\.]+)\s+(\d+(e-|\.)\d*)/;
	
	##### parse through the data
    foreach my $line (@blast_file) { # LOOP through the file
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
			#print 'a....'.$1.':'.$2.':'.$3."\n";
        } elsif(!$self->{PARAMS}{REFSEQ_ONLY}) {
            if($line =~ /$gb/) {
                $acc = $1;															# SAVE accession number
                $bs = $2;															# SAVE raw score
                $ev = $3;															# SAVE e-value	
                #print 'b....'.$1.':'.$2.':'.$3."\n";
            } elsif($line =~ /$sp/) {
                $acc = $1;															# SAVE accession number
                $bs = $2;															# SAVE raw score
                $ev = $3;															# SAVE evalue
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


