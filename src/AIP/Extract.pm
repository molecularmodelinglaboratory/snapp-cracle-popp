#!/usr/bin/perl -w

=head aip.pm info

# Filename:	Extract.pm
# Author:	Stephen Bush
# Date:		01.28.07 [v0.2]
# Modified: 03.03.09 [v0.4]

=cut

package AIP::Extract;
require Exporter;

our @ISA      = qw( Exporter AIP::Core );
our @EXPORT   = qw(  ); # insert sub names here
use vars qw( );
our $version  = 1.00;

#use diagnostics;
use strict;
use warnings;
use Carp;

use AIP::Core;
use File::Copy;
use File::Path;
use Switch '__';

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
	#my $self = shift;
	
    my $sibling = shift;
	
    my $self = AIP::Core->new();
    $self->{DIR} = $sibling->{DIR};
    $self->{FILES} = $sibling->{FILES};
    $self->{PARAMS} = $sibling->{PARAMS};
    $self->{TARGET} = $sibling->{TARGET};
    $self->{DATE} = $sibling->{DATE};
    $self->{HANDLE} = $sibling->{HANDLE};
    $self->{FLAGS} = $sibling->{FLAGS};
	$self->{PARAM} = {
		THRESHOLD => 0,
		UPPER_THRESHOLD => 0
	};
	
	#$self->{FLAGS}{K}				= 0;
	#$self->{FLAGS}{THRESHOLD}		= 0;
	#$self->{FLAGS}{CLIQUES_DONE}	= 0;
	#$self->{FLAGS}{CORRELATION_II}	= 0;
	#$self->{FLAGS}{FIND}			= 0;
	#$self->{FLAGS}{GRAPH}			= 0;
	#$self->{FLAGS}{MAXIMAL}			= 0;
	#$self->{FLAGS}{NETWORK}			= 0;
	#$self->{FLAGS}{SUBMAXIMAL}		= 0;
	#$self->{FLAGS}{TRUTHTABLE}		= 0;
	#$self->{COUNT}					= {};
	#$self->{REVERSE_HEADER}			= {};
	#$self->{HEADER}				= ();
	#$self->{CORRELATION_I}			= undef;
	#$self->{CORRELATION_II}			= undef;
	$self->{HASH} = {};
    
	bless($self,$class);
	return $self;
}

# DESTROY
sub DESTROY {

}

# set_original
sub set_original {
	my $self = shift;
	if (@_) {
		$self->{ORIGINAL} = shift;
		$self->{ORG_HASH} = $self->matrix_to_hash($self->{ORIGINAL},0.50);
	}
	return $self->{ORG_HASH};
}

## _EDGES
sub _edges {
	my $self = shift;
	my $count = 0;
	
	undef($self->{TREE});
	for(my $a = 0; $a <= $#{$self->{HEADER}}; $a++) {
		$count = 0;
		$self->{TREE}{$self->{HEADER}[$a]} = ();
		$self->{REVERSE_HEADER}{$self->{HEADER}[$a]} = $a;
		for(my $z = 0; $z <= $#{$self->{HEADER}}; $z++) {
			if($self->{GRID}[$a][$z] >= $self->{FLAGS}{THRESHOLD} && $a != $z) {
				$self->{TREE}{$self->{HEADER}[$a]}{$self->{HEADER}[$z]} = $self->{GRID}[$a][$z];
				$count++;
			}
		}
		if($count == 0) { delete $self->{TREE}{$self->{HEADER}[$a]}; }
	}
	
	return $self->{TREE};
}

## _EDGES_DIMACS
sub _edges_dimacs {
	my $self = shift;
	my $outfile = shift;
	my $count = 0;
	my $edges;

	##### create DIMACS output format
	for(my $a = 0; $a <= $#{$self->{HEADER}}; $a++) {
		$edges .= "c <$a> <$self->{HEADER}->[$a]>\n";
		for(my $z = $a; $z <= $#{$self->{HEADER}}; $z++) {
			if($self->{GRID}[$a][$z] eq 'NA') { 
				die("\n$a\t$z\t$@\n\n"); }
			if($self->{GRID}[$a][$z] >= $self->{FLAGS}{THRESHOLD} && $a != $z) {
				$a++; $z++;
				$edges .= "e $a $z\n";
				$count++;
				$a--; $z--;
			}
		}
	} 
	my $temp = $#{$self->{HEADER}} + 1;
	$edges = "p edge $temp $count\n".$edges;
	
	##### open output file
	open(DIMACS,">$outfile");
	print DIMACS $edges;
	close(DIMACS);
	
	return $outfile;
}

## CLIQUE_MAXIMAL
sub clique_maximal {
	my $self = shift;
	my $graph_files;
	my $temp_k = $self->{FLAGS}{K}; $self->{FLAGS}{K} = 0;
	my $temp_submaximal = $self->{FLAGS}{SUBMAXIMAL}; $self->{FLAGS}{SUBMAXIMAL} = 0;
	my $temp_file = "$self->{DIR}{THRESHOLD}$self->{FILES}{CLIQUE_MAXIMAL}";
	my $total_graph_count;
	
	print "Generating maximal cliques...\n";
	$temp_file = $self->_edges_dimacs($temp_file);
	push(@{$graph_files},$self->_calculate_clique($self->{HEADER},$temp_file));	#
	$total_graph_count += $#{$graph_files->[$#{$graph_files}]} + 1;
	
	$self->{FLAGS}{CLIQUES_DONE} = 1;											# SET flag for cliques done
	
	##### generating graphviz
	my $count = 1;
	if($self->{FLAGS}{FIND}) {
		if(exists($self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}})) {
			my $seed_cliques = $self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}};
			$self->{TABLE_OF_CONTENTS} = {};
			$self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}} = $seed_cliques;
			$self->{DIR}{FIND} = "$self->{DIR}{THRESHOLD}find$self->{FLAGS}{FIND}/";
			unless(-d $self->{DIR}{FIND}) { mkdir $self->{DIR}{FIND}, 0777 or croak(@_); }
			if($self->{FLAGS}{GRAPH}) {
				print "Generating maximal clique graphs with $self->{FLAGS}{FIND}.\n";
				$total_graph_count = $#{$seed_cliques} + 1;
				for(my $a = 0; $a <= $#{$seed_cliques}; $a++) {
					print "\rGenerating maximal clique graph [$count] of [$total_graph_count]...   ";
					$self->_publish_graphviz($seed_cliques->[$a],$self->{DIR}{FIND});
					$count++;
				} print "\n";
			}
		} else { print "No maximal cliques containing $self->{FLAGS}{FIND} were found.\n"; }
	} elsif($self->{FLAGS}{GRAPH}) {
		for(my $a = 0; $a <= $#{$graph_files}; $a++) {
			for(my $z = 0; $z <= $#{$graph_files->[$a]}; $z++) {
				print "\rGenerating maximal clique graph [$count] of [$total_graph_count]...   ";
				$self->_publish_graphviz($graph_files->[$a][$z]);
				$count++;
			}
		} print "\n";
	}
	
	$self->{FLAGS}{K} = $temp_k;
	$self->{FLAGS}{SUBMAXIMAL} = $temp_submaximal;
	return 0;
}

## _CLIQUE_SUBMAXIMAL
sub clique_submaximal {
	my $self 		= shift;
	my $temp_k 		= $self->{FLAGS}{K}; $self->{FLAGS}{K} = 0;
	my $edges 		= {
		MAIN		=> '',
		SECONDARY	=> '' };
	my $node		= {};
	my @cipher;
	my $num_edges	= 0;
	my $num_nodes	= 0;
	my $tier		= 0;
	my $graph_files;
	my $temp_dir;
	my $outfile;
	my $total_graph_count = 0;
	
	my $header_count = $#{$self->{HEADER}};
	for(my $i = 0; $i <= $header_count; $i ++) {								# LOOP through all of our acc
		my $seed = $self->{HEADER}[$i];									# - SAVE our acc to an easier variable
		print "\rCalculating subgraph [$i] of [$header_count] $seed...    ";		# - NOTIFY the user
		if(!exists($node->{$seed})) {									# - IF we don't have run through this acc before
			$num_nodes++;												# -- INCREMENT the number of nodes
			push(@cipher,$seed);										# -- PUSH 
			$node->{$seed} = $num_nodes;									# -- ADD it to the cipher haip
		}
		$edges->{MAIN} .= "c [$seed] \n";									# - ADD comment to DIMACS file -- first node
		foreach my $set (keys %{$self->{CORRELATION_II}{$seed}}) {				# - LOOP through SEED DEPTH
			if(!exists($node->{$set})) {									# -- IF we don't have run through this acc before
				$num_nodes++;											# --- INCREMENT the number of nodes
				push(@cipher,$set);										# --- PUSH 
				$node->{$set} = $num_nodes;								# --- ADD it to the cipher haip
			}
			$edges->{MAIN} .= "e $node->{$seed} $node->{$set}\n";				# -- ADD edge (numeric node) to DIMACS file
			$num_edges++;												# -- INCREMENT the number of edges
			$edges->{SECONDARY} .= "c [$set]\n";							# -- ADD comment to DIMACS file -- next node
			foreach my $subset (keys %{$self->{CORRELATION_II}{$set}}) {		# -- LOOP through SEED DEPTH + 1
				if($subset eq $seed) { next; }							# --- IF we already have this pair 
				if(!exists($node->{$subset})) {							# --- IF we don't have run through this acc before
					$num_nodes++;										# ---- INCREMENT the number of nodes
					push(@cipher,$subset);								# ---- PUSH 
					$node->{$subset} = $num_nodes; }						# ---- ADD it to the cipher haip
				$edges->{SECONDARY} .= "e $node->{$set} $node->{$subset}\n";	# --- ADD edge to DIMACS file
				$num_edges++;											# --- INCREMENT number of edges
			}
		}
		$edges->{MAIN} = "p edge $num_nodes $num_edges\n".$edges->{MAIN};		# - ADD file heading to DIMACS file
		$self->{DIR}{TEMP} = sprintf("%sTier_%04d/",						# - FORMULATE a directory name
			$self->{DIR}{THRESHOLD},
			$tier);
		mkdir "$self->{DIR}{TEMP}", 0777;									# - MAKE the directory
		$outfile = $self->{DIR}{TEMP}.$seed.".edges";						# - SAVE the filename
		
		open(DIMACS,">$outfile");										# - OPEN DIMACS file
		print DIMACS $edges->{MAIN}.$edges->{SECONDARY};						# - PRINT edges to DIMACS file
		close(DIMACS);													# - CLOSE DIMACS file
		
		push(@{$graph_files},$self->_calculate_clique(\@cipher,$outfile));		# - CALCULATE the cliques and save the filenames
		$total_graph_count += $#{$graph_files->[$#{$graph_files}]} + 1;
		
		$tier++;														# - INCREMENT the tier number
		$num_nodes 	= 0;												# - RESET the number of nodes
		$num_edges 	= 0;												# - RESET the number of edges
		undef $node;													# - RESET the cipher
		undef @cipher;													# - RESET the decipher
		$edges 		= {
			MAIN		=> '',
			SECONDARY	=> ''  };
	} print "\n";
	
	$self->{FLAGS}{CLIQUES_DONE} = 1;										# SET flag for cliques done
	
	##### publish graphviz
	my $count = 1;
	if($self->{FLAGS}{FIND}) {
		if(exists($self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}})) {
			my $seed_cliques = $self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}};
			$self->{TABLE_OF_CONTENTS} = {};
			$self->{TABLE_OF_CONTENTS}{$self->{FLAGS}{FIND}} = $seed_cliques;
			$self->{DIR}{FIND} = "$self->{DIR}{THRESHOLD}find$self->{FLAGS}{FIND}/";
			unless(-d $self->{DIR}{FIND}) { mkdir $self->{DIR}{FIND}, 0777 or croak(@_); }
			if($self->{FLAGS}{GRAPH}) {
				print "Generating submaximal clique graphs with $self->{FLAGS}{FIND}.\n";
				$total_graph_count = $#{$seed_cliques} + 1;
				for(my $a = 0; $a <= $#{$seed_cliques}; $a++) {
					print "\rGenerating submaximal clique graph [$count] of [$total_graph_count]...   ";
					$self->_publish_graphviz($seed_cliques->[$a],$self->{DIR}{FIND});
					$count++;
				} print "\n";
			}
		} else { print "No submaximal cliques containing $self->{FLAGS}{FIND} were found.\n"; }
	} elsif($self->{FLAGS}{GRAPH}) {
		for(my $a = 0; $a <= $#{$graph_files}; $a++) {
			for(my $z = 0; $z <= $#{$graph_files->[$a]}; $z++) {
				print "\rGenerating submaximal clique graph [$count] of [$total_graph_count]...   ";
				$self->_publish_graphviz($graph_files->[$a][$z]);
				$count++;
			}
		} print "\n";
	}
	
	
	$self->{FLAGS}{K} = $temp_k;
	return 0;
}

####
# MATRIX_TO_HASH
sub matrix_to_hash {
    my $self = shift;
    my $grid = shift;
    my $threshold = shift;
    my $head = shift(@{$grid});
	my $hash = {};
    
    # j = i + 1 serves to eliminate the diagonal
    # grid should be square, i.e. length works for both checkpoints
    my $grid_length = $#{$grid};
    for(my $i = 0; $i <= $grid_length; $i++) {
		#print $i."\t";
        for(my $j = $i + 1; $j <= $grid_length; $j++) {
            #print $j."\t".$grid->[$i][$j]."\n"; die;
			next if($grid->[$i][$j] < $threshold);
            $hash->{$head->[$i]} = {} if(!$hash->{$head->[$i]});
            $hash->{$head->[$i]}{$head->[$j]} = $grid->[$i][$j];
            $hash->{$head->[$j]}{$head->[$i]} = $grid->[$i][$j];
        }    
    }
    unshift(@{$grid},$head);
	
	$self->{HASH} = $hash;
	
    return $hash;    
}

####
# FIND_NETWORK
sub find_subnet {
	my $self = shift;
    my $find = @_ ? shift : $self->{TARGET}{ID};
	
	my $subnet = { $find => {} };
	
    # Depth of 2
    
    # I wonder if we might be missing a few here...
    # but since we're saving
    foreach my $i (keys %{$self->{HASH}{$find}}) {
        #print 'a:'.$i."\n";
        # We are looping through every node $i connected to $find
        # Then we loop through every node $j connected to $i
        # Each time we save the edge connection $find->$i
        #   or $i->$j
        # Skipping $find->$i is extraneous, in $i->$j we wouldn't
        #   make $find->$j
        # UNLESS $find == $i, but j = i +1 in matrix_to_hash covers this
        # Skipping $i->$j makes sense since $j will = $find at some point
        # Do we really need to duplicate?
        #   If we're looping through the entire hash every time...no
        #   Better safe than sorry   
        #next if(exists $subnet->{$find}{$i}); # reasoning: if we have either, we've iterated through one of em
        $subnet->{$i} = {} if(!exists $subnet->{$i}); # but if we don't, there might still be something in $i
        $subnet->{$find}{$i} = $self->{ORG_HASH}{$find}{$i};
        $subnet->{$i}{$find} = $self->{ORG_HASH}{$find}{$i};
        
        foreach my $j (keys %{$self->{HASH}{$i}}) {
            next if(exists $subnet->{$i}{$j});
            $subnet->{$j} = {} if(!exists $subnet->{$j});
            $subnet->{$i}{$j} = $self->{ORG_HASH}{$i}{$j};
            $subnet->{$j}{$i} = $self->{ORG_HASH}{$i}{$j};
			#print $i.'::'.$j.'::'.$subnet->{$j}{$i}."\n";
        }
    }
    
    # Because of looping through $i,
    # If $j is not in $find, our target may not connect to all nodes in the subnet
    # However, our count is based upon the number of edges from $find
	#$subnet->{ctot}	= scalar keys %{$subnet->{$find}}; 
	#$subnet->{c25p} = $subnet->{ctot} * 0.25;
	#$subnet->{c50p} = $subnet->{ctot} * 0.50;
	#$subnet->{c75p} = $subnet->{ctot} * 0.75;
    
    return $subnet;
}

## GRAPHVIZ SUBNETWORKS
sub graphviz_subnet {
	my $self = shift;
    my $subnet = shift;
    my $label = shift;
	my $find = @_ ? shift : $self->{TARGET}{ID};
    
	##### ARRANGE by seed highest correlation
	my @order = sort {
        $subnet->{$find}{$b} <=> $subnet->{$find}{$a} } 
		keys %{$subnet->{$find}};
	unshift(@order,$find);
	
	##### PREPARE a Table of Contents and Graphviz file
	my $graph = '';												# INITIALIZE graph
	my $subnet_omit = {};											# RESET omit
	my $subnet_count = {};
    
	my $count_total = scalar keys %{$subnet->{$find}};
	my $count_quarter = $count_total * 0.25;							# CALCULATE one fourth of number
	my $count_half = $count_total * 0.50;						# CALCULATE half
	my $count_threequarter	= $count_total * 0.75;						# CALCULATE three quarters
	
	foreach my $z (@order) {									# LOOP through our ordered list
		$subnet_omit->{$z} = 1;								# - ADD node to omit list
		
		$subnet_count->{$z} = scalar keys %{$subnet->{$z}};
        
        # switches with comparisons are...sticky
		my $node_type;
		if($subnet_count->{$z} < $count_quarter) {
			$node_type = '{ node [style=filled,color="#ffffff",fontcolor="#000000"] '.$z.' }';
        }
		elsif($subnet_count->{$z} < $count_half) {
			$node_type = '{ node [style=filled,color="#99ffff",fontcolor="#000000"] '.$z.' }';
        }
		elsif($subnet_count->{$z} < $count_threequarter) {
			$node_type = '{ node [style=filled,color="#99ccff",fontcolor="#000000"] '.$z.' }';
        }
		elsif($subnet_count->{$z} < $count_total) {
			$node_type = '{ node [style=filled,color="#6666ff",fontcolor="#FFFFFF"] '.$z.' }';
        }
		elsif($subnet_count->{$z} == $count_total) {
			if($z eq $find) {
				$node_type = '{ node [shape="circle",style=filled,color="#000099",fontcolor="#FFFFFF"] '.$z.' }';
            }
			else {
				$node_type = '{ node [style=filled,color="#000099",fontcolor="#FFFFFF"] '.$z.' }';
            }
        }
		elsif($subnet_count->{$z} > $count_total) {
            # This is possible, just highly unlikely
            # Since we find all nodes connected to all nodes from our target
            #   we may end up with nodes which are better hubs than our target
            $node_type = '{ node [style=filled,color="#990099",fontcolor="#FFFFFF"] '.$z.' }';
        }
        else{}
		$graph = $node_type."\n".$graph;
		
		foreach my $y (keys %{$subnet->{$z}}) {				# - LOOP through the node's list
			
            # Checking for $y in omit is merely to make sure we don't
            #   log the same edge twice
            #   i.e. we've already looped through $y as a $z
            if(!exists($subnet_omit->{$y})) {					# -- IF we haven't already looped through this node
				$graph .= $z.' -- '.$y.' ';						# --- ADD the node to the graph
				$graph .= $self->_graphviz_edge($subnet->{$z}{$y});
			}
		}
	}
    
    return $self->_generate_graphviz(\$graph,$label,$find);
}


## NETWORK
sub graphviz_network {
	my $self = shift;
	my $thresh = @_ ? shift : 0.50;
	my $use_org = @_ ? shift : 0;
	my $file = $self->{DIR}{GRAPH}.'aip_network.graphviz';
	my $graph = #'graph { '."\n".
		#'label="'.$file.'";'."\n".
		#'overlap=false; splines=true;'."\n".
		'{ node [shape="circle",style=filled,color="#000099",fontcolor="#FFFFFF"] '.$self->{TARGET}{ID}.' }'."\n";
	my $non_priority = '';
	
	my $hash = $use_org ? $self->{ORG_HASH} : $self->{HASH};
	
	print 'Creating network file...'."\n";
	foreach my $i (keys %{$hash}) {
		foreach my $j (keys %{$hash->{$i}}) {
			next if($hash->{$i}{$j} < $thresh);
			if($hash->{$i}{$j} < 0.70) {
				$non_priority .= $i.' -- '.$j;
				$non_priority .= $self->_graphviz_edge($hash->{$i}{$j});
			}
			else {
				$graph .= $i.' -- '.$j;
				$graph .= $self->_graphviz_edge($hash->{$i}{$j});
			}
		}
	} #$graph .= '}';

	$thresh = $thresh * 100;
	$graph .= $non_priority;

	#open(GRAPH,">$file");
	#print GRAPH $graph;
	#close(GRAPH);
	
	return $self->_generate_graphviz(\$graph,$use_org.'n'.$thresh,$self->{TARGET}{ID});
	
}

sub _graphviz_edge {
	my $self = shift;
	my $t = shift;
	my $l = '';
	
	if($t >= 0.95) {
		$l = '[color="#FF0000"];'."\n"; }
	elsif($t >= .90) {
		$l = '[color="#FFCC00"];'."\n"; }
	elsif($t >= .85) {
		$l = '[color="#00CC00"];'."\n"; }
	elsif($t >= .80) {
		$l = '[color="#00CCCC"];'."\n"; }
	elsif($t >= .75) {
		$l = '[color="#0000CC"];'."\n"; }
	elsif($t >= .70) {
		$l = '[color="#CC00CC"];'."\n"; }
	elsif($t >= .50) {
		$l = '[color="#CCCCCC"];'."\n"; }
	else{}
	
	return $l;
}

## CALCULATE CLIQUE
sub _calculate_clique() {
	my $self 	= shift;
	my $cipher	= shift || $self->{HEADER};
	my $infile	= shift || $self->{FILES}{CLIQUE_II};
	my $outfile	= $infile;	
	my $edges	= '';
	my $count	= 0;
	my $files 	= ();
	my @clique_list;
	my @nodes;
	my @line_split;
	my $proper_index;
	my %current_cliques;
	my $gv_file;
	
	##### dealing with the outfile
	$files->[0] = $outfile;
	$outfile =~ s/(.*)\.\w*/$1/;
	$files->[1] = $outfile.".clique";
	$files->[2] = ();

	##### run cliquer
	my @args;
	if($self->{FLAGS}{ALL_CLIQUES}) { 
		@args = ("~/lib/cliquer-1.2/cl -q -q -a -m $self->{FLAGS}{K} $infile > $files->[1]"); }
	else { 
		@args = ("~/lib/cliquer-1.2/cl -q -q -a $infile > $files->[1]"); }
	system(@args) == 0
		or croak("system @args failed: $?");
		
	#### get cliquer output
	open(CLIQUER,"<$files->[1]");										# OPEN list of clique nodes
	my @cliquer = <CLIQUER>;											# READ in the file
	close(CLIQUER);												# CLOSE the file
	foreach my $line (@cliquer) {									# LOOP through each clique returned
		chomp($line);												# REMOVE newlines
		@line_split = split(/:\s*/,$line);								# SPLIT clique data and clique
		$current_cliques{$line_split[1]} = 1;							# ADD clique to those we have
		$line_split[0] =~ s/size=(\d*).*/$1/;							# GET the clique size (k)
		if($count == 0 && $self->{FLAGS}{SUBMAXIMAL}) {					# IF first clique from file and submaximal
			if($line_split[0] <= 3) { 								# IF k is 3 or less
				rmtree($self->{DIR}{TEMP}); 							#	-- delete associated files
				##pring WARNING "$outfile returned a clique of > 3.\n";		# 	-- WARNING
				last; }											#	-- SKIP ALL (all are maximal)
			$self->{DIR}{CLIQUE} = sprintf("%sk%03d/",					# SAVE new directory name based on clique size
				$self->{DIR}{THRESHOLD},$line_split[0]);
			if(!-d $self->{DIR}{CLIQUE}) { mkdir $self->{DIR}{CLIQUE}, 0777; }	# IF directory doesn't exist, create it
			$files->[0] =~ s|$self->{DIR}{TEMP}||;						# REMOVE the dir from variables
			$files->[1] =~ s|$self->{DIR}{TEMP}||;						# REMOVE the dir from variables
			$outfile	=~ s|$self->{DIR}{TEMP}|$self->{DIR}{CLIQUE}|;		# CHANGE the dir for variable
			move("$self->{DIR}{TEMP}$files->[0]","$self->{DIR}{CLIQUE}$files->[0]");	# COPY from temp dir to clique dir
			move("$self->{DIR}{TEMP}$files->[1]","$self->{DIR}{CLIQUE}$files->[1]");	# COPY from temp dir to clique dir
			rmtree($self->{DIR}{TEMP});	}							# DELETE the temp dir
		@nodes = split(/\s/,$line_split[1]);							# SPLIT the clique string into an array
		for(my $i = 0; $i <= $#nodes; $i++) {							# LOOP through nodes
			$proper_index = $nodes[$i] - 1;							# DECIPHER clique node into our node list
			$nodes[$i] = $cipher->[$proper_index]; }					# TRANSLATE node from int to acc
		$gv_file = $self->_generate_graphviz(\@nodes,$outfile,$count);		# CREATE a graphviz file from clique data
		push(@{$files->[2]},$gv_file);								# ADD file to our list
		for(my $i = 0; $i <= $#nodes; $i++) {							# LOOP through nodes
			push(@{$self->{TABLE_OF_CONTENTS}{$nodes[$i]}},$gv_file); }		# ADD the file to table of contents
		$count++;													# INCREMENT count of cliques per file
		undef @line_split;											# DELETE our line
		undef @nodes;												# DELETE our nodes
	}
	
	return $files->[2];											# RETURN our list of graphviz files
}

## GENERATE GRAPHVIZ
sub _graphviz_clique {
	my $self	= shift;
	my $clique	= shift or croak("No clique input.\n");
	my $file	= shift;
	my $page	= shift || 0;
	my $graph	= '';
	$file =~ s/(.*)\..*/$1/;
	if($page)	{ $file .= sprintf("_%02d.graphviz",$page); }
	else		{ $file .= '.graphviz'; }
	
	$graph = "graph {
		label=\"$file\";
		overlap=false; splines=true;
		{ node [style=filled,color=\"#000099\",fontcolor=\"#FFFFFF\"] $self->{TARGET}{ID}}\n";
	for(my $a = 0; $a <= $#{$clique}; $a++) {
		for(my $z = $a + 1; $z <= $#{$clique}; $z++) {
			$graph .= "$clique->[$a] -- $clique->[$z] ";
			switch($self->{CORRELATION_I}{$clique->[$a]}{$clique->[$z]}) {
				case __ < .70 { $graph .= "[color=\"#CCCCCC\"];\t"; }	# Light Gray
				case __ < .75 { $graph .= "[color=\"#0099FF\"];\t"; }	# Light Blue
				case __ < .80 { $graph .= "[color=\"#00CCCC\"];\t"; }	# Cyan
				case __ < .85 { $graph .= "[color=\"#00CC00\"];\t"; }	# Green
				case __ < .90 { $graph .= "[color=\"#FFCC00\"];\t"; }	# Orange
				case __ < .95 { $graph .= "[color=\"#FF0000\"];\t"; }	# Red
				case __ < 1.0 { $graph .= "[color=\"#000000\"];\t"; }	# Black
				else{}
			} $graph .= "\n";
		}
	} $graph .= "}";
	
	open(GRAPH,">$file");
	print GRAPH $graph;
	close(GRAPH);
	
	return $file;	
}

## GENERATE GRAPHVIZ
sub _generate_graphviz {
    my $self = shift;
    my $graph_data = shift;
    my $label = shift;
	my $find = @_ ? shift : $self->{TARGET}{ID};
    
    my $graph_file = $self->{DIR}{GRAPH}.$find.'_'.$label;
    
    my $graph_final = 'graph { '."\n".
		"\t".'label="'.$graph_file.'";'."\n".
		"\t".'overlap=false; splines=true;'."\n".
		$$graph_data.
        '}';
    
    open(my $gfh,'>',$graph_file.'.data');
	print $gfh $graph_final;
	close($gfh);
    
    my $gcmd = ['neato -Tpng -o'.$graph_file.'.png '.$graph_file.'.data &'];
    return $gcmd;
}

## PUBLISH GRAPHVIZ
sub _publish_graphviz {
	my $self	= shift;
	my $infile	= shift || $self->{FILES}{GRAPH};
	my $dir		= shift || 0;
	
	if(!-f $infile) {
		$infile = "$self->{DIR}{THRESHOLD}$infile";
    }
	my $outfile = shift || $infile;
	$outfile =~ s/(.*)\.graphviz/$1\.png/;
	
	if(-d $dir) {
		$outfile =~ s|.*(\/\w*\.png)|$dir$1|; }

	##### call graphviz
	my @args = ("neato","-Tpng","-o$outfile","$infile");
	system(@args) == 0
	 or croak "system @args failed: $?";
	 
	return "$self->{DIR}{THRESHOLD}$outfile";
}


## HEADER
sub _header {
	my $self = shift;
	$self->{HEADER} = shift(@{$self->{GRID}});	
	
	##### put our seed first
	for(my $acc = 0; $acc <= $#{$self->{HEADER}}; $acc++) {
		if($self->{HEADER}->[$acc] eq $self->{TARGET}{ID}) { splice @{$self->{HEADER}}, $acc, 1; last; }
	} unshift(@{$self->{HEADER}},$self->{TARGET}{ID});
}

## ????
#sub initialize_extract {
#	my $self = shift;
#	my $fh;
	#my $infile = @_ ? shift : $self->{FILES}{CORRELATION};
	
	#if(-f "$self->{DIR}{THRESHOLD}$infile") {
	#	$infile = "$self->{DIR}{THRESHOLD}$infile";
	#}
	
	#$self->SUPER::_read_grid($infile);
	#$self->_header();
#	
#For some reason, qr is acting like a subroutine and not acting as if
#the variables are references.
#	
#	my $row = 0;
#	my $col = 0;
#	
#	my $regex = qr/
#		(:?
#			(:?
#				([\w]+)\W?
#				(?{ $self->{MATRIX}[$row][$col++] = $1; })
#			)+
#			\n
#			(?{ $row++; $col = 0; })
#		)
#	/xi;
#	
#	open($fh,'<'.$self->{DIR}{DATA}.$self->{FILES}{PCORR1});
#	my @lines = <$fh>;
#	close($fh);
#	
#	my $dangerously_huge_slurp = join('',@lines);
#	@lines = undef;
#	
#	$dangerously_huge_slurp = s/$regex//g;
#
#}


#####
## FLAGS

sub maximal {
	my $self = shift;
	if (@_) { $self->{FLAGS}{MAXIMAL} = shift; }
	return $self->{FLAGS}{MAXIMAL};
}
sub submaximal {
	my $self = shift;
	if (@_) { $self->{FLAGS}{SUBMAXIMAL} = shift; }
	return $self->{FLAGS}{SUBMAXIMAL};
}
sub network {
	my $self = shift;
	if (@_) { $self->{FLAGS}{NETWORK} = shift; }
	return $self->{FLAGS}{NETWORK};
}
sub threshold {
	my $self = shift;
	if (@_) { $self->{FLAGS}{THRESHOLD} = shift; }
	return $self->{FLAGS}{THRESHOLD};
}
sub k {
	my $self = shift;
	if (@_) { $self->{FLAGS}{K} = shift; }
	return $self->{FLAGS}{K};
}
sub graphviz {
	my $self = shift;
	if (@_) { $self->{FLAGS}{GRAPH} = shift; }
	return $self->{FLAGS}{GRAPH};
}
sub graph_seed_only {
	my $self = shift;
	$self->{FLAGS}{FIND} = $self->{TARGET}{ID};
	return $self->{FLAGS}{FIND};
}

sub extractor_recovery {
	my $self = shift;
	
	##### check recovery status
	my $verified = $self->verify_project();
	if(!$verified) { croak("Invalid project.\n\n"); }
	elsif ($verified == 4) { print "No recovery needed.\n"; return 0; }
	$self->{FLAGS}{RECOVERY} = 1;
	
	my $collector = AIP::Collector->new();
	
	$collector->acc($self->{TARGET}{ID});
	$collector->directory($self->{DIR}{MAIN});
	if($collector->collector_recovery()) { croak("Unable to recover.\n\n"); }
	
	return 0;
}

sub find {
	my $self = shift;
	if(@_) { $self->{FLAGS}{FIND} = shift; }
	return $self->{FLAGS}{FIND};
}

1;
__END__