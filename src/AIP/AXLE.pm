#!/usr/bin/perl -w

# Filename:	AXLE.pm
# Author:	Stephen Bush <sjbush@unc.edu>
# Date:		01.28.07 [v0.2]
# Modified: 03.03.09 [v0.4]
# AIP XML Log Engine (AXLE)

package AIP::AXLE;
require Exporter;

our @ISA      = qw( Exporter );
our @EXPORT   = qw( ); # insert sub names here
use vars qw( );
our $version  = 0.1;

use strict;
use Carp;
#use AIP::Core;

####
# CONSTRUCTOR & DESTRUCTOR

## NEW
sub new {
    my $class = shift;
    my $self =  {
        LOG => {
            HEAD => '<?xml version="1.0" encoding="UTF-8" ?>'."\n".
                '<?xml-stylesheet type="text/xsl" href="aip_axle.xsl"?>'."\n".
		'<AIP>'."\n",
            FOOT => "\n".'</AIP>'."\n",
            PARAM => '',
            TARGET => '',
            DATE => '',
            SUBNET => ''
        },
        DIR => shift,
        QCOMMIT => []
    };
    push(@{$self->{QCOMMIT}},'HEAD');
    bless($self,$class);
    return $self;    
}

## DESTORY
sub destroy {
    my $self->shift;
    
    $self->finalize;
}

####
# LOG FUNC

# LOG_INIT
sub log_init {
    my $self = shift;
    my $a = shift;
    my $b = shift;
    my $c = shift;
   
    $self->log_param($b);
    $self->log_target($a);
    $self->log_date($c);

    $self->{FILES}{XML} = $self->{TARGET}{ID}.'.xml';
    
    if(-f $self->{DIR}{XML}.$self->{FILES}{XML}) {
        unlink($self->{DIR}{XML}.$self->{FILES}{XML});
    }

    return 1;
}

# LOG_TARGET
sub log_target {
    my $self = shift;
    $self->{TARGET} = shift;
    
    $self->{LOG}{TARGET} = "\n".
        '<target is_protein="'.$self->{PARAM}{IS_PROTEIN}.'">'."\n".
        "\t".'<id>'.$self->{TARGET}{ID}.'</id>'."\n".
        "\t".'<sequence>'.$self->{TARGET}{SEQ}.'</sequence>'."\n".
    '</target>';
    
    push(@{$self->{QCOMMIT}},'TARGET');
    return 1;
}

# LOG_PARAM
sub log_param {
    my $self = shift;
    $self->{PARAM} = shift;
    
    $self->{LOG}{PARAM} = "\n".'<parameters>'."\n".
        "\t".'<db>'.$self->{PARAM}{DB}.'</db>'."\n".
        "\t".'<program>'.$self->{PARAM}{PROGRAM}.'</program>'."\n".
        "\t".'<numiterations>'.$self->{PARAM}{NUMITERATIONS}.'</numiterations>'."\n".
        "\t".'<evalue>'.$self->{PARAM}{EVALUE}.'</evalue>'."\n".
        "\t".'<hypothetical>'.$self->{PARAM}{HYPOTHETICAL}.'</hypothetical>'."\n".
    '</parameters>';
    
    push(@{$self->{QCOMMIT}},'PARAM');
    return 1;
}

# LOG_DATE
sub log_date {
    my $self = shift;
    $self->{DATE} = shift;
    
    $self->{LOG}{DATE} = "\n".'<date>'."\n".
        "\t".'<year>'.$self->{DATE}{YEAR}.'</year>'."\n".
        "\t".'<month>'.$self->{DATE}{MONTH}.'</month>'."\n".
        "\t".'<day>'.$self->{DATE}{DAY}.'</day>'."\n".
        "\t".'<hour>'.$self->{DATE}{HOUR}.'</hour>'."\n".
        "\t".'<minute>'.$self->{DATE}{MINUTE}.'</minute>'."\n".
        "\t".'<second>'.$self->{DATE}{SECOND}.'</second>'."\n".
    '</date>';

    push(@{$self->{QCOMMIT}},'DATE');
    return 1;
}

# LOG_KEY
sub log_key {
    my $self = shift;
    my $key = shift;
    
    $self->{LOG}{KEY} = "\n".'<key>'."\n";
    
    foreach my $i (sort keys %{$key}) {
        $self->{LOG}{KEY} .= "\t".'<member id="'.$i.'">'.
            $key->{$i}.'</member>'."\n";
    }
    
    $self->{LOG}{KEY} .= "\n".'</key>'."\n";
    
    push(@{$self->{QCOMMIT}},'KEY');
    return 1;
}

# LOG_SUBNET
sub log_subnet {
    my $self = shift;
    my $subnet = shift;
    my $label = shift;
    my $lower = shift;
    my $upper = shift;
    my $find = @_ ? shift : $self->{TARGET}{ID};
    
    # Note: Redo this later...it reads through the subnetwork twice
    # ...not very effic, but it should be small, so not a huge deal
    # right now.
    # grab edge count
    my $count_100 = scalar keys %{$subnet->{$find}};
	my $count_75p = $count_100 * 0.75;
	my $count_50p = $count_100 * 0.50;
	my $count_25p = $count_100 * 0.25;
	
    my $nodes = {
        edges_g100 => {},
        edges_100 => {},
        edges_75p => {},
        edges_50p => {},
        edges_25p => {},
        edges_0p => {}
    };
    foreach my $i (keys %{$subnet}) {
        my $count = keys %{$subnet->{$i}};
        if($count > $count_100) {
            $nodes->{edges_g100}{$i} = $count;
        }
        elsif($count == $count_100) {
            $nodes->{edges_all}{$i} = $count;
        }
        elsif($count < $count_100) {
            $nodes->{edges_75p}{$i} = $count;
        }
        elsif($count < $count_75p) {
            $nodes->{edges_50p}{$i} = $count;
        }
        elsif($count < $count_50p) {
            $nodes->{edges_25p}{$i} = $count;
        }
        elsif($count < $count_25p) {
            $nodes->{edges_0p}{$i} = $count;
            
        }
    }
    
    $self->{LOG}{$label} = "\n".'<subnet name="'.$label.'" '.
        'lower="'.$lower.'" '.
        'upper="'.$upper.'">'."\n";
    
    foreach my $i (keys %{$nodes}) {
        $self->{LOG}{$label} .= "\t".'<level name="'.$i.'">'."\n";
        foreach my $j (sort{$a le $b} keys %{$nodes->{$i}}) {
            $self->{LOG}{$label} .= "\t\t".'<node>'."\n".
	    	"\t\t\t".'<count>'.$nodes->{$i}{$j}.'</count>'."\n".
		"\t\t\t".'<id>'.$j.'</id>'."\n".
		"\t\t".'</node>'."\n";
        }
        $self->{LOG}{$label} .= "\t".'</level>'."\n";
    }
    $self->{LOG}{$label} .= '</subnet>';

    push(@{$self->{QCOMMIT}},$label);
    return 1;
}

# LOG_SUBNETS
sub log_subnets {
    my $self = shift;
    my $subnets = shift;
    my $find = @_ ? shift : $self->{TARGET}{ID};
    
    my $re = qr/t(\d{2})[\_\-](\d{1,3})/;
    foreach my $i (keys %{$subnets}) {
        $i =~ m/$re/;
        my $lower = $1;
        my $upper = $2;
        
        print $i.' :: '.$lower.' :: '.$upper."\n";
        $self->log_subnet($subnets->{$i},$i,$lower,$upper,$find);
    }
    return 1;
}

####
# WRITE FUNC

# COMMIT
sub commit {
    my $self = shift;
    my $done = '';
    
    my $fh;
    open($fh,'>>'.$self->{DIR}{XML}.$self->{FILES}{XML});
    while(my $i = shift(@{$self->{QCOMMIT}})) {
        print $fh $self->{LOG}{$i};
        $done .= $i.', ';
    }
    print $fh "\n".'<!-- ['.$done.'] commited on '.$self->{DATE}{STAMP}.' -->'."\n";    
    close($fh);
    return 1;
}

# FINALIZE
sub finalize {
    my $self = shift;
    
    my $fh;
    open($fh,'>>'.$self->{DIR}{XML}.$self->{FILES}{XML});
    if($self->{QCOMMIT}) {
        print $fh "\n".'<!-- The following was not committed before finalizing... -->'."\n";
        $self->commit; # This works b/c of the comment at the end of commit
    }
    print $fh $self->{LOG}{FOOT};
    close($fh);
    return 1;
}

1;
__END__
