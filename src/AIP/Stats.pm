#!/usr/bin/perl -w

# Filename: Stats.pm
# Author: Stephen Bush

# Simple statistical functions
# Purpose: provide pearson correlation capability in perl.

package AIP::Stats;
require Exporter;

our @ISA      = qw( Exporter );
our @EXPORT   = qw( pcoeff pcoeff_matrix ); # insert sub names here
use vars qw( );
our $version  = 1.00;

use strict;
use warnings;
use Carp;

sub pcoeff {
    my $x = shift;
    my $y = shift;

	my $sum_sq_x = 0;
	my $sum_sq_y = 0;
	my $sum_x = 0;
	my $sum_y = 0;
	my $sum_coproduct = 0;
	my $index = $#{$x}; #last index
	my $total = $index + 1;

	for my $i (0..$index) {
		$sum_sq_x += $x->[$i]**2;
		$sum_sq_y += $y->[$i]**2;
		$sum_x += $x->[$i];
		$sum_y += $y->[$i];
		$sum_coproduct += $x->[$i] * $y->[$i];
	}

	my $mean_x = $sum_x / $total;
	my $mean_y = $sum_y / $total;
	
	my $var_x = (1 / $index) * $sum_sq_x - ($total / $index) * $mean_x**2;
	my $var_y = (1 / $index) * $sum_sq_y - ($total / $index) * $mean_y**2;
	my $sd_x = sqrt($var_x);
	my $sd_y = sqrt($var_y);
	my $covar = (1 / $index) * ( $sum_coproduct - (($sum_x * $sum_y) / $total));
	my $correlation = $covar / ($sd_x * $sd_y);
	
	my $ss_x = $sum_sq_x - $sum_x**2/$total;
	my $ss_y = $sum_sq_y - $sum_y**2/$total;
	my $ss_co = $sum_coproduct - ($sum_x * $sum_y) / $total;
	my $tcorr = $ss_co / sqrt( $ss_x * $ss_y);
	
	my $cc = 0;
	for my $j (0..$index) {
		$cc += (($x->[$j] - $mean_x) / $sd_x) *  (($y->[$j] - $mean_y) / $sd_y);
	}
	$cc *= (1/$total);
	
	#if ($tcorr != $correlation) { die('no good'."\n"); }
	
	return $cc;
	
	#for my $i (1..$N) {
	#	$sweep = ($i - 1.0) / $i;
	#	$delta_x = $x->[$i] - $mean_x;
	#	$delta_y = $y->[$i] - $mean_y;
	#	$sum_sq_x += $delta_x * $delta_x * $sweep;
	#	$sum_sq_y += $delta_y * $delta_y * $sweep;
	#	$sum_coproduct += $delta_x * $delta_y * $sweep;
	#	$mean_x += $delta_x / $i;
	#	$mean_y += $delta_y / $i;
	#}
	#$N++; # total count
	#my $pop_sd_x = sqrt( $sum_sq_x/$N );
	#my $pop_sd_y = sqrt( $sum_sq_y/$N );
	#my $cov_x_y = $sum_coproduct/$N;
	#my $correlation = $cov_x_y / ($pop_sd_x * $pop_sd_y);
	#
	#return $correlation;
	
    #my $t = @{$d->[0]};
    #my $ss = [0,0,0,0,0];
    #
    #for(my $i = 0; $i < $t; $i++) {
    #    $ss->[0] += $d->[0][$i];
    #    $ss->[1] += $d->[0][$i]**2;
    #    $ss->[2] += $d->[1][$i];
    #    $ss->[3] += $d->[1][$i]**2;
    #    $ss->[4] += $d->[0][$i] * $d->[1][$i];
    #}
    #
    #return ($ss->[4] - ($ss->[0]*$ss->[2])/$t) /
    #    (
    #        (($ss->[1] - (($ss->[0]**2)/$t))**(1/2)) *
    #        (($ss->[3] - (($ss->[2]**2)/$t))**(1/2))
    #    );
}

sub pcoeff_matrix {
    my $m = shift;
	my $t = []; # transvers
	my $pc = []; # pearson

	#map {
	#	map {
	#		print $_.' ';
	#	} @{$_};
	#	print "\n";
	#} @{$m};
	
	
    my $size = scalar @{$m} - 1; # index last element

	# build transverse array
	foreach my $i (0..$size) {
		foreach my $j (0..$size) {
			$t->[$j] = [] if (!exists $t->[$j]);
			$t->[$j][$i] = $m->[$i][$j];
		}
	}
	#print  'm'.$#{$m}.' '.$#{$m->[0]}."\n";
	#print  't'.$#{$t}.' '.$#{$t->[0]}."\n";

    foreach my $i (0..$size)  {
        foreach my $j (0..$size) {
            #map { $temp->[$_] = $m->[$_][$i];  } (0 .. $size) ;
            #map { print sprintf("%.2f,",$temp->[$_]); } (0 .. $size);
            #print "\n";
            #print $i.' '.$j."\n";
            
            $pc->[$i][$j] = pcoeff(
                $m->[$i],
                $t->[$j]
            );
            #$pc->[$j][$i] = $pc->[$i][$j];
        }
    }
    
    return $pc;
}

1;

__END__


