#!/usr/bin/perl
use strict;
use warnings;

# Script to process scores from breakpoint_calculator
# writes out a single tab text file with bp counts and simulation details
# (c) 2012 Aaron Darling
# Licensed under the GPL.

use strict;
require simujobparams;
use simujobparams;

if( @ARGV < 1 || @ARGV > 2 ){
	print STDERR "Usage: bp_scoring.pl <tabtext output file>\n";
	exit;
}

my $total_alignments = $simujobparams::x_count * $simujobparams::y_count;
my $output_basename = $ARGV[0];
print "x_var is: $simujobparams::x_variable, y_var is: $simujobparams::y_variable\n"; 

open(OUTDATA, ">$output_basename");
print OUTDATA "nt_sub_rate\tindel_rate\tsmall_ht_rate\tlarge_ht_rate\tinv_rate\ttrue_2way_bp\ttrue_3way_bp\tcalc_2way_bp\tcalc_3way_bp\n"; 
my $x_valI = 0;
my $y_valI = 0;
my $alignjob = -1;
for( my $repI = 0; $repI < $simujobparams::repetitions; $repI++ ){
	print "Scanning replicate: $repI\n";
	for( my $x_val = 0; $x_val < $simujobparams::x_count; $x_val++ ){
		for( my $y_val = 0; $y_val < $simujobparams::y_count; $y_val++ ){
			$alignjob++;
			next unless -e "alignjob.$alignjob/true_bp_2way.txt";
			do "alignjob.$alignjob/simujobparams.pm";
			my $outline = "";
			open(my $TRUE2, "alignjob.$alignjob/true_bp_2way.txt");
			open(my $TRUE3, "alignjob.$alignjob/true_bp_3way.txt");
			open(my $CALC2, "alignjob.$alignjob/calc_bp_2way.txt");
			open(my $CALC3, "alignjob.$alignjob/calc_bp_3way.txt");
			$outline .= $simujobparams::nt_sub_scale."\t";
			$outline .= $simujobparams::indel_rate."\t";
			$outline .= $simujobparams::small_ht_rate."\t";
			$outline .= $simujobparams::large_ht_rate."\t";
			$outline .= $simujobparams::inv_rate;
			$outline = add_token($TRUE2, $outline);
			$outline = add_token($TRUE3, $outline);
			$outline = add_token($CALC2, $outline);
			$outline = add_token($CALC3, $outline);
			print OUTDATA $outline."\n";;
		}
	}
}

sub add_token {
	my $FILE = shift;
	my $outline = shift;
	my $line = <$FILE>;
	chomp $line;
	my @toks = split( /\s+/, $line);
	$outline .= "\t".$toks[0];
	return $outline
;
}
