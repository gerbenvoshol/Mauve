#!/usr/bin/perl -w

# A script to generate a gradient graph of alignment scores using zimg
# Usage: gradientgen <data directory> <score file> <gradient output> [pixel scale]

use strict;

sub round { return int($_[0]) + (($_[0]*10%10 >= 5) ? 1 : 0); }

if( @ARGV != 3 && @ARGV != 4 ){
	print STDERR "Usage: Rgradientgen <data directory> <score file> <gradient output>\n";
	exit;
}

my $data_dir = $ARGV[0];
my $score_filename = $ARGV[1];
my $gradient_filename = "/tmp/scores.txt";
my $gradient_output = $ARGV[2];

#my $data_rows = 0;
#my $data_cols = 0;	# the number of x values
#my $x_increment = 0;
#my $y_increment = 0;

# read in the job description
my $run_description_file = "$data_dir/run_description.txt";
open( RUN_FILE, "$run_description_file" ) || die "Unable to open run description file $run_description_file\n";
my $cur_line = <RUN_FILE>;
chomp $cur_line;
(my $data_cols, my $data_rows, my $x_increment, my $y_increment ) = split( /\s+/, $cur_line );
close RUN_FILE;

open( SCOREFILE, $score_filename ) || die "Unable to open score file $score_filename\n";

my @data_values;
my @data_count;
for( my $dataI = 0; $dataI < $data_rows * $data_cols; $dataI++ ){
	$data_values[ $dataI ] = 0;
	$data_count[ $dataI ] = 0;
}
while( my $cur_line = <SCOREFILE> ){
	chomp $cur_line;
	if( length( $cur_line ) == 0 ){
		next;
	}
	my @cur_values = split( /\s+/, $cur_line );
	my $cur_x = round( $cur_values[0] / $x_increment );
	my $cur_y = round( $cur_values[1] / $y_increment );
	$data_values[ $cur_y * $data_cols + $cur_x ] += $cur_values[2];
	$data_count[ $cur_y * $data_cols + $cur_x ]++;
}

close( SCOREFILE );

# write out the zimg gradient input file
open( GRADIENTFILE, ">$gradient_filename" );

my $valueI = 0;
#for( my $rowI = $data_rows; $rowI > 0; $rowI-- ){
for( my $rowI = 0; $rowI < $data_rows; $rowI++ ){
	for( my $colI = 0; $colI < $data_cols; $colI++ ){
		my $valueI = $rowI * $data_cols + $colI;
		my $cur_value = $data_values[ $valueI ];
		if( $data_count[ $valueI ] > 0 ){
			$cur_value /= $data_count[ $valueI ];
		}
		print GRADIENTFILE $cur_value;
		print GRADIENTFILE "\n";
	}
}

close GRADIENTFILE;

chomp(my $cwd = `pwd`);
chdir($data_dir);
my $r_cmd = "R CMD BATCH ~/bin/rgradientplot.R";
system( $r_cmd );


