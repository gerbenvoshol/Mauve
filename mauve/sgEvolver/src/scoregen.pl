#!/usr/bin/perl -w

# Script to extract scores from a set of alignments and 
# generate an accuracy plot in postscript format
# uses the R package to generate postscript, requires R from http://r-project.org
# Usage: scoregen <score output file> [backbone]
# if collating scores on backbone regions use the backbone option
#

use strict;
require simujobparams;
use simujobparams;

if( @ARGV < 1 || @ARGV > 2 ){
	print STDERR "Usage: scoregen <postscript output base filename> \n";
	exit;
}

my $total_alignments = $simujobparams::x_count * $simujobparams::y_count;
my $score_infilename = "scores.txt";
my $plot_basename = $ARGV[0];
my $tmp_scorefile = "scorelist.txt";
my $fname_ext = "dat";
my $evolved_name = "evolved.dat";
print "x_var is: $simujobparams::x_variable, y_var is: $simujobparams::y_variable\n"; 


my @score_measures = ("Sum of pairs accuracy",
	"Sum of pairs positive predictive value",
	"Sum of pairs LCB sensitivity",
	"Sum of pairs LCB positive predictive value",
	"Absolute distance, min",
	"Absolute distance, first quartile",
	"Absolute distance, second quartile",
	"Absolute distance, third quartile",
	"Absolute distance, max",
	"Indel SP sensitivity",
	"Indel SP PPV",
	"mean",
	"variance"
);

my @output_filenames = (
"$plot_basename.nt_sn.ps",
"$plot_basename.nt_ppv.ps",
"$plot_basename.lcb_sn.ps",
"$plot_basename.lcb_ppv.ps",
"$plot_basename.bp_localization.ps",
"$plot_basename.indel_sn.ps",
"$plot_basename.indel_ppv.ps"
);

my @data0;
my @count0;
my @data1;
my @count1;
my @data2;
my @count2;
my @data3;
my @count3;
my @data4;
my @count4;
my @data5;
my @count5;
my @data6;
my @count6;
my @data7;
my @count7;
my @data8;
my @count8;
my @data9;
my @count9;
my @data10;
my @count10;
my @data11;
my @count11;

# create an array with the score values in it
my @data_values = (\@data0, \@data1, \@data2, \@data3, \@data4, \@data5, \@data6, \@data7, \@data8, \@data9, \@data10, \@data11);
my @data_count = (\@count0, \@count1, \@count2, \@count3, \@count4, \@count5, \@count6, \@count7, \@count8, \@count9, \@count10, \@count11);

# initialize a multidimensional array
for( my $mI = 0; $mI < @score_measures; $mI++ )
{
	for( my $dataI = 0; $dataI < $total_alignments; $dataI++ ){
		$data_values[$mI][$dataI] = 0;
		$data_count[$mI][$dataI] = 0;
	}
}

my $alignjob = -1;
for( my $repI = 0; $repI < $simujobparams::repetitions; $repI++ ){
print "Scanning replicate: $repI\n";
for( my $x_val = 0; $x_val < $simujobparams::x_count; $x_val++ ){
	for( my $y_val = 0; $y_val < $simujobparams::y_count; $y_val++ ){
	        $alignjob++;
		my $tmp_score = "alignjob.$alignjob/$score_infilename";
	        next if !(-e $tmp_score );
		
	        open( SCOREFILE, "$tmp_score" );

		# indexed as given in $score_measures
		my @scores = (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	        while( my $cur_line = <SCOREFILE> ){
			for( my $mI = 0; $mI < @score_measures; $mI++ )
			{
				if( $cur_line =~ /$score_measures[$mI]/ )
				{
					my @vals = split( /\:/, $cur_line);
					if($vals[1] =~ /\d+/)
					{
						$scores[$mI] = $vals[1];
						chomp $scores[$mI];
					}
				}
			}
	        }
	        close SCOREFILE;
		if( $scores[0] == -1 ){
			print "Warning: alignjob.$alignjob did not have a score.\n";
			next;
		}
		for( my $mI = 0; $mI < @score_measures; $mI++ )
		{
	        	$data_values[$mI][ $y_val * $simujobparams::y_count + $x_val ] += $scores[$mI];
	        	$data_count[$mI][ $y_val * $simujobparams::y_count + $x_val ]++;
		}
	}
}
}

for( my $mI = 0; $mI < 4; $mI++ )
{
	writeScoreFile( $tmp_scorefile, $mI );
	createGradientPlot( $output_filenames[$mI], "rgradientplot.R" );
}

# now create plots for the breakpoint localization
open( BPSCORES, ">$tmp_scorefile" );
for( my $dI = 0; $dI < $simujobparams::x_count * $simujobparams::y_count; $dI++ )
{
	if($data_count[4][$dI]==0){
		print BPSCORES "1000000\t1000000\t1000000\t1000000\t1000000";
	}
	my $yi = int($dI / $simujobparams::y_count);
	my $xi = $dI % $simujobparams::y_count;
	print BPSCORES "\n" if( $xi == 0 && $yi > 0 );
	for( my $sI = 4; $sI < 9; $sI++ )
	{
		print BPSCORES "\t" unless( $xi == 0 && $sI == 4 );
		if($data_count[$sI][$dI] > 0)
		{
			print BPSCORES ($data_values[$sI][$dI]/$data_count[$sI][$dI]);
		}else{
			print BPSCORES "0";
		}
	}
} 

createGradientPlot($output_filenames[4],"rbplocalization.R");

writeScoreFile( $tmp_scorefile, 9 );
createGradientPlot( $output_filenames[5], "rgradientplot.R" );
writeScoreFile( $tmp_scorefile, 10 );
createGradientPlot( $output_filenames[6], "rgradientplot.R" );


exit 0;

sub round { return int($_[0]) + (($_[0]*10%10 >= 5) ? 1 : 0); }

sub writeScoreFile
{
	my $score_fname = shift;
	my $mI = shift;
	# write out a temporary file with all the score values
	# in a form parseable by R 
	open( GRADIENTFILE, ">$score_fname" );

	my $valueI = 0;
	for( my $rowI = 0; $rowI < $simujobparams::x_count; $rowI++ ){
	        for( my $colI = 0; $colI < $simujobparams::y_count; $colI++ ){
	                my $valueI = $rowI * $simujobparams::y_count + $colI;
	                my $cur_value = $data_values[$mI][ $valueI ];
	                if( $data_count[$mI][ $valueI ] > 0 ){
	                        $cur_value /= $data_count[$mI][ $valueI ];
	                }
	                print GRADIENTFILE $cur_value;
	                print GRADIENTFILE "\n";
	        }
	}

	close GRADIENTFILE;
}

sub createGradientPlot
{
	my $output_fname = shift;
	my $plot_command = shift;
	$plot_command = "rgradientplot.R" if( $plot_command eq "" );
	open( RGP, "which $plot_command |" );
	my $rgradientplot = <RGP>;
	chomp $rgradientplot;
	close RGP;
	if( $rgradientplot eq "" || $rgradientplot =~ /No\ / )
	{
		$rgradientplot = "$simujobparams::tools_dir/$plot_command";
	}
	
	my $r_cmd = "/usr/bin/env R CMD BATCH $rgradientplot";
	my $retval = system( $r_cmd );
	if( $retval == 0 ){
		`mv scoreplot.ps $output_fname`;
		`rm $tmp_scorefile`;
	}else{
		print STDERR "An error occurred executing R CMD BATCH $rgradientplot\n";
		print STDERR "Check the file $plot_command"."out for details\n";
	}
}


