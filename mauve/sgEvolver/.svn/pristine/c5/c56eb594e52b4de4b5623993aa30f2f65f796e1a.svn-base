#!/usr/bin/env perl

# This script automates the process of generating a set of evolved
# sequences using specific evolution parameters for seq-gen and
# sgEvolver

#
# This script reads configuration parameters from a simujobparams.pm file
# The simujobparams.pm file must be manually edited and contains additional
# documentation.
#
# (c) 2004, 2005 Aaron Darling, All rights reserved
#

use strict;
use POSIX;
require simujobparams;


# do the work:
# generate job descriptions that can be used to 
# perform simulation, alignment, and scoring

my $total_runs = 0;
my $x_max = $simujobparams::x_start + $simujobparams::x_increment * $simujobparams::x_count;
my $y_max = $simujobparams::y_start + $simujobparams::y_increment * $simujobparams::y_count;

# perform a mutation rate hack for nucleotide substition:
# seq-gen doesn't like a 0 mutation rate so just set the
# rate to be abysmally small.
my $nt_sub_hack_rate = 0.000000000000001; # 10e-13 is less than 1 mutation in an entire human genome
if( $simujobparams::x_variable eq "nt_substitution" && $simujobparams::x_start <= 0 ){
	$simujobparams::x_start = $nt_sub_hack_rate;
}
if( $simujobparams::y_variable eq "nt_substitution" && $simujobparams::y_start <= 0 ){
	$simujobparams::y_start = $nt_sub_hack_rate;
}

$|++;	# this magical incantation causes stdout to be unbuffered
print "Generating simulation directories";

my $x_valI = 0;
my $y_valI = 0;
for( my $repI = 0; $repI < $simujobparams::repetitions; $repI++ ){
	my $x_val = $simujobparams::x_start;
	for( $x_valI = 0; $x_valI < $simujobparams::x_count; $x_valI++ ){ 
		$x_val += $simujobparams::x_increment if( $x_valI > 0 );
# modulate nt substitution rate on x-axis
		$simujobparams::nt_sub_scale = $x_val if $simujobparams::x_variable eq "nt_substitution";
		$simujobparams::indel_rate = $x_val if $simujobparams::x_variable eq "indels";
		$simujobparams::small_ht_rate = $x_val if $simujobparams::x_variable eq "small_ht";
		$simujobparams::large_ht_rate = $x_val if $simujobparams::x_variable eq "large_ht";
		$simujobparams::inv_rate = $x_val if $simujobparams::x_variable eq "inversions";
		$simujobparams::tree_scale = $x_val if $simujobparams::x_variable eq "tree_scale";

		my $y_val = $simujobparams::y_start;
		for( $y_valI = 0; $y_valI < $simujobparams::y_count; $y_valI++ ){
			$y_val += $simujobparams::y_increment if( $y_valI > 0 );
	                
			$simujobparams::nt_sub_scale = $y_val if $simujobparams::y_variable eq "nt_substitution";
	                $simujobparams::indel_rate = $y_val if $simujobparams::y_variable eq "indels";
                	$simujobparams::small_ht_rate = $y_val if $simujobparams::y_variable eq "small_ht";
        	        $simujobparams::large_ht_rate = $y_val if $simujobparams::y_variable eq "large_ht";
	                $simujobparams::inv_rate = $y_val if $simujobparams::y_variable eq "inversions";
			$simujobparams::tree_scale = $y_val if $simujobparams::y_variable eq "tree_scale";
			
			if( $total_runs < $simujobparams::skip_runs ){
				$total_runs++;
				next;
			}
			next if( $simujobparams::last_run >= 0 && $total_runs > $simujobparams::last_run );

		# get some random seed values
			my $seqgen_random = int(rand(0x7FFFFFFF));
			my $seqgen_donor_random = int(rand(0x7FFFFFFF));
			my $sgEvolver_random = int(rand(0x7FFFFFFF));

		# make a directory for the job files
			my $job_dir = "alignjob.$total_runs";
			mkdir( $job_dir );
			chdir( $job_dir );
			$total_runs++;

		# stat the ancestral donor to get its file size
		# select random subsequences for ancestral and donor
			my @stat_info = stat($simujobparams::ancestral_donor);
			my $size_modulus = $stat_info[7] - 2 * $simujobparams::sequence_length;
			my $ancestral_start = POSIX::floor( rand( $size_modulus ) );
			my $donor_modulus = $stat_info[7] - $ancestral_start - 2 * $simujobparams::sequence_length;
			my $donor_start = $ancestral_start + $simujobparams::sequence_length + POSIX::floor( rand( $donor_modulus ) );

		# write alignment params to a file
			open( PARAMS, ">simujobparams.pm"  );
			print PARAMS "package simujobparams;\n\n";
			
			print PARAMS "\$tools_dir=\"$simujobparams::tools_dir\";\n";
			print PARAMS "\$mauve_dir=\"$simujobparams::mauve_dir\";\n";
			print PARAMS "\$mavid_dir=\"$simujobparams::mavid_dir\";\n";
			print PARAMS "\$lagan_dir=\"$simujobparams::lagan_dir\";\n";
			print PARAMS "\$tba_dir=\"$simujobparams::tba_dir\";\n";		
	
			print PARAMS "\$seq_length=$simujobparams::sequence_length;\n";
			print PARAMS "\$nt_sub_scale=$simujobparams::nt_sub_scale;\n";
			print PARAMS "\$gamma_shape=$simujobparams::gamma_shape;\n";
			print PARAMS "\$indel_rate=$simujobparams::indel_rate;\n";
			print PARAMS "\$small_ht_rate=$simujobparams::small_ht_rate;\n";
			print PARAMS "\$small_ht_size=$simujobparams::small_ht_size;\n";
			print PARAMS "\$large_ht_rate=$simujobparams::large_ht_rate;\n";
			print PARAMS "\$large_ht_min=$simujobparams::large_ht_min;\n";
			print PARAMS "\$large_ht_max=$simujobparams::large_ht_max;\n";
			print PARAMS "\$inv_rate=$simujobparams::inv_rate;\n";
			print PARAMS "\$inv_size=$simujobparams::inv_size;\n";
			print PARAMS "\$nt_a_freq=$simujobparams::nt_a_freq;\n";
			print PARAMS "\$nt_c_freq=$simujobparams::nt_c_freq;\n";
			print PARAMS "\$nt_g_freq=$simujobparams::nt_g_freq;\n";
			print PARAMS "\$nt_t_freq=$simujobparams::nt_t_freq;\n";

			print PARAMS "\$tree_scale=$simujobparams::tree_scale;\n";
			print PARAMS "\$score_subset=\"$simujobparams::score_subset\";\n";

			print PARAMS "\$seq_count=\"$simujobparams::seq_count\";\n";
			print PARAMS "\@seqnames = (";
			for( my $i = 0; $i < $simujobparams::seq_count; $i++ )
			{
				print PARAMS "," if( $i > 0 );
				if( $i < scalar(@simujobparams::seqnames) )
				{
					print PARAMS "\"".$simujobparams::seqnames[$i]."\"";
				}else{
					print PARAMS "\"Tx$i\"";
				}
			}
			print PARAMS ")\;\n";

			print PARAMS "\$tree_filename=\"$simujobparams::tree_filename\";\n";
			print PARAMS "\$ancestral_donor=\"$simujobparams::ancestral_donor\";\n";
			print PARAMS "\$ancestral_seq_name=\"$simujobparams::ancestral_seq_name\";\n";
			print PARAMS "\$donor_seq_name=\"$simujobparams::donor_seq_name\";\n";
			print PARAMS "\$seqgen_out_name=\"$simujobparams::seqgen_out_name\";\n";
			print PARAMS "\$evolved_seqs_name=\"$simujobparams::evolved_seqs_name\";\n";
			print PARAMS "\$evolved_seqs_fname=\"$simujobparams::evolved_seqs_fname\";\n";

			print PARAMS "\$ancestral_start=$ancestral_start;\n";
			print PARAMS "\$donor_start=$donor_start;\n";

			print PARAMS "\$seqgen_random=$seqgen_random;\n";
			print PARAMS "\$seqgen_donor_random=$seqgen_donor_random;\n";
			print PARAMS "\$sgEvolver_random=$sgEvolver_random;\n";

			close PARAMS;
		
		# generate a tree file
			my $base_tree = $simujobparams::phylogeny;
			$base_tree = starTree( $base_tree, $simujobparams::seq_count, 1 ) if( $simujobparams::phylogeny eq "star" );
			my $scaled_tree = scaleTree( $base_tree, $simujobparams::tree_scale );
			open( TREEOUT, ">$simujobparams::tree_filename" );
			print TREEOUT $scaled_tree;
			close TREEOUT;
			
		# placate the user by printing a .
			print ".";
		# leave the directory
			chdir( ".." );
                # kludge for seq-gen's inability to accept a 0 mutation rate
	                if( $simujobparams::y_variable eq "nt_substitution" && $y_val == $nt_sub_hack_rate ){
        	                $y_val = 0;
                	}
		}
		if( $simujobparams::x_variable eq "nt_substitution" && $x_val == $nt_sub_hack_rate ){
			$x_val = 0;
		}
	}
}


createPbsScripts();
createCondorScripts();
createSGEScripts();

my $run_description_file = "run_description.txt";
open( RUN_FILE, ">$run_description_file" );
# print out the number of cols, rows, and the increments used for each
print RUN_FILE "$x_valI\t$y_valI\t$simujobparams::x_increment\t$simujobparams::y_increment\n";
close RUN_FILE;

print "\n\nJob generation complete.\n";
print "Edit the file mauveAlign.condor to test an aligner other than mauveAligner.\n";
print "Start jobs with something like 'condor_submit_dag -maxjobs 100 jobs.dag'.\n\n";
print "Or qsub sge.sh on Sun Grid Engine\n";

exit(0);

sub createPbsScripts
{

open( QSUB_FILE, ">qsub.sh" );
for( my $jobI = 0; $jobI < $total_runs; $jobI++ )
{
	print QSUB_FILE "cd alignjob.$jobI\n";
	print QSUB_FILE "qsub -wd -V -l jobfs=1GB,vmem=1GB,walltime=5:00:00 -k n -N $jobI$simujobparams::aligner runjob.sh\n";
	print QSUB_FILE "cd ..\n\n";
	open( JOBSH, ">alignjob.$jobI/runjob.sh" );
	print JOBSH "#!/bin/sh\n";
	print JOBSH "MYCWD=`pwd`\n";
	print JOBSH "cp * \$PBS_JOBFS/\n";
	print JOBSH "cd \$PBS_JOBFS\n";
	print JOBSH $simujobparams::tools_dir."/simujobrun.pl $simujobparams::aligner\n";
	print JOBSH "cd \$MYCWD\n";
	print JOBSH "cp \$PBS_JOBFS/* \$MYCWD/\n";
	close JOBSH;
	`chmod 755 alignjob.$jobI/runjob.sh`;
}
close QSUB_FILE;

`chmod 755 qsub.sh`;

}

sub createSGEScripts
{
my $sge_scratch_dir = $simujobparams::SGEscratch;
open( SGE_FILE, ">sge.sh" );
print SGE_FILE "\#!/bin/bash\n";
print SGE_FILE "\#\$ -cwd\n";
print SGE_FILE "\#\$ -V\n";
print SGE_FILE "\#\$ -S /bin/bash\n";
print SGE_FILE "\#\$ -t 1-$total_runs\n";
print SGE_FILE "LOCKFILE=`pwd`/lock\n";
print SGE_FILE "mkdir -p $sge_scratch_dir/aln\$SGE_TASK_ID\n";
print SGE_FILE "cd alignjob.\$SGE_TASK_ID\n";
print SGE_FILE "export CURDIR=\`pwd\`\n";
print SGE_FILE "while [ ! `mktemp -q \$LOCKFILE` ]; do\n";
print SGE_FILE "        sleep 10s\n";
print SGE_FILE "done\n";
print SGE_FILE "cp * $sge_scratch_dir/aln\$SGE_TASK_ID\n";
print SGE_FILE "rm \$LOCKFILE\n";
print SGE_FILE "cd $sge_scratch_dir/aln\$SGE_TASK_ID\n";
print SGE_FILE $simujobparams::tools_dir."/simujobrun.pl $simujobparams::aligner\n";
print SGE_FILE "mv * \$CURDIR\n";
print SGE_FILE "rm -rf ../aln\$SGE_TASK_ID\n";
close SGE_FILE;
`chmod 755 sge.sh`;
}


sub createCondorScripts
{

# print a condor dagman job list
# submit these jobs with the command "condor_submit_dag -maxjobs ## jobs.dag"
open( DAGMAN_FILE, ">jobs.dag" );
for( my $jobI = 0; $jobI < $total_runs; $jobI++ ){
	print DAGMAN_FILE "job simu$jobI mauveAlign.condor\n";
	print DAGMAN_FILE "vars simu$jobI aligndirname=\"alignjob.$jobI\"\n\n";
}
close DAGMAN_FILE;

chomp(my $current_dir = `pwd`);

	# write a condor job submission script			
my $condor_job_file = "mauveAlign.condor";
open( JOB_FILE, ">$condor_job_file" );
print JOB_FILE qq{
####################                   
# Condor job submission for mauveAligner
####################                                                   

universe       = vanilla
Executable     = $simujobparams::tools_dir/simujobrun.pl
# arguments should be one of none, mauve, progressiveMauve, mlagan, slagan, or mavid
arguments      = $simujobparams::aligner
Requirements   = Memory >= 256 && OpSys == "LINUX" && Arch =="INTEL"
Rank = Memory >= 512
#Image_Size     = 900 Meg                                                

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $simujobparams::tree_filename, simujobparams.pm
notification = Never
coresize = 1000
copy_to_spool = False

initialdir = \$(aligndirname)

Error   = simujobrun.err
Output  = simujobrun.out
Log     = $current_dir/simujobrun.log

Queue
};

close JOB_FILE;

}			

sub scaleTree
{
        my $tree = shift;
        my $scale = shift;
        my $out = "";
        # split tree text on : character, which always occurs before a branch length
        my @toks = split( /:/, $tree );
        return $tree if @toks == 1;
        for( my $tokI = 0; $tokI < @toks; $tokI++ )
        {
                $out .= ":" unless $tokI == 0;
                if( $toks[$tokI] =~ /(^.*)([\),])(.*)/ )
                {
                        my $len = $1;
                        $len *= $scale;
                        $out .= $len;
                        $out .= $2;
                        $out .= $3;
                }else{
                        $out .= $toks[$tokI];
                }
        }
        return $out;
}

sub starTree
{
	my $tree = shift;
	my $seq_count = shift;
	my $scale = shift;
	$tree = "";
	for( my $seqI = 0; $seqI < $seq_count - 1; $seqI++ )
	{
		$tree .= "(";
	}
	for( my $seqI = 0; $seqI < $seq_count; $seqI++ )
	{
		$tree .= "," unless( $seqI == 0 );
		$tree .= "Tx$seqI:$scale";
		$tree .= "):0" unless( $seqI == 0 );
	}
	$tree .= ";";
	return $tree;
}
