#!/usr/bin/env perl

# This script automates the process of generating a set of evolved
# sequences using specific evolution parameters for seq-gen and
# sgEvolver, aligning those sequences, and scoring the alignments


use strict;
use POSIX;
require simujobparams;

if( @ARGV < 1 || @ARGV > 2 ){
	die "Usage: simujobrun.pl <mauve|progressiveMauve|mavid|mlagan|tba|none> [debug]";
}

# Set the location of the evolution and scoring tools here:
# this is the directory where seq-gen, sgEvolver, and scoreAlignment reside
my $xmfa_alignment = "";

# an array of files that need to be deleted when cleaning up
my @delete_files = ();

# check whether we're debugging...
my $debug = 0;
$debug = 1 if @ARGV > 1 && $ARGV[1] eq "debug";

# ROADMAP:
# 1) simulate evolution according to parameters
# 2) align genomes using the selected aligner
# 3) score alignments


		# write any command lines executed to a file
open( COMLINES, ">command_lines.txt"  );


		# generate an ancestral base
open( ANCESTRAL, ">$simujobparams::ancestral_seq_name" );
print ANCESTRAL " 1 $simujobparams::seq_length\n";
print ANCESTRAL "Ancestral\n";
close ANCESTRAL;
	
my $rval;

		# extract the specified ancestral sequence
my $extract_cl = "dd if=$simujobparams::ancestral_donor bs=1 skip=$simujobparams::ancestral_start count=$simujobparams::seq_length";
executeCommand( $extract_cl, ">$simujobparams::ancestral_seq_name", "ancestral_dd.err" );

		# remove any stale output file
my $rm_cl = "rm -f $simujobparams::seqgen_out_name";
executeCommand( $rm_cl, "rm.out", "rm.err" );

		# generate a data set
my $seqgen_cl = $simujobparams::tools_dir."seq-gen -mHKY -a $simujobparams::gamma_shape -i 0 -z $simujobparams::seqgen_random -f ";
$seqgen_cl .= "$simujobparams::nt_a_freq $simujobparams::nt_c_freq $simujobparams::nt_g_freq $simujobparams::nt_t_freq ";
$seqgen_cl .= "-t 4 -k 1 -wa -on -s $simujobparams::nt_sub_scale $simujobparams::tree_filename ";
$seqgen_cl .= "< $simujobparams::ancestral_seq_name ";
executeCommand( $seqgen_cl, ">$simujobparams::seqgen_out_name", "seqgen_ancestral.err" );

		# generate a donor base
open( DONOR, ">$simujobparams::donor_seq_name" );
print DONOR " 1 $simujobparams::seq_length\n";
print DONOR "Donor\n";
close DONOR;


		# extract the specified donor sequence
$extract_cl = "dd if=$simujobparams::ancestral_donor bs=1 skip=$simujobparams::donor_start count=$simujobparams::seq_length";
executeCommand( $extract_cl, ">$simujobparams::donor_seq_name", "donor_dd.err" );

		# generate a donor data set
$seqgen_cl = $simujobparams::tools_dir."seq-gen -mHKY -a $simujobparams::gamma_shape -i 0 -z $simujobparams::seqgen_donor_random  -f ";
$seqgen_cl .= "$simujobparams::nt_a_freq $simujobparams::nt_c_freq $simujobparams::nt_g_freq $simujobparams::nt_t_freq ";
$seqgen_cl .= "-t 4 -k 1 -wa -on -s $simujobparams::nt_sub_scale $simujobparams::tree_filename ";
$seqgen_cl .= "< $simujobparams::donor_seq_name ";
executeCommand( $seqgen_cl, ">$simujobparams::seqgen_out_name", "seqgen_donor.err" );


		# continue evolution with sgEvolver
my $sgevolver_cl = $simujobparams::tools_dir."sgEvolver --indel-freq=".$simujobparams::indel_rate." --small-ht-freq=$simujobparams::small_ht_rate --small-ht-size=$simujobparams::small_ht_size".
" --large-ht-freq=$simujobparams::large_ht_rate --inversion-freq=$simujobparams::inv_rate --large-ht-min=$simujobparams::large_ht_min --large-ht-max=$simujobparams::large_ht_max".
" --random-seed=$simujobparams::sgEvolver_random".
" --inversion-size=$simujobparams::inv_size $simujobparams::tree_filename $simujobparams::seqgen_out_name $simujobparams::evolved_seqs_name $simujobparams::evolved_seqs_fname";
$rval = executeCommand( $sgevolver_cl, "sgEvolver.out", "sgEvolver.err" );

		# die if sgEvolver failed
die "Failure in sgEvolver" if( $rval != 0 );

die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_name";
die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_fname";

		# delete seq-gen related files
push( @delete_files, $simujobparams::seqgen_out_name );
push( @delete_files, $simujobparams::ancestral_seq_name );
push( @delete_files, $simujobparams::donor_seq_name );

deleteFiles( @delete_files );

#
# STEP 2:  align genomes using the selected aligner
#

	# collect timing information
	# NOTE:  timing info measures wall clock time and includes
	# any time the job was preempted--not good when using a condor
	# pool.  It should be valid when run on a single, otherwise unused machine
	# We can't directly measure cpu time (e.g. using /usr/bin/time) because
	# the various aligners may have wrapper scripts that use little CPU time
my $start_time = time();

alignMauve() if( $ARGV[0] eq "mauve" );
alignMavid() if( $ARGV[0] eq "mavid" );
alignMlagan() if( $ARGV[0] eq "mlagan" );
alignSlagan() if( $ARGV[0] eq "slagan" );
alignProgressiveMauve() if( $ARGV[0] eq "progressiveMauve" );
alignTBA() if( $ARGV[0] eq "tba" );
exit(0) if( $ARGV[0] eq "none" );
my $end_time = time();

open( ALN_TIME_FILE, ">alignment_time.txt" );
print ALN_TIME_FILE ($end_time - $start_time);
print ALN_TIME_FILE " seconds\n";
close ALN_TIME_FILE;

#
# STEP 3: score the alignment
#

my $score_lcb_arg = "";
unless( $ARGV[0] eq "mauve" || $ARGV[0] eq "progressiveMauve" )
{
	$score_lcb_arg = "--disable-lcb-scoring";
}

# project to a subset if the user requested it
unless( $simujobparams::score_subset eq "" )
{
	my $projector_path = $simujobparams::tools_dir."alignmentProjector";
	my $project_evolved_cl = "$projector_path $simujobparams::evolved_seqs_name projected_$simujobparams::evolved_seqs_name "
				." $simujobparams::evolved_seqs_fname projected_$simujobparams::evolved_seqs_fname "
				.$simujobparams::score_subset; 
	my $project_aligned_cl = "$projector_path $xmfa_alignment projected_$xmfa_alignment "
				." $simujobparams::evolved_seqs_fname projected_$simujobparams::evolved_seqs_fname "
				.$simujobparams::score_subset;
	executeCommand( $project_evolved_cl, "project_evolved.out", "project_evolved.err" );
	executeCommand( $project_aligned_cl, "project_aligned.out", "project_aligned.err" ); 
	push( @delete_files, $simujobparams::evolved_seqs_name );
	push( @delete_files, "$simujobparams::evolved_seqs_fname" );
	push( @delete_files, $xmfa_alignment );
	$simujobparams::evolved_seqs_name = "projected_$simujobparams::evolved_seqs_name";
	$xmfa_alignment = "projected_$xmfa_alignment";
	$simujobparams::evolved_seqs_fname = "projected_$simujobparams::evolved_seqs_fname";
}

my $scoreAlignment_path = $simujobparams::tools_dir."scoreAlignment2";

my $score_cl = "$scoreAlignment_path $simujobparams::evolved_seqs_name $xmfa_alignment $simujobparams::evolved_seqs_fname $score_lcb_arg";
executeCommand( $score_cl, "scores.txt", "scoring.err" );

		# extract "backbone" and do a backbone scoring
my $bb_file = "backbone.dat";
push( @delete_files, $bb_file );
my $extractbb_cl = $simujobparams::tools_dir."extractBackbone $simujobparams::evolved_seqs_fname $simujobparams::evolved_seqs_name 50 50 $bb_file";
#executeCommand( $extractbb_cl, "extract_bb.out", "extract_bb.err" );

$score_cl = "$scoreAlignment_path $bb_file $xmfa_alignment $simujobparams::evolved_seqs_fname $score_lcb_arg";
#executeCommand( $score_cl, "bb_scores.txt", "bb_scoring.err" );

		# delete evolved sequence data
push( @delete_files, $simujobparams::evolved_seqs_name );
push( @delete_files, $simujobparams::evolved_seqs_fname );

deleteFiles( @delete_files );
exit(0);

sub executeCommand {
  my $command = shift;
  my $stdout_file = shift;
  my $stderr_file = shift;
  $command .= " >$stdout_file 2>$stderr_file";
  print "Executing $command\n";
  print COMLINES "$command\n";
  my $rval = system($command);
  `echo "Exited with code $rval" >> $stderr_file` if $rval != 0;
  return $rval;
}

sub deleteFiles {
  if( $debug != 1 ){
    foreach(@_){
      `rm -f $_`;
      print COMLINES "rm -f $_\n";
    }
  }
}

sub alignMauve {
  $xmfa_alignment = "aligned.dat";
  push( @delete_files, $xmfa_alignment );
  push( @delete_files, "aligned.mums" );
  push( @delete_files, "aligned.mauve" );
  my $mauve_cl = $simujobparams::mauve_dir."mauveAligner --mums --output aligned.mums $simujobparams::evolved_seqs_fname";
  executeCommand( $mauve_cl, "mums.out", "mums.err" );

  $mauve_cl = $simujobparams::mauve_dir."mauveAligner --match-input aligned.mums --output aligned.mauve ".
  "--output-alignment=$xmfa_alignment $simujobparams::evolved_seqs_fname";
  $rval = executeCommand( $mauve_cl, "aligner.out", "aligner.err" );

  # die if mauve failed
  die "Failure in aligner" if( $rval != 0 );
}

sub alignMavid {
  $xmfa_alignment = "mavid.dat";

  # dust the sequence for mavid
  my $dust_cl = "$simujobparams::mavid_dir/mdust $simujobparams::evolved_seqs_fname";
  executeCommand( $dust_cl, "$simujobparams::evolved_seqs_fname.masked", "mdust.err" );

  my $mavid_cl = $simujobparams::mavid_dir."/mavid $simujobparams::tree_filename $simujobparams::evolved_seqs_fname";
  executeCommand( $mavid_cl, "mavid.out", "mavid.err" );


  # step 2.1 convert mavid.mfa to xmfa format
  my $xmfa_cl = "$simujobparams::tools_dir/mfa2xmfa mavid.mfa $xmfa_alignment";
  executeCommand( $xmfa_cl, "mfa2xmfa.out", "mfa2xmfa.err" );

  push( @delete_files, "$simujobparams::evolved_seqs_fname.masked" );
  push( @delete_files, "mavid.mfa" );
  push( @delete_files, "mavid.phy" );
  push( @delete_files, "$xmfa_alignment" );

}

sub alignProgressiveMauve {
  $xmfa_alignment = "proaligned.dat";
  push( @delete_files, $xmfa_alignment );
  push( @delete_files, "proaligned.mums" );
  push( @delete_files, "proaligned.mauve" );
  my $mauve_cl = $simujobparams::mauve_dir."progressiveMauve --mums --output=proaligned.mums $simujobparams::evolved_seqs_fname";
  executeCommand( $mauve_cl, "promums.out", "promums.err" );

  $mauve_cl = $simujobparams::mauve_dir."progressiveMauve --match-input=proaligned.mums --output=$xmfa_alignment ".
  " $simujobparams::evolved_seqs_fname";
  $rval = executeCommand( $mauve_cl, "proaligner.out", "proaligner.err" );

  # die if mauve failed
  die "Failure in aligner" if( $rval != 0 );

}

sub alignMlagan {
	$xmfa_alignment = "mlagan.dat";
	push( @delete_files, $xmfa_alignment );

# split up the input sequences into one file per sequence
	my $mfa_to_multi_cl = "$simujobparams::tools_dir/mfaToMultiFiles $simujobparams::evolved_seqs_fname";
	my $rval = executeCommand( $mfa_to_multi_cl, "mfa_to_multi.out", "mfa_to_multi.err" );
	die "Failure in mfaToMultiFiles" if( $rval != 0 );
# read the guide tree and remove branch lengths and commas
	open( TREEFILE, "test.tree" );
	my $tree = <TREEFILE>;
	$tree = makeTBAtree($tree);
	close TREEFILE;
# set the path for lagan
	$ENV{"LAGAN_DIR"} = $simujobparams::lagan_dir;

	push( @delete_files, "mlagan.out" );
	my $mlagan_cl = "$simujobparams::lagan_dir/mlagan Taxon* -tree \"$tree\"";
        $rval = executeCommand( $mlagan_cl, "mlagan.out", "mlagan.err" );
        die "Failure in mlagan" if( $rval != 0 );
	`rm Taxon*`;
	`rm *anch*`;	
	my $xmfa_cl = "$simujobparams::tools_dir/mfa2xmfa mlagan.out $xmfa_alignment";
	$rval = executeCommand( $xmfa_cl, "mfa2xmfa.out", "mfa2xmfa.err" );
	die "Failure in mfa2xmfa" if( $rval != 0 );
}

sub alignSlagan {

}

sub alignTBA {
	$xmfa_alignment = "tbaaligned.dat";
# split up the input sequences into one file per sequence
	my $mfa_to_multi_cl = "mfaToMultiFiles $simujobparams::evolved_seqs_fname";
	my $rval = executeCommand( $mfa_to_multi_cl, "mfa_to_multi.out", "mfa_to_multi.err" );
	die "Failure in mfaToMultiFiles" if( $rval != 0 );

# read the guide tree and remove branch lengths and commas
	open( TREEFILE, "test.tree" );
	my $tree = <TREEFILE>;
	$tree = makeTBAtree($tree);
	close TREEFILE;
	my $pathenv = $ENV{"PATH"};
	print "path is: $pathenv\n\n";
	$ENV{"PATH"} = $simujobparams::tba_dir.":$pathenv";
# run all_bz
	my $all_bz_cl = "all_bz \"$tree\"";
	$rval = executeCommand( $all_bz_cl, "all_bz.out", "all_bz.err" );
	die "Failure in all_bz" if( $rval != 0 );	
# run tba
	my $tba_cl = "tba \"$tree\" *.*.maf tba.maf";
	$rval = executeCommand( $tba_cl, "tba.out", "tba.err" );
	die "Failure in TBA" if( $rval != 0 );
	push( @delete_files, "*.*.maf" );
	push( @delete_files, "tba.maf" );
	my $maf2xmfa_cl = "maf2xmfa tba.maf $xmfa_alignment";
	for( my $i = 0; $i < $simujobparams::seq_count; $i++ )
	{
		$maf2xmfa_cl .= " ".$simujobparams::seqnames[$i];
		push( @delete_files, $simujobparams::seqnames[$i] );
	}
	$rval = executeCommand( $maf2xmfa_cl, "maf2xmfa.out", "maf2xmfa.err" );
	die "Failure in maf2xmfa" if( $rval != 0 );
}

sub makeTBAtree
{
        my $tree = shift;
        my $out = "";
        # split tree text on : character, which always occurs before a branch length
        my @toks = split( /:/, $tree );
        return $tree if @toks == 1;
        for( my $tokI = 0; $tokI < @toks; $tokI++ )
        {
                if( $toks[$tokI] =~ /(^.*)([\),])(.*)/ )
                {
                        my $len = $1;
                        $out .= $2;
                        $out .= $3;
                }else{
                        $out .= $toks[$tokI];
                }
        }
	# replace commas with spaces
	$out =~ s/,/ /g;
	# get rid of semicolon
	$out =~ s/\;//g;
        return $out;
}

