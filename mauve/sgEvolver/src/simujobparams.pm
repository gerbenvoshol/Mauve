package simujobparams;


#
# This configuration file specifies simulation parameters for the simujobgen
# and simujobrun perl scripts.  To perform the desired simulation, edit the 
# parameters below
#
# After editing this configuration file, run simujobgen.pl in the same directory
# as this config file.   simujobgen.pl will generate a directory for each 
# requested simulation that contains a specification
# file describing the simulation to be done.  Use the 'simujobrun.pl' script to
# actually perform the simulation, alignment, and scoring of the alignment
#
# (c) 2004, 2005 Aaron Darling, All rights reserved
#

# set these paths to the locations of sgEvolver (tools_dir), 
# mauveAligner (mauve_dir), mavid (mavid_dir), and Multi/Shuffle-
# LAGAN (lagan_dir).  mavid and lagan are not required.
# the scripts expect seq-gen, scoreAlignment, extractBackbone,
# and rgradientplot.R to be installed in $tools_dir
$tools_dir="/home/koadman/bin/";
$mauve_dir="/home/koadman/bin/";
$mavid_dir="/home/koadman/bin/";
$lagan_dir="/u/d/a/darling/lagan121/";
$tba_dir="/u/d/a/darling/tba";

# set this value to the name of the aligner to be tested
# currently implemented values are mauve, progressiveMauve, tba, mlagan, and mavid
$aligner = "tba";

# set this to the number of replicate simulations to perform.
# Each simulation will have the same average mutation parameters 
# but uses a different random seed

$repetitions = 3; 


# set this to the length of genomes to simulate
# the source sequence file must be at least twice this length

$sequence_length = 1000000;


# the path to a file containing source sequence to use as the ancestor
# at the root of the tree.  Two random, non-overlapping, stretches of 
# $sequence_length characters will be chosen from this file.  One gets
# used as the ancestral sequence for the simulated genomes, the other
# gets used as ancestral sequence in a second identical simulation from 
# which sequence to be inserted gets drawn (the donor sequences).
#
# this source file must be a raw sequence file.  e.g. it contains no
# FastA headers, no newlines, nothing but characters in {A,C,G,T,a,c,g,t}

$ancestral_donor = "/home/koadman/sequence/o157.raw";

#
# the scale and rate parameters scale the branch lengths of the 
# phylogenetic tree to arrive at the average number of events
# per site.  for example, an nt_sub_scale of .05 and a branch 
# length of .3 will result in .015 nucleotide substitutions per
# base along that branch.  average rates for other mutation events 
# are calculated similarly

$nt_sub_scale = 0.02;
$indel_rate = 0;
$small_ht_rate = 0;
$small_ht_size = 200;
$large_ht_rate = 0;
$large_ht_min = 10000;
$large_ht_max = 60000;
$inv_rate = 0;
$inv_size = 50000;

#
# all mutation rates can be scaled simultaneously by scaling the tree branch lengths
#
$tree_scale = 1;

#
# specify which parameters to modulate when generating an alignment profile
# two parameters can be modulated and their values will ultimately be
# plotted on the X and Y axes of an aligner profile
#
# x_variable and y_variable can take possible values of:
# {nt_substitution,indels,small_ht,large_ht,inversions,tree_scale}

$x_variable = "nt_substitution";
$x_start = 0;        # the starting mutation rate for the x axis variable

$x_increment = 0.05; # the mutation rate is incremented by this much for 
                        # each successive simulation
$x_count = 20;    # increment the mutation rate this many times
                     # the total number of simulations will be 
                     # $x_count * $y_count * $repetitions, and the number of
                     # points on the alignment profile will be $x_count * $y_count
$y_variable = "indels";
$y_start = 0;
$y_increment = 0.01;
$y_count = 20;


# gamma shape parameter for an uneven distribution
# of nucleotide substitution sites

$gamma_shape = 1;

# background nucleotide frequencies
$nt_a_freq = .22;     # passed to Seq-gen
$nt_c_freq = .28;     # passed to Seq-gen
$nt_g_freq = .28;     # passed to Seq-gen
$nt_t_freq = .22;     # passed to Seq-gen

# need at least one sequence name entry per sequence...
@seqnames = ("Taxon1","Taxon2","Taxon3","Taxon4","Taxon5","Taxon6","Taxon7","Taxon8","Taxon9");
$seq_count = 2;    # set this to the number of organisms in the phylogeny

# the phylogenetic tree relating the organisms
# The tree topology and branch lengths must be given in Newick format
# inferred tree for 9 enterobacteria
$phylogeny = "(((($seqnames[0]:0.12979,($seqnames[1]:0.06944,$seqnames[2]:0.06473)".
             ":0.09078):0.00963,($seqnames[3]:0.09309,$seqnames[4]:0.09051)".
             ":0.07744):0.03698,$seqnames[5]:0.19844):0.2646,($seqnames[6]".
             ":0.11243,($seqnames[7]:0.06886,$seqnames[8]:0.06788):0.04328):0.35968);";

# an example phylogeny for two taxa
# $phylogeny="($seqnames[0]:.5,$seqnames[1]:.5);";

# score only a subset of the sequences, index starts at 0
# leave empty to score all sequences
# $score_subset = "0 1 2";
$score_subset = "";

#
# the rest of these parameters usually need not be modified
#
$tree_filename = "test.tree";
$ancestral_seq_name = "ancestral.dat";
$donor_seq_name = "donor.dat";
$seqgen_out_name = "seqgen.dat";
$evolved_seqs_name = "evolved.dat";	# file for the fully evolved alignment output
$evolved_seqs_fname = "evolved_seqs.fas";	# file for the ungapped sequences

$skip_runs = 0;	# the first data set to generate.  
                # useful when adding additional runs to an existing simulation data set
                # For example, when another replicate of the experiment is needed

$last_run = -1;	# the last data set to generate, set to negative for a run to completion


