#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <stack>
#include <cstring>
#include <ctime>
#include "getopt.h"
#include "libGenome/gnSequence.h"
#include "libMems/gnAlignedSequences.h"
#include "Mutator.h"
#include "Alignment.h"
#include <boost/algorithm/string.hpp>

extern "C" {
#include "twister.h"
}

using namespace std;
using namespace genome;
using namespace mems;

#ifdef _MSC_VER
/** Visual Studio doesn't define strtoll() */
int64 strtoll( const char* a, char ** b, int c){
	return strtol( a, b, c );
}
#endif

class MutationEvent {
public:
	MutationEvent( const char* mut_name, double mut_frequency, Mutator& mut_mutator );
	MutationEvent( const MutationEvent& me );
	MutationEvent& operator=( const MutationEvent& me );
	string name;
	double frequency;
	Mutator& mutator;
private:
	MutationEvent();
};

MutationEvent::MutationEvent( const char* mut_name, double mut_frequency, Mutator& mut_mutator ) :
name( mut_name ),
frequency( mut_frequency ),
mutator( mut_mutator )
{
}

MutationEvent::MutationEvent( const MutationEvent& me ) :
name( me.name ),
frequency( me.frequency ),
mutator( me.mutator )
{
}

MutationEvent& MutationEvent::operator=( const MutationEvent& me )
{
	name = me.name;
	frequency = me.frequency;
	return *this;
}


class AlignmentEvolver {
public:
	AlignmentEvolver( vector< MutationEvent >& mutation_events );
	void evolve( const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment, int64 random_seed = -1 );
protected:
	vector< MutationEvent >& mut_events;
};

AlignmentEvolver::AlignmentEvolver( vector< MutationEvent >& mutation_events ) :
mut_events( mutation_events )
{
}

void AlignmentEvolver::evolve( const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment, int64 random_seed )
{	
	Alignment& ev = evolved_alignment;
	
	// scale the mutation frequencies to the tree length,
	// assuming the branch lengths given are expected number of
	// substitutions per site.
	
	double tree_length = tree.getHeight();
	
	// scale mutation rates by tree length
	uint mutI;
	for( mutI = 0; mutI < mut_events.size(); mutI++ ){
		mut_events[ mutI ].frequency *= tree_length;
	}

	// get the total mutation rate so we can sample mutation event types from a distribution
	double total_mutation_rate = 0;
	for( mutI = 0; mutI < mut_events.size(); mutI++ )
		total_mutation_rate += mut_events[ mutI ].frequency;
	
	double total_branch_len = 0;
	vector< double > branch_lengths;
	node_id_t nodeI = 0;
	for( nodeI = 0; nodeI < tree.size(); nodeI++ ){
		branch_lengths.push_back( tree[ nodeI ].distance );
		total_branch_len += tree[ nodeI ].distance;
	}

	// calculate the total number of events
//  # of events should be per base pair per unit time
// if the user wants mutation rates to be independent of nucleotide substitution rates an alternate tree can be used for seg-gen	
	gnSeqI total_events = (gnSeqI)(total_mutation_rate * total_branch_len * ev.length());
	
	// normalize branch lengths to a distribution
	for( nodeI = 0; nodeI < tree.size(); nodeI++ )
		branch_lengths[ nodeI ] /= total_branch_len;
	
	// normalize event frequencies to a distribution
	for( mutI = 0; mutI < mut_events.size(); mutI++ )
		mut_events[ mutI ].frequency /= total_mutation_rate;
	
	// sample branches for events
	vector< vector< uint > > branch_events( branch_lengths.size() );
	gnSeqI eventI;
	double randy, sum;
	uint typeI;
	// this should be random enough, just don't evolve two trees in the same second...
	uint seed = random_seed == -1 ? time(NULL) : random_seed;
//	seed = 1063558787;
	cout << "Random seed: " << seed << endl;
	SetSeed( seed );
	for( eventI = 0; eventI < total_events; eventI++ ){
		randy = rndu();
		sum = 0;
		for( nodeI = 0; nodeI < branch_lengths.size(); nodeI++ ){
			sum += branch_lengths[ nodeI ];
			if( sum >= randy )
				break;
		}
		// watch out for floating point imprecision
		if( nodeI == branch_lengths.size() )
			nodeI--;
		
		// so the event will occur on branch 'nodeI'.
		// now sample the event type
		randy = rndu();
		sum = 0;
		for( typeI = 0; typeI < mut_events.size(); typeI++ ){
			sum += mut_events[ typeI ].frequency;
			if( sum >= randy )
				break;
		}
		// watch out for floating point imprecision
		if( typeI == mut_events.size() )
			typeI--;

		branch_events[ nodeI ].push_back( typeI );
	}

	// now perform mutations in tree-branching order
	stack< node_id_t > node_stack;
	stack< node_id_t > child_stack;
	node_stack.push( tree.root );

	while( node_stack.size() > 0 ){
		// if we just recursed to a new child:
		if( node_stack.size() != child_stack.size() ){
			// perform mutations on this node
			for( eventI = 0; eventI < branch_events[ node_stack.top() ].size(); eventI++ ){
//				if( eventI == 319 && node_stack.top() == 4 )
//					cerr << "boom\n";
				if( debugCheckingLevel() == 1 && eventI % 100 == 0 ){
					debugChecks( 1 );
					evolved_alignment.checkLengths();
				}
				if( eventI % 500 == 0 )
					cout << "branch: " << node_stack.top() << " event: " << eventI << " " << 
						mut_events[ branch_events[ node_stack.top() ][ eventI ] ].name << "\n";
				mut_events[ branch_events[ node_stack.top() ][ eventI ] ].mutator.mutate( 
						node_stack.top(), tree, evolved_alignment );
				if( debugCheckingLevel() == 1 && eventI % 100 == 0 )
					debugChecks( 0 );
			}
			
			// push a child if possible
			if( tree[ node_stack.top() ].children.size() > 0 ){
				node_stack.push( tree[ node_stack.top() ].children[ 0 ] );
				child_stack.push( 1 );
			}else	// done with this leaf node
				node_stack.pop();
			continue;
		}
		if( child_stack.top() != tree[ node_stack.top() ].children.size() ){
			// push a child
			node_stack.push( tree[ node_stack.top() ].children[ child_stack.top() ] );
			child_stack.top()++;
			continue;
		}
		// no more children under this parent node
		node_stack.pop();
		child_stack.pop();
	}
	
	if( debugCheckingLevel() > 0 )
		debugChecks( 1 );
	// finally apply all the inversion events
	evolved_alignment.applyInversions( tree );
}

void print_usage( const char* pname ){
	cerr << "Usage:" << endl;
	cerr << pname << " [options] <tree file> <ancestral genome> <ancestral pan-genome> <output alignment file> <evolved sequence output file>" << endl;
	cerr << "Options:" << endl;
	cerr << "\t    --stop-codon-bias=<number> degree of bias against introducing stop codons, between 0 and 1, default 0" << endl;
	cerr << "\t    --indel-freq=<number> Frequency of indel events, relative to nucleotide substitutions" << endl;
	cerr << "\t    --small-ht-freq=<number> Frequency of small horizontal transfer and deletion events" << endl;
	cerr << "\t    --large-ht-freq=<number> Frequency of large horizontal transfer and deletion events" << endl;
	cerr << "\t    --inversion-freq=<number> Frequency of inversion events" << endl;
	cerr << "\t    --indel-size=<number> Poisson parameter for size of indels [3]" << endl;
	cerr << "\t    --small-ht-size=<number> Average length of small horizontal transfer events (sampled from an exponential distribution) [200]" << endl;
	cerr << "\t    --large-ht-min=<number> Minimum length of large horizontal transfer events (sampled uniformly between min and max) [10000]" << endl;
	cerr << "\t    --large-ht-max=<number> Maximum length of large horizontal transfer events (sampled uniformly between min and max) [60000]" << endl;
	cerr << "\t    --inversion-size=<number> Average size of inversion events, samples will be taken from an exponential distribution [50000]" << endl;
	cerr << "\t    --random-seed=<number> Use the specified random seed instead of time(NULL)" << endl;
	cerr << endl;
}

bool is_stop_codon(const vector< vector<int> >& cds_table, const string& seq, gnSeqI site){
	for(int i=0; i<cds_table.size(); i++){
		if(site < cds_table[i][0] || site > cds_table[i][1])
			continue;	// not intersecting
		if(cds_table[i][2] == 1){
			int codon_start = ((site - cds_table[i][0]) / 3) * 3;
			string codon = seq.substr(cds_table[i][0] + codon_start - 1, 3);
			if(boost::iequals(codon, "TGA") || boost::iequals(codon, "TAA") || boost::iequals(codon, "TAG") ){
				return true;
			}
		}else{
			int left = cds_table[i][0];
			int right = cds_table[i][1];
			int codon_start = ((cds_table[i][1] - site) / 3) * 3;
			string codon = seq.substr(cds_table[i][1] - codon_start - 3, 3);
			if(boost::iequals(codon, "TCA") || boost::iequals(codon, "TTA") || boost::iequals(codon, "CTA") ){
				return true;
			}
		}

	}
	return false;
}

gnAlignedSequences evolveNucleotides(PhyloTree<TreeNode>& tree, gnSequence& gns, const vector< vector<int> >& cds_table, double stop_codon_bias){
	// evolve under F81 model -- e.g. maintain background composition 
	// calculate background nucleotide composition
	vector<double> acgt(4,0);
	string seq = gns.ToString();
	for(int i=0; i<seq.length(); i++){
		switch(seq[i]){
			case 'A':
			case 'a':
				acgt[0]++;
				break;
			case 'C':
			case 'c':
				acgt[1]++;
				break;
			case 'G':
			case 'g':
				acgt[2]++;
				break;
			case 'T':
			case 't':
				acgt[3]++;
				break;
		}
	}
	double total = acgt[0] + acgt[1] + acgt[2] + acgt[3];
	for(int i=0; i<4; i++){
		acgt[i] /= total;
	}

	// walk the tree and evolve sequences
	stack<node_id_t> node_stack;
	vector<string> evolved_sequences(tree.size(), gns.ToString());
	vector<string> names(tree.size());
	int ancestor_count = 0;
	node_stack.push(tree.root);
	while(!node_stack.empty()) {
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		// add the children to the stack
		for(int i = 0; i < tree[cur_node].children.size(); i++){
			node_stack.push(tree[cur_node].children[i]);
		}
		if(tree[cur_node].name.size() > 0){
			names[cur_node] = tree[cur_node].name;
		}else{
			stringstream ss;
			ss << "Ancestor " << ancestor_count++;
			names[cur_node] = ss.str();
		}
		if(tree[cur_node].parents.size() == 0)
			continue;


		// sample the number of mutation events from the Poisson distribution
		int event_count = poissonSample(gns.length() * tree[cur_node].distance);
		// set the sequence to be identical to the parent
		evolved_sequences[cur_node] = evolved_sequences[tree[cur_node].parents[0]];
		// now for each mutation, sample its location uniformly and apply it
		int new_nt[4] = {'A','C','G','T'};
		for(int i = 0; i < event_count; i++){
			int site = uniformSample(0, gns.length());
			int nt = categoricalSample(acgt);
			// TODO: add check for stop codons here!!
			if(is_stop_codon(cds_table, evolved_sequences[cur_node], site+1)){
				// don't touch an existing stop codon...for now.
			}else{
				char cur_nt = evolved_sequences[cur_node][site];
				evolved_sequences[cur_node][site] = new_nt[nt];
				if(is_stop_codon(cds_table, evolved_sequences[cur_node], site+1)){
					// accept the new codon with probability 1 - stop_codon_bias
					double r = rndu();
					if( r > 1 - stop_codon_bias )
						evolved_sequences[cur_node][site] = cur_nt;
				}
			}
		}
	}
	gnAlignedSequences gnas;
	gnas.sequences = evolved_sequences;
	gnas.names = names;
	return gnas;
}

void parse_gff(istream& gff_in, vector< vector<int> >& cds_table){
	string line;
	while(getline(gff_in, line)){
		if(line[0] == '#')
			continue;
		istringstream line_str(line);
		vector<string> keys(9);
		for(int i=0; i<9; i++){
			line_str >> keys[i];
		}
		if(keys[2] != "CDS")
			continue;
		cds_table.push_back(vector<int>(3));
		cds_table.back()[0] = atoi(keys[3].c_str());
		cds_table.back()[1] = atoi(keys[4].c_str());
		cds_table.back()[2] = keys[6] == "+" ? 1 : 0;
	}
}

#define NELEMS(a) ( sizeof( a ) / sizeof( *a ) )

/**
 * 
 */
int main( int argc, char* argv[] ){

	const char* m_argv[] = {
		"sgEvolver",

		"--indel-freq=0.01",
		"--small-ht-freq=0",
		"--large-ht-freq=0",
		"--inversion-freq=.00001",
		"--large-ht-min=500",
		"--large-ht-max=1000",
		"--inversion-size=4000",
		"--stop-codon-bias=0.95",
		"--ancestral-gff=NC_000915.gff",
		"--accessory-gff=NC_000921.gff",
		"--random-seed=42",

		"--debug-checking",
		"test.tree",
		"NC_000915.gbk",
		"NC_000921.gbk",
		"evolved.dat",
		"evolved.fas",
	};

	int m_argc = NELEMS( m_argv );
//	argc = m_argc; argv = (char**)m_argv;

try{
	if( argc <= 0 ){
		print_usage( "sgEvolver" );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0] );
		return -1;
	}
	
	
	double indel_frequency = 0;
	double small_ht_frequency = 0;
	double large_ht_frequency = 0;
	double inversion_frequency = 0;
	gnSeqI indel_size = 1;
	gnSeqI small_ht_size = 200;
	gnSeqI large_ht_min = 10000;
	gnSeqI large_ht_max = 60000;
	gnSeqI inversion_size = 50000;
	string tree_filename;
	string ancestral_filename;
	string ancestral_pan_filename;
	string ancestral_gff_filename;
	string ancestral_pan_gff_filename;
	string output_filename;
	string evolved_filename;
	bool print_version = false;
	int64 random_seed = -1;
	double stop_codon_bias = 0; // 0 is no bias against stop codons, 1 is total avoidance of introduced stops
	char* tmp;

	// parse command line with gnu getopt
	int opt;
	int long_opt;
	int ac = argc;
	char** av = argv;
	
	const char* short_args= "";
	struct option long_opts[] = {
		{"indel-freq", required_argument, &long_opt, 1 },
		{"small-ht-freq", required_argument, &long_opt, 2 },
		{"large-ht-freq", required_argument, &long_opt, 3 },
		{"inversion-freq", required_argument, &long_opt, 4 },
		{"indel-size", required_argument, &long_opt, 5 },
		{"small-ht-size", required_argument, &long_opt, 6 },
		{"large-ht-min", required_argument, &long_opt, 7 },
		{"large-ht-max", required_argument, &long_opt, 8 },
		{"inversion-size", required_argument, &long_opt, 9 },
		{"version", no_argument, &long_opt, 10 },
		{"debug-checking", optional_argument, &long_opt, 11 },
		{"random-seed", optional_argument, &long_opt, 12 },
		{"stop-codon-bias", optional_argument, &long_opt, 13 },
		{"ancestral-gff", optional_argument, &long_opt, 14 },
		{"accessory-gff", optional_argument, &long_opt, 15 },
		{0, 0, 0, 0}	// for correct termination of option list
						// getopt_long can segfault without this
	};

	int indexptr;
	while( (opt = getopt_long( ac, av, short_args, long_opts, &indexptr )) != EOF ){
		switch( opt ){
			case 0:
				switch( long_opt ){
					case 1:
						indel_frequency = strtod( optarg, &tmp );
						break;
					case 2:
						small_ht_frequency = strtod( optarg, &tmp );
						break;
					case 3:
						large_ht_frequency = strtod( optarg, &tmp );
						break;
					case 4:
						inversion_frequency = strtod( optarg, &tmp );
						break;
					case 5:
						indel_size = strtoll( optarg, &tmp, 10 );
						break;
					case 6:
						small_ht_size = strtoll( optarg, &tmp, 10 );
						break;
					case 7:
						large_ht_min = strtoll( optarg, &tmp, 10 );
						break;
					case 8:
						large_ht_max = strtoll( optarg, &tmp, 10 );
						break;
					case 9:
						inversion_size = strtoll( optarg, &tmp, 10 );
						break;
					case 10:
						print_version = true;
						break;
					case 11:
						if( optarg != NULL )
							debugCheckingLevel( atoi( optarg ) );
						else
							debugCheckingLevel( 1 );
						break;
					case 12:
						random_seed = atol( optarg );
						break;
					case 13:
						stop_codon_bias = atof( optarg );
						break;
					case 14:
						ancestral_gff_filename = optarg;
						break;
					case 15:
						ancestral_pan_gff_filename = optarg;
						break;
					default:
						print_usage( argv[0] );
						return -1;
						break;
				}
		}
	}
	
	// print the version if the user requested it
	if( print_version ){
		cerr << "sgEvolver " << " build date " << __DATE__ << " at " << __TIME__ << endl;
	}

	if( optind + 4 >= argc ){
		cerr << "You must specify a tree file name, an input ancestral genome file name, an input ancestral pan-genome file name, an output alignment file name, and an evolved sequences output file name\n";
		return -1;
	}
	tree_filename = av[ optind ];
	ancestral_filename = av[ optind + 1 ];
	ancestral_pan_filename = av[ optind + 2 ];
	output_filename = av[ optind + 3 ];
	evolved_filename = av[ optind + 4 ];
	
	// open the tree file
	ifstream tree_file( tree_filename.c_str() );
	if( !tree_file.is_open() ){
		cerr << "Unable to open " << tree_filename << endl;
		return -1;
	}

	// parse the tree file
	PhyloTree<TreeNode> tree;
	tree.readTree( tree_file );

	// parse the ancestral (pan-)genome files
	gnSequence ancestor;
	ancestor.LoadSource(ancestral_filename);
	gnSequence ancestor_pangenome;
	ancestor_pangenome.LoadSource(ancestral_pan_filename);

	// try to load GFF annotations, if available
	vector< vector<int> > ancestral_cds;
	vector< vector<int> > ancestral_pan_cds;
	if(ancestral_gff_filename.length() > 0){
		ifstream gff_in(ancestral_gff_filename.c_str());
		if( !gff_in.is_open() ){
			cerr << "Unable to open " << ancestral_gff_filename << endl;
			return -1;
		}
		parse_gff(gff_in, ancestral_cds);
	}
	if(ancestral_pan_gff_filename.length() > 0){
		ifstream gff_in(ancestral_gff_filename.c_str());
		if( !gff_in.is_open() ){
			cerr << "Unable to open " << ancestral_pan_gff_filename << endl;
			return -1;
		}
		parse_gff(gff_in, ancestral_pan_cds);
	}

	// do some nucleotide substitution
	gnAlignedSequences destination;
	gnAlignedSequences donor;
	destination = evolveNucleotides(tree, ancestor, ancestral_cds, stop_codon_bias);
	donor = evolveNucleotides(tree, ancestor_pangenome, ancestral_pan_cds, stop_codon_bias);

	string line;

	// do the benjarath length consistency check
	for( int seqI = 1; seqI < destination.sequences.size(); seqI++ ){
			cerr << "dest seq " << seqI -1 << " length " << destination.sequences[ seqI-1 ].length() << endl;
			cerr << "dest seq " << seqI << " length " << destination.sequences[ seqI ].length() << endl;
		if( destination.sequences[ seqI ].length() != destination.sequences[ seqI - 1 ].length() )
		{
			throw "Error parsing input alignment.  Check that all sequences in the alignment have equal lengths and that the correct number of sequences is specified.";
		}
	}
        for( int seqI = 1; seqI < destination.sequences.size(); seqI++ ){
			cerr << "donor seq " << seqI -1 << " length " << donor.sequences[ seqI-1 ].length() << endl;
			cerr << "donor seq " << seqI << " length " << donor.sequences[ seqI ].length() << endl;
		if( donor.sequences[ seqI ].length() != donor.sequences[ seqI - 1 ].length() )
		{
			throw "Error parsing input alignment.  Check that all sequences in the alignment have equal lengths and that the correct number of sequences is specified.";
		}
        }

	// initialize the mutators with the parameters given on the command line
	IndelInserter ii( donor, indel_size );
	IndelDeleter id( donor, indel_size );
	SmallHTInserter shti( donor, small_ht_size );
	SmallHTDeleter shtd( donor, small_ht_size );
	LargeHTInserter lhti( donor, large_ht_min, large_ht_max );
	LargeHTDeleter lhtd( donor, large_ht_min, large_ht_max );
	Inverter inv( donor, inversion_size );

	vector< MutationEvent > mutation_events;
	mutation_events.push_back( MutationEvent( "indel_insert", indel_frequency, ii ) );
	mutation_events.push_back( MutationEvent( "indel_deletion", indel_frequency, id ) );
	mutation_events.push_back( MutationEvent( "small_ht_insert", small_ht_frequency, shti ) );
	mutation_events.push_back( MutationEvent( "small_ht_deletion", small_ht_frequency, shtd ) );
	mutation_events.push_back( MutationEvent( "large_ht_insert", large_ht_frequency, lhti ) );
	mutation_events.push_back( MutationEvent( "large_ht_deletion", large_ht_frequency, lhtd ) );
	mutation_events.push_back( MutationEvent( "inversion", inversion_frequency, inv ) );
	
	// start evolving the sequence alignment
	AlignmentEvolver ae( mutation_events );
	
	Alignment evolved_alignment( destination, donor );
	ae.evolve( tree, evolved_alignment, random_seed );
	
	// open the alignment output file and write the alignment
	ofstream output_file( output_filename.c_str() );
	if( !output_file.is_open() ){
		cerr << "Error opening " << output_filename << endl;
		return -1;
	}
	evolved_alignment.writeAlignment( output_file, tree );
	
	// open the evolved sequences file and write the sequences
	ofstream evolved_file( evolved_filename.c_str() );
	if( !evolved_file.is_open() ){
		cerr << "Error opening " << evolved_filename << endl;
		return -1;
	}
	evolved_alignment.writeEvolvedSequences( evolved_file, tree );

}catch( gnException& gne ) {
	cerr << "Unhandled gnException: " << gne << endl;
	return -1;
}catch( exception& e ) {
	cerr << "Unhandled exception: " << e.what() << endl;
	return -1;
}catch( char const* message ) {
	cerr << message << endl;
	return -1;
}

//catch(...){
//	cerr << "Unknown exception occurred.\n";
//	return -1;
//}

	return 0;
}
