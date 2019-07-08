#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Mutator.h"
#include <vector>
#include <cmath>

extern "C"{
#include "twister.h"
}

using namespace std;
using namespace genome;
using namespace mems;

/**
 * Sample an integer from a uniform distribution bounded by min and max
 * @param min 	The minimum sample value
 * @param max 	The maximum sample value
 */
gnSeqI uniformSample( gnSeqI min, gnSeqI max ){
	if( max == min )
		return 0;

	gnSeqI sample = (gnSeqI)(rndu() * (max - min));
	return sample + min;
}

/**
 * Sample an exponential distribution with mean theta
 * @param theta	the exponential mean
 */
double exponentialSample( double theta ){
	double z_samp = rndu();
	z_samp = -log( z_samp ) * theta;
	return z_samp;
}

/**
 * Sample a non-negative integer from a poisson distribution with mean p
 * @param p	the poisson mean
 */
gnSeqI poissonSample( double p ){
	gnSeqI count = 0;
	double lambda = 1;
	double sum = 0;
	for( ;; count++ ){
		sum += exponentialSample( lambda );
		if( sum >= p )
			break;
	}
	return count;
}

/**
 * Sample a from the categorical distribution
 * @param d the categorical distribution
 */
int categoricalSample( const vector<double>& d ){
	int count = 0;
	double sum = 0;
	double r = rndu();
	for( ; count < d.size(); count++ ){
		sum += d[count];
		if( sum >= r )
			break;
	}
	return count;
}

double prob_inframe = 0.95;
gnSeqI getIndelSize( int size ){
	gnSeqI len = 0;
	do{
		len = poissonSample( size );
		len *= 3;
	}while( len == 0 );
	// with some probability, create a potential frameshift
	if(rndu() > prob_inframe){
		len -= rndu() > 0.5 ? 1 : 2;
	}
	return len;
}

void IndelInserter::getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len ) {
	source_len = getIndelSize(size);
	source_len = source_len < donor.alignedSeqsSize() ? source_len : donor.alignedSeqsSize() - 1;
	source_start = uniformSample( 0, donor.alignedSeqsSize() - source_len - 1 );
	dest = uniformSample( 0, dest_len );
}

void SmallHTInserter::getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len ) {
	do{
		source_len = (gnSeqI)exponentialSample( size ) + 1;	// force it between 1 and sequence length
	}while( source_len == 0 );
	source_len = source_len < donor.alignedSeqsSize() ? source_len : donor.alignedSeqsSize() - 1;
	source_start = uniformSample( 0, donor.alignedSeqsSize() - source_len - 1 );
	dest = uniformSample( 0, dest_len );
}

void LargeHTInserter::getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len ) {
	source_len = uniformSample( min_size, max_size );
	source_len = source_len < donor.alignedSeqsSize() ? source_len : donor.alignedSeqsSize() - 1;
	source_len = source_len == 0 ? 1 : source_len;
	source_start = uniformSample( 0, donor.alignedSeqsSize() - source_len - 1 );
	dest = uniformSample( 0, dest_len );
}

void IndelDeleter::getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len ) {
	len = getIndelSize(size);
	len = len < dest_len ? len : dest_len;
	start = uniformSample( 0, dest_len - len );
}

void SmallHTDeleter::getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len ) {
	do{
		len = (gnSeqI)exponentialSample( size );	// don't let it be 0
	}while( len == 0 );
	len = len < dest_len ? len : dest_len;
	start = uniformSample( 0, dest_len - len );
}

void LargeHTDeleter::getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len ) {
	len = uniformSample( min_size, max_size );
	len = len < dest_len ? len : dest_len;
	start = uniformSample( 0, dest_len - len );
}

void Inverter::getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len ) {
	do{
		len = (gnSeqI)exponentialSample( size );	// don't let it be 0
	}while( len == 0 );
	// treat as circular
	gnSeqI s = uniformSample( 0, dest_len );
	gnSeqI circ = (s+len)%dest_len;
	start = s < circ ? s : circ;
	len = abs((int64)circ-(int64)s);
}


// take sequence from donor and insert it into sequences below nodeI
void Inserter::mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment ) {
	gnSeqI source_start;
	gnSeqI source_length;
	gnSeqI dest, dest_col;
	// nodeI and all sequences below it in the tree should have the same length
	gnSeqI dest_len = evolved_alignment.sequenceLength( nodeI );
	getLocation( source_start, source_length, dest, dest_len );
	dest_col = evolved_alignment.getColumnIndex( nodeI, dest );
	recursiveInsert( tree.root, nodeI, tree, evolved_alignment, dest_col, source_start, source_length );
	if( debugChecks() )
		evolved_alignment.checkLengths();
}

void Inserter::recursiveInsert( node_id_t cur_node, node_id_t insert_node, const PhyloTree<TreeNode>& t, 
		Alignment& evolved_alignment, gnSeqI point, gnSeqI source_left, gnSeqI source_length ){
	if( cur_node == insert_node ){
		// insert the actual sequence in this one
		evolved_alignment.applyInsertion( cur_node, point, source_left, source_length );
	}else
		evolved_alignment.applyGapInsertion( cur_node, point, source_length );
		
	for( uint childI = 0; childI < t[ cur_node ].children.size(); childI++ ){
		// make the insert node follow through the tree down from the initial insertion
		node_id_t new_ins_node = cur_node == insert_node ? t[ cur_node ].children[ childI ] : insert_node;
		recursiveInsert( t[ cur_node ].children[ childI ], new_ins_node, t, evolved_alignment, point, source_left, source_length );
	}
}

void Deleter::mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment ) {
	gnSeqI start;
	gnSeqI length;
	// nodeI and all sequences below it in the tree should have the same length
	gnSeqI dest_len = evolved_alignment.sequenceLength( nodeI );
	getLocation( start, length, dest_len );
	if( length == 0 )
		return;	// there isn't anything left to delete

	// convert sequence position to alignment column
	evolved_alignment.getDeletionColumns( nodeI, start, length );
	recursiveDelete( tree.root, nodeI, tree, evolved_alignment, start, length );
	if( debugChecks() )
		evolved_alignment.checkLengths();
}

void Deleter::recursiveDelete( node_id_t cur_node, node_id_t insert_node, const PhyloTree<TreeNode>& t, Alignment& evolved_alignment, 
		gnSeqI left_end, gnSeqI length ){
	if( cur_node == insert_node )
		// delete sequence in this one
		evolved_alignment.applyDeletion( cur_node, left_end, length );
		
	for( uint childI = 0; childI < t[ cur_node ].children.size(); childI++ ){
		// make the insert node follow through the tree down from the initial insertion
		node_id_t new_ins_node = cur_node == insert_node ? t[ cur_node ].children[ childI ] : insert_node;
		recursiveDelete( t[ cur_node ].children[ childI ], new_ins_node, t, evolved_alignment, left_end, length );
	}
}

void Inverter::mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment ) {
	gnSeqI start;
	gnSeqI length;
	// nodeI and all sequences below it in the tree should have the same length
	gnSeqI dest_len = evolved_alignment.sequenceLength( nodeI );
	getLocation( start, length, dest_len );
	if( length == 0 )
		return;	// there isn't anything left to invert

	// convert sequence position to alignment column
	evolved_alignment.getDeletionColumns( nodeI, start, length );
	evolved_alignment.addInversion( nodeI, start, length );
	cerr << "Inversion\t" << start << "\t" << start+length << "\n";
// don't do recursive inversion -- applyInversions() will take care of applying it
// to each sequence	
}
