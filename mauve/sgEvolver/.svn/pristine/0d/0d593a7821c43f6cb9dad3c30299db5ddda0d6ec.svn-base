#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __Mutator_h__
#define __Mutator_h__

#include "libMems/gnAlignedSequences.h"
#include "libMems/PhyloTree.h"
#include <vector>
#include "gnStringSource.h"
#include "libGenome/gnSequence.h"
#include "Alignment.h"

/**
 * The Mutator class defines a generic interface for classes that
 * perform mutations on a sequence alignment
 */
class Mutator {
public:
	/**
	 * The Mutator constructor specifies a sequence alignment from which to draw insertion sequences
	 * @param donor_seqs Sequence alignment for insertion sequences
	 */
	Mutator( mems::gnAlignedSequences& donor_seqs );
	virtual ~Mutator(){}
	/**
	 * Perform a mutation on a sequence alignment.
	 * Given a sequence alignment and a phylogenetic tree relating the sequences in
	 * the alignment, this function will perform a mutation on the sequence alignment
	 * at a particular point in the alignment's phylogenetic history.
	 * @param nodeI		The node in the phylogenetic tree to perform the mutation at
	 * @param tree		The phylogenetic tree corresponding to the sequence alignment
	 * @param evolved_alignment		An alignment of sequences corresponding to each node
	 * 								(including ancestral nodes) of the phylogenetic tree.
	 *								If there are multiple entries in the vector, the sequence
	 *								is not collinear.
	 */
	virtual void mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment ) = 0;
protected:
	mems::gnAlignedSequences& donor;
	std::vector< genome::gnSequence > donor_table;
private:
	Mutator();
};

gnSeqI poissonSample( double p );
double exponentialSample( double theta );
gnSeqI uniformSample( gnSeqI min, gnSeqI max );
int categoricalSample( const std::vector<double>& d );

inline
Mutator::Mutator( mems::gnAlignedSequences& donor_seqs ) : 
donor( donor_seqs )
{
	for( uint seqI = 0; seqI < donor.sequences.size(); seqI++ ){
		gnStringSource* gnss = new gnStringSource( donor.sequences[ seqI ] );
		genome::gnSequence gns( *gnss->GetSpec() );
		donor_table.push_back( gns );
	}
}

/**
 * The Inserter class provides basic functionality for classes
 * that perform insertion mutations
 */
class Inserter : public Mutator {
public:
	Inserter( mems::gnAlignedSequences& donor_seqs );
	void mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment );
protected:
	void recursiveInsert( node_id_t cur_node, node_id_t insert_node, const PhyloTree<TreeNode>& t, 
		Alignment& evolved_alignment, gnSeqI point, gnSeqI source_left, gnSeqI source_length );
	virtual void getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len ) = 0;
};

inline
Inserter::Inserter( mems::gnAlignedSequences& donor_seqs ) : 
Mutator( donor_seqs )
{}

/**
 * The indel inserter samples insertion locations uniformly and insertion lengths from
 * a poisson distribution
 */
class IndelInserter: public Inserter {
public:
	IndelInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI size );
	IndelInserter( const IndelInserter& id );
	IndelInserter& operator=( const IndelInserter& id );

protected:
	void getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len );
	gnSeqI size;
};

/**
 * The small horizontal transfer (HT) inserter samples insertion locations uniformly and 
 * insertion lengths from an exponential distribution
 */
class SmallHTInserter: public Inserter {
public:
	SmallHTInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI size );
	SmallHTInserter( const SmallHTInserter& id );
	SmallHTInserter& operator=( const SmallHTInserter& id );

protected:
	void getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len );
	gnSeqI size;
};

/**
 * The small horizontal transfer (HT) inserter samples insertion locations uniformly and 
 * insertion lengths from a bounded uniform distribution
 */
class LargeHTInserter: public Inserter {
public:
	LargeHTInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI min_size, gnSeqI max_size );
	LargeHTInserter( const LargeHTInserter& id );
	LargeHTInserter& operator=( const LargeHTInserter& id );

protected:
	void getLocation( gnSeqI& source_start, gnSeqI& source_len, gnSeqI& dest, gnSeqI dest_len );
	gnSeqI min_size;
	gnSeqI max_size;
};

/**
 * The Deleter class provides basic functionality for classes
 * that perform deletion mutations
 */
class Deleter : public Mutator {
public:
	Deleter( mems::gnAlignedSequences& donor_seqs );
	void mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment );
protected:
	void recursiveDelete( node_id_t cur_node, node_id_t insert_node, const PhyloTree<TreeNode>& t,
		Alignment& evolved_alignment, gnSeqI left_end, gnSeqI len );
	virtual void getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len ) = 0;
};

inline
Deleter::Deleter( mems::gnAlignedSequences& donor_seqs ) : 
Mutator( donor_seqs )
{}

/**
 * The indel deleter samples deletion locations uniformly and deletion lengths from
 * a poisson distribution
 */
class IndelDeleter: public Deleter {
public:
	IndelDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI size );
	IndelDeleter( const IndelDeleter& id );
	IndelDeleter& operator=( const IndelDeleter& id );

protected:
	void getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len );
	gnSeqI size;
};

/**
 * The small Horizontal Transfer (HT) deleter samples deletion locations uniformly 
 * and deletion lengths from an exponential distribution
 */
class SmallHTDeleter: public Deleter {
public:
	SmallHTDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI size );
	SmallHTDeleter( const SmallHTDeleter& id );
	SmallHTDeleter& operator=( const SmallHTDeleter& id );

protected:
	void getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len );
	gnSeqI size;
};

/**
 * The small Horizontal Transfer (HT) deleter samples deletion locations uniformly 
 * and deletion lengths from a bounded uniform distribution
 */
class LargeHTDeleter: public Deleter {
public:
	LargeHTDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI min_size, gnSeqI max_size );
	LargeHTDeleter( const LargeHTDeleter& id );
	LargeHTDeleter& operator=( const LargeHTDeleter& id );

protected:
	void getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len );
	gnSeqI min_size;
	gnSeqI max_size;
};

/**
 * The Inverter class mutates sequences by performing inversions.
 * The location of inversion is chosen uniformly along the sequence,
 * and the length of inversion is chosen according to an exponential distribution
 * with a mean length specified by the 'size' constructor parameter.
 */
class Inverter : public Mutator {
public:
	Inverter( mems::gnAlignedSequences& donor_seqs, gnSeqI size );
	Inverter( const Inverter& i );
	Inverter& operator=( const Inverter& i );

	void mutate( node_id_t nodeI, const PhyloTree<TreeNode>& tree, Alignment& evolved_alignment );
protected:
	void recursiveInvert( node_id_t cur_node, node_id_t insert_node, const PhyloTree<TreeNode>& t,
		Alignment& evolved_alignment, gnSeqI left_end, gnSeqI length );
	void getLocation( gnSeqI& start, gnSeqI& len, gnSeqI dest_len );
	gnSeqI size;
};

inline
IndelInserter::IndelInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI size ) : 
Inserter( donor_seqs )
{
	this->size = size;
}
inline
IndelInserter::IndelInserter( const IndelInserter& i ) :
Inserter( i.donor ) 
{
	*this = i;
}
inline
IndelInserter& IndelInserter::operator=( const IndelInserter& i ) {
	size = i.size;
	return *this;
}

inline
IndelDeleter::IndelDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI size ) : 
Deleter( donor_seqs )
{
	this->size = size;
}
inline
IndelDeleter::IndelDeleter( const IndelDeleter& i ) :
Deleter( i.donor ) 
{
	*this = i;
}
inline
IndelDeleter& IndelDeleter::operator=( const IndelDeleter& i ) {
	size = i.size;
	return *this;
}

inline
SmallHTInserter::SmallHTInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI size ) : 
Inserter( donor_seqs )
{
	this->size = size;
}
inline
SmallHTInserter::SmallHTInserter( const SmallHTInserter& i ) :
Inserter( i.donor ) 
{
	*this = i;
}
inline
SmallHTInserter& SmallHTInserter::operator=( const SmallHTInserter& i ) {
	size = i.size;
	return *this;
}

inline
SmallHTDeleter::SmallHTDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI size ) : 
Deleter( donor_seqs )
{
	this->size = size;
}
inline
SmallHTDeleter::SmallHTDeleter( const SmallHTDeleter& i ) :
Deleter( i.donor ) 
{
	*this = i;
}
inline
SmallHTDeleter& SmallHTDeleter::operator=( const SmallHTDeleter& i ) {
	size = i.size;
	return *this;
}

inline
LargeHTInserter::LargeHTInserter( mems::gnAlignedSequences& donor_seqs, gnSeqI min_size, gnSeqI max_size ) : 
Inserter( donor_seqs )
{
	this->min_size = min_size;
	this->max_size = max_size;
}
inline
LargeHTInserter::LargeHTInserter( const LargeHTInserter& i ) :
Inserter( i.donor ) 
{
	*this = i;
}
inline
LargeHTInserter& LargeHTInserter::operator=( const LargeHTInserter& i ) {
	min_size = i.min_size;
	max_size = i.max_size;
	return *this;
}

inline
LargeHTDeleter::LargeHTDeleter( mems::gnAlignedSequences& donor_seqs, gnSeqI min_size, gnSeqI max_size ) : 
Deleter( donor_seqs )
{
	this->min_size = min_size;
	this->max_size = max_size;
}
inline
LargeHTDeleter::LargeHTDeleter( const LargeHTDeleter& i ) :
Deleter( i.donor ) 
{
	*this = i;
}
inline
LargeHTDeleter& LargeHTDeleter::operator=( const LargeHTDeleter& i ) {
	min_size = i.min_size;
	max_size = i.max_size;
	return *this;
}

inline
Inverter::Inverter( mems::gnAlignedSequences& donor_seqs, gnSeqI size ) : 
Mutator( donor_seqs )
{
	this->size = size;
}
inline
Inverter::Inverter( const Inverter& i ) :
Mutator( i.donor ) 
{
	*this = i;
}
inline
Inverter& Inverter::operator=( const Inverter& i ) {
	size = i.size;
	return *this;
}


#endif // __Mutator_h__
