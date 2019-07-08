#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __Alignment_h__
#define __Alignment_h__

#include "libMems/gnAlignedSequences.h"
#include "libMems/Aligner.h"
#include "libGenome/gnSequence.h"
#include <vector>
//#include "GappedIntervalSequence.h"
#include "GisSplayTree.h"
#include "libMems/PhyloTree.h"

#define GappedIntervalSequence GisSplayTree

/** 
 *  Sets or checks the level of extra debug checking to perform.
 *  0 indicates no checking
 *  1 indicates intermittent checks
 *  2 specifies data structure consistency checks after every mutation operation
 *  This can be set using the --debug-checking=# command line switch
 *  Setting this to 2 will make the program _very_ slow
 */
int debugCheckingLevel( int level_request = -1 );

/** 
 * An auxiliary function to help debug checking.  
 */
bool debugChecks( int set = -1 );

/**
 * This class describes an interval of a mutated sequence
 * The interval of sequence is either a range of sequence from another source
 * or is a string of gaps
 * The GappedIntervalSequence will be a search tree of these that tracks sequence
 * coordinates through insertions and deletions, allowing fast lookup of 
 * ungapped sequence coordinates in a gapped sequence
 */
class SourceInterval {
public:
	SourceInterval( uint sourceI, int64 leftI, gnSeqI len );
	SourceInterval( gnSeqI len );
	SourceInterval( const SourceInterval& ti );
	
	gnSeqI GetLength() const;
	void CropStart( gnSeqI crop_len );
	void CropEnd( gnSeqI crop_len );
	gnSeqI GetSeqLength() const;
	
	int64 leftEnd() const;
	bool isGap() const;
	uint sourceID() const;
private:
	uint source_id;		/**< An ID pointing back to the sequence this interval refers to */
	int64 left_end;		/**< The left end of this interval in the source sequence */
	gnSeqI length;		/**< The length of this interval */
	bool gap;			/**< True if this is an interval of gaps */
};

/** 
 * This class is used to track relationships between LCBs during the inversion process.
 */
//class LCB{
//public:
//	vector< int64 > left_end;	/**< The left end position of the LCB in each sequence */
//	vector< int64 > right_end;  /**< The right end position of the LCB in each sequence */
//	vector< uint > left_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
//	vector< uint > right_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
//	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
//};

/**
 * Compares left_end coordinate of two LCBs
 */
class LCBLeftComparator {
public:
	LCBLeftComparator( uint seq ){
		m_seq = seq;
	}
	LCBLeftComparator( LCBLeftComparator& lmc ){
		m_seq = lmc.m_seq;
	}
	boolean operator()(const mems::LCB& a, const mems::LCB& b) const{
		return genome::absolut( a.left_end[ m_seq ] ) < genome::absolut( b.left_end[ m_seq ] );
	}
	~LCBLeftComparator() {
	}
protected:
	uint m_seq;
private:
	LCBLeftComparator();
};


class Alignment {
public:
	Alignment( const mems::gnAlignedSequences& base_alignment, const mems::gnAlignedSequences& donor_seqs );
	Alignment( const Alignment& a );
	Alignment& operator=( const Alignment& a );
	
	gnSeqI sequenceLength( uint seqI ) const;
	/**
	 * Translates the sequence coordinate baseI to a contig local coordinate
	 */
	void getBase( uint seqI, gnSeqI& baseI, uint& alignI );
	void addInversion( uint seq, gnSeqI left_end, gnSeqI length );
	void applyDeletion( uint seq, gnSeqI left_end, gnSeqI length );
	void applyInsertion( uint seq, gnSeqI point, gnSeqI source_left, gnSeqI source_length );
	void applyGapInsertion( uint seq, gnSeqI point, gnSeqI length );
	void applyInversions( const PhyloTree<TreeNode>& t );
	gnSeqI getColumnIndex( gnSeqI seqI, gnSeqI baseI );
	void getDeletionColumns( uint seq, gnSeqI& left_end, gnSeqI& length );

	/**
	 * Write out the alignment in XMFA (extended Multi-FastA) format
	 * @param os	The output stream to write the alignment
	 * @param t		The phylogenetic tree for the aligned sequences
	 * @param write_ancestors	If set to true, ancestral sequences will be written.  False by default.
	 */
	void writeAlignment( std::ostream& os, const PhyloTree<TreeNode>& t, boolean write_ancestors = false ) const;

	/**
	 * Write the fully evolved sequences to a Multi-FastA file
	 * @param os	The output stream to write the sequences into
	 * @param t		The phylogenetic tree for the sequences
	 * @param write_ancestors	If set to true, ancestral sequences will be written.  False by default.
	 */
	void writeEvolvedSequences( std::ostream& os, const PhyloTree<TreeNode>& t, boolean write_ancestors = false ) const;
	gnSeqI length() const{ return intervals[ 0 ].length(); }
	void checkLengths() const;
protected:
	std::vector< mems::gnAlignedSequences > alignment;	/**< A vector of alignments */
	const mems::gnAlignedSequences& donor;	/**< The alignment of donor sequences used for insertions */
	/** 
	 * A GappedIntervalSequence for each sequence maps alignment coordinates to 
	 * original sequence coordinates and vice-versa 
	 */
	std::vector< GappedIntervalSequence< SourceInterval > > intervals;
	/**
	 * A list of inversion events to be applied to the sequences after all
	 * insertion and deletion events have completed
	 * inversion_list indexed by: [ sequence ][ inversionI ]< left_end, length >
	 */
	std::vector< std::list< std::pair< gnSeqI, gnSeqI > > > inversion_list;
	std::vector< mems::LCB > lcb_list;	/**< The locally collinear blocks (regions of orthologous sequence without rearrangements) */
	std::vector< uint > first_lcbs;		/**< The left-most LCB in each sequence */
	std::vector< std::string > ungapped_seqs;	/**< the fully evolved sequences */
};

CREATE_EXCEPTION( EvolutionError );

#endif	//__Alignment_h__

