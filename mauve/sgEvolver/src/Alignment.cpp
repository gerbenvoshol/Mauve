#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Alignment.h"
#include <vector>
#include "gnStringSource.h"
#include "gnGapSource.h"
#include "libGenome/gnFASSource.h"
#include <algorithm>
using namespace std;
using namespace genome;
using namespace mems;

int debugCheckingLevel( int level_request )
{
	static int debug_checking_level = 0;
	if( level_request >= 0 ){
		debug_checking_level = level_request;
		GappedIntervalSequence< SourceInterval >::debug_checking = level_request;
	}
	return debug_checking_level;
}
bool debugChecks( int set ){
	static bool debug_checks = false;	
	if( set == 0 )
		debug_checks = false;
	else if( set == 1 )
		debug_checks = true;
	return debug_checks;
}

SourceInterval::SourceInterval( uint sourceI, int64 leftI, gnSeqI len ) :
	source_id( sourceI ),
	left_end( leftI ),
	length( len ),
	gap( false )
{}

SourceInterval::SourceInterval( gnSeqI len ) :
	source_id( 0 ),
	left_end( 0 ),
	length( len ),
	gap( true )
{}

SourceInterval::SourceInterval( const SourceInterval& ti ) :
	source_id( ti.source_id ),
	left_end( ti.left_end ),
	length( ti.length ),
	gap( ti.gap )
{}

gnSeqI SourceInterval::GetLength() const{
	return length;
}

gnSeqI SourceInterval::GetSeqLength() const{
	return gap ? 0 : length;
}

bool SourceInterval::isGap() const{
	return gap;
}

int64 SourceInterval::leftEnd() const{
	return left_end;
}

uint SourceInterval::sourceID() const{
	return source_id;
}

void SourceInterval::CropStart( gnSeqI crop_len ){
	if( left_end >= 0 )
		left_end += crop_len;
	length -= crop_len;
}

void SourceInterval::CropEnd( gnSeqI crop_len ){
	if( left_end < 0 )
		left_end -= crop_len;
	length -= crop_len;
}

void Alignment::checkLengths() const{
	for( uint testI = 1; testI < intervals.size(); testI++ ){
		intervals[ testI ].checkTree();
		if( intervals[ testI ].length() != intervals[ testI - 1 ].length() ){
			cerr << "when I say oh, you say shit\n";
			cerr << "seq " << testI << " has length: " << intervals[ testI ].length() << endl;
			cerr << "but seq " << testI - 1 << " has length: " << intervals[ testI - 1 ].length() << endl;
		}
	}
}
const uint TARGET_ID = 1;
const uint DONOR_ID = 2;

Alignment::Alignment( const gnAlignedSequences& base_alignment, const gnAlignedSequences& donor_seqs ) :
	donor( donor_seqs )
 {
 	alignment.push_back( base_alignment );
	inversion_list = vector< list< pair< gnSeqI, gnSeqI > > >( base_alignment.sequences.size(), list< pair< gnSeqI, gnSeqI > >() );
	intervals = vector< GappedIntervalSequence< SourceInterval > >( base_alignment.sequences.size() );
	for( uint seqI = 0; seqI < intervals.size(); seqI++ ){
		SourceInterval si( TARGET_ID, 0, alignment[0].sequences[ seqI ].length() );
		intervals[ seqI ].insert( si );
	}
}

Alignment::Alignment( const Alignment& a ) :
	donor( a.donor )
{
	*this = a;
}

Alignment& Alignment::operator=( const Alignment& a ) {
//	donor = a.donor;
	alignment = a.alignment;
	intervals = a.intervals;
	return *this;
}

gnSeqI Alignment::sequenceLength( uint seqI ) const {
	return intervals[ seqI ].sequenceLength();
}

/**
 * Return the index of the baseI'th non-gap character
 */
gnSeqI Alignment::getColumnIndex( gnSeqI seqI, gnSeqI baseI ) {
	GappedIntervalSequence< SourceInterval >::iterator iter = intervals[ seqI ].find_seqindex( baseI );
	return intervals[ seqI ].getStart( iter ) + baseI - intervals[ seqI ].getSequenceStart( iter );
}

/**
 * Adds an inversion record to a particular sequence in the alignment
 * the left_end and length refer to alignment columns, NOT
 * sequence coordinates
 */
void Alignment::addInversion( uint seq, gnSeqI left_end, gnSeqI length ) {
	if( alignment.size() != 1 )
		Throw_gnExMsg( EvolutionError(), "Inversions were applied before evolution finished\n" );
	inversion_list[ seq ].push_back( pair< gnSeqI, gnSeqI >( left_end, length ) );
}

void Alignment::getDeletionColumns( uint seq, gnSeqI& left_end, gnSeqI& length ){
	gnSeqI new_len = left_end + length;
	length = getColumnIndex( seq, new_len );
	// determine the left column
	left_end = getColumnIndex( seq, left_end );
	
	length = length - left_end;
}

/**
 * Applies a deletion to a particular sequence in the alignment
 * the left_end and length refer to alignment columns, NOT
 * sequence coordinates
 */
void Alignment::applyDeletion( uint seq, gnSeqI left_end, gnSeqI length ) {
	if( alignment.size() != 1 )
		Throw_gnExMsg( EvolutionError(), "Inversions were applied before evolution finished\n" );
	
	// apply the deletion to the gapped interval sequence	
	intervals[ seq ].erase( left_end, length );
	
	SourceInterval si( length );
	intervals[ seq ].insert( si, left_end );
	
	// now apply the deletion to the inversion records
	std::list< pair< gnSeqI, gnSeqI > >::iterator inv_iter = inversion_list[ seq ].begin();
	boolean deleted_begin = false;
	for( ; inv_iter != inversion_list[ seq ].end(); inv_iter++ ){
		if( deleted_begin ){
			inv_iter--;
			deleted_begin = false;
		}
		gnSeqI i_left = inv_iter->first;
		gnSeqI i_right = inv_iter->first + inv_iter->second;
		gnSeqI right_end = left_end + length;
		if( i_left >= left_end && i_left < right_end ){
			if( i_right <= right_end ){
				// the inversion is completely within the deletion.  delete the inversion
				std::list< pair< gnSeqI, gnSeqI > >::iterator del_iter = inv_iter;
				if( inv_iter == inversion_list[ seq ].begin() ){
					inversion_list[ seq ].erase( del_iter );
					inv_iter = inversion_list[ seq ].begin();
					deleted_begin = true;
					if( inversion_list[ seq ].size() == 0 )
						break;	// get out if there's nothing left
				}else{
					inv_iter--;
					inversion_list[ seq ].erase( del_iter );
				}
				continue;
			}
			
			// crop the inversion coordinates
			inv_iter->second -= right_end - i_left;
			inv_iter->first = right_end;
		}else if( i_right >= left_end && i_right < right_end ){
			// shorten the inversion
			inv_iter->second -= i_right - left_end;
		}else if( i_left < left_end && i_right >= right_end ){
			// take a chunk out of the middle of the inversion
			inv_iter->second -= length;
		}
	}
}

void Alignment::applyInsertion( uint seq, gnSeqI point, gnSeqI source_left, gnSeqI source_length ) {
	if( alignment.size() != 1 )
		Throw_gnExMsg( EvolutionError(), "Inversions were applied before evolution finished\n" );

	// insert the new sequence
	SourceInterval si( DONOR_ID, source_left, source_length );
	intervals[ seq ].insert( si, point );
	
	// apply the insertion coordinates to the inversion records
	std::list< pair< gnSeqI, gnSeqI > >::iterator inv_iter = inversion_list[ seq ].begin();
	for( ; inv_iter != inversion_list[ seq ].end(); inv_iter++ ){
		gnSeqI i_left = inv_iter->first;
		gnSeqI i_right = inv_iter->first + inv_iter->second;
		if( i_left <= point && point < i_right ){
			inv_iter->second += source_length;
		}
	}
}

void Alignment::applyGapInsertion( uint seq, gnSeqI point, gnSeqI length ) {
	if( alignment.size() != 1 )
		Throw_gnExMsg( EvolutionError(), "Inversions were applied before evolution finished\n" );

	// insert the gap sequence
	SourceInterval si( length );
	intervals[ seq ].insert( si, point );
	
	// apply the insertion coordinates to the inversion records
	std::list< pair< gnSeqI, gnSeqI > >::iterator inv_iter = inversion_list[ seq ].begin();
	for( ; inv_iter != inversion_list[ seq ].end(); inv_iter++ ){
		gnSeqI i_left = inv_iter->first;
		gnSeqI i_right = inv_iter->first + inv_iter->second;
		if( i_left <= point && point < i_right ){
			inv_iter->second += length;
		}
	}
}


/**
 * Follows a phylogenetic tree to determine which nodes of the tree need to have
 * regions of sequence inverted.  Performs the inversion on the appropriate sequences
 */
void recursiveInvert( node_id_t cur_node, node_id_t invert_node, const PhyloTree<TreeNode>& t,
	vector< LCB >& lcb_list, vector< uint >& first_lcbs, uint left_lcb, uint right_lcb ){

	if( cur_node == invert_node ){
		
		// LCBs should be properly broken up.  
		// left_lcb and right_lcb are the left-most and right-most LCBs to be inverted
		// Inversion takes place here:
		// Need to reverse ordering and change sign of the LCBs
		uint cur_left = left_lcb;
		uint left_tmp = lcb_list[ cur_left ].left_adjacency[ cur_node ];
		uint right_tmp = lcb_list[ cur_left ].right_adjacency[ cur_node ];
		while( cur_left != right_lcb ){
			uint adj = right_tmp;
			right_tmp = lcb_list[ adj ].right_adjacency[ cur_node ];
			lcb_list[ adj ].right_adjacency[ cur_node ] = cur_left;
			lcb_list[ cur_left ].left_adjacency[ cur_node ] = adj;
			// change sign of cur_left
			lcb_list[ cur_left ].left_end[ cur_node ] *= -1;
			lcb_list[ cur_left ].right_end[ cur_node ] *= -1;
			// advance to the next lcb
			cur_left = adj;
		}

		// change sign of right_lcb
		lcb_list[ right_lcb ].left_end[ cur_node ] *= -1;
		lcb_list[ right_lcb ].right_end[ cur_node ] *= -1;

		// re-link the ends
		// right_tmp needs to point to left_lcb
		lcb_list[ left_lcb ].right_adjacency[ cur_node ] = right_tmp;
		if( right_tmp != (uint)-1 )
			lcb_list[ right_tmp ].left_adjacency[ cur_node ] = left_lcb;
		

		// left_tmp needs to point to right_lcb
		lcb_list[ right_lcb ].left_adjacency[ cur_node ] = left_tmp;
		if( left_tmp != (uint)-1 ){
			lcb_list[ left_tmp ].right_adjacency[ cur_node ] = right_lcb;
		}else{
			// update first_lcbs
			first_lcbs[ cur_node ] = right_lcb;
		}
	}
	
	for( uint childI = 0; childI < t[ cur_node ].children.size(); childI++ ){
		// make the invert node follow through the tree down from the initial inversion
		node_id_t new_ins_node = cur_node == invert_node ? t[ cur_node ].children[ childI ] : invert_node;
		recursiveInvert( t[ cur_node ].children[ childI ], new_ins_node, t, lcb_list, first_lcbs, left_lcb, right_lcb );
	}

}

/**
 * filter all gap characters from a sequence
 */
void filterGaps( const string& gapped_seq, string& ungapped_seq ){
	ungapped_seq.clear();
	for( gnSeqI charI = 0; charI < gapped_seq.length(); charI++ )
		if( gapped_seq[ charI ] != '-' )
			ungapped_seq += gapped_seq[ charI ];
}

/**
 * Applying the inversions:
 *
 *	1) Track alignment columns for inversions
 *	2) Perform inversions using LCB structs referencing alignment columns
 *	3) Collapse gaps in the LCB structs
 *  4) Generate correctly rearranged, ungapped sequence
 *     and update the sequence coordinates in LCB structs 
 *
 * Caution:  No user serviceable parts inside
 */
void Alignment::applyInversions( const PhyloTree<TreeNode>& t ) {
	if( alignment.size() != 1 )
		Throw_gnExMsg( EvolutionError(), "Inversions can only be applied once\n" );
	
	uint seq_count = alignment[ 0 ].sequences.size();
	uint seqI;
	uint lcbI;
	std::list< pair< gnSeqI, gnSeqI > >::iterator inv_iter;

	// Step B.2: apply the inversions to LCBs
	
	// initialize the LCB list with a single LCB
	lcb_list = vector< LCB >( 1 );
	first_lcbs = vector< uint >( seq_count, 0 );
	for( seqI = 0; seqI < seq_count; seqI++ ){
		lcb_list[ 0 ].left_end.push_back( 1 );
		lcb_list[ 0 ].right_end.push_back( intervals[ 0 ].length() );
	}
	lcb_list[ 0 ].left_adjacency = vector< uint >( seq_count, (uint)-1 );
	lcb_list[ 0 ].right_adjacency = vector< uint >( seq_count, (uint)-1 );
	lcb_list[ 0 ].lcb_id = 0;

	// start inverting
	for( seqI = 0; seqI < seq_count; seqI++ ){
		inv_iter = inversion_list[ seqI ].begin();
		for( ; inv_iter != inversion_list[ seqI ].end(); inv_iter++ ){
			gnSeqI left_end = inv_iter->first;
			gnSeqI length = inv_iter->second;
			// check wether it's necessary to break LCBs for this inversion,
			// first on the left_end, then on left_end + length (right_end)
			gnSeqI breakpoint = left_end;
			lcbI = first_lcbs[ seqI ];
			uint left_lcb, right_lcb;
			gnSeqI block_sum = 0;	// tally of how far into the sequence we've scanned
			for( uint i = 0; i < 2; i++ ){
				gnSeqI block_len;
				while( (block_len = absolut( lcb_list[ lcbI ].right_end[ seqI ] ) - absolut( lcb_list[ lcbI ].left_end[ seqI ] ) + 1 )
						+ block_sum <= breakpoint ){
					block_sum += block_len;
					lcbI = lcb_list[ lcbI ].right_adjacency[ seqI ];
				}
				// get the distance from the start of the LCB where we will break
				gnSeqI break_start_offset;
				if( lcb_list[ lcbI ].left_end[ seqI ] < 0 )
					break_start_offset = block_len + block_sum - breakpoint;
				else
					break_start_offset = breakpoint - block_sum;

				if( i == 0 )
					left_lcb = lcbI;
				else
					right_lcb = lcbI;
				if( block_sum < breakpoint ){
					// break the LCB up into two pieces
					LCB new_lcb = lcb_list[ lcbI ];
					new_lcb.lcb_id = lcb_list.size();
					// set the start positions and adjacencies
					uint seqJ;
					for( seqJ = 0; seqJ < seq_count; seqJ++ ){
						int64 break_coord = 0;
						// calculate the breakpoint coordinate
						if( new_lcb.left_end[ seqJ ] < 0 ){
							break_coord = new_lcb.right_end[ seqJ ] + break_start_offset;
						}else{
							break_coord = new_lcb.left_end[ seqJ ] + break_start_offset;
						}
						// set the breakpoint coordinate and
						// link up the new LCB into the list
						if( new_lcb.left_end[ seqJ ] < 0 ){
							uint left_adj = lcb_list[ lcbI ].left_adjacency[ seqJ ];
							if( left_adj != (uint)-1 )
								lcb_list[ left_adj ].right_adjacency[ seqJ ] = lcb_list.size();
							else
								first_lcbs[ seqJ ] = lcb_list.size();
							lcb_list[ lcbI ].left_adjacency[ seqJ ] = lcb_list.size();
							new_lcb.right_adjacency[ seqJ ] = lcb_list[ lcbI ].lcb_id;

							new_lcb.right_end[ seqJ ] = break_coord;
							lcb_list[ lcbI ].left_end[ seqJ ] = break_coord;
							if( new_lcb.left_end[ seqI ] < 0 )
								lcb_list[ lcbI ].left_end[ seqJ ]--;
							else
								new_lcb.right_end[ seqJ ]++;
							
						}else{
							uint right_adj = lcb_list[ lcbI ].right_adjacency[ seqJ ];
							if( right_adj != (uint)-1 )
								lcb_list[ right_adj ].left_adjacency[ seqJ ] = lcb_list.size();
							lcb_list[ lcbI ].right_adjacency[ seqJ ] = lcb_list.size();
							new_lcb.left_adjacency[ seqJ ] = lcb_list[ lcbI ].lcb_id;
							new_lcb.left_end[ seqJ ] = break_coord;
							lcb_list[ lcbI ].right_end[ seqJ ] = break_coord;
							if( new_lcb.left_end[ seqI ] < 0 )
								lcb_list[ lcbI ].right_end[ seqJ ]--;
							else
								new_lcb.left_end[ seqJ ]++;
						}
					}
					if( debugChecks() ){
						for( seqJ = 0; seqJ < seq_count; seqJ++ ){
							int64 n_l_j = new_lcb.left_end[ seqJ ];
							int64 n_r_j = new_lcb.right_end[ seqJ ];
							int64 n_l_0 = new_lcb.left_end[ 0 ];
							int64 n_r_0 = new_lcb.right_end[ 0 ];
							int64 n_j_diff = absolut( n_r_j ) - absolut( n_l_j ) + 1;
							int64 n_0_diff = absolut( n_r_0 ) - absolut( n_l_0 ) + 1;
							int64 o_l_j = lcb_list[ lcbI ].left_end[ seqJ ];
							int64 o_r_j = lcb_list[ lcbI ].right_end[ seqJ ];
							int64 o_l_0 = lcb_list[ lcbI ].left_end[ 0 ];
							int64 o_r_0 = lcb_list[ lcbI ].right_end[ 0 ];
							int64 o_j_diff = absolut( o_r_j ) - absolut( o_l_j ) + 1;
							int64 o_0_diff = absolut( o_r_0 ) - absolut( o_l_0 ) + 1;
							if( ( n_j_diff != n_0_diff ) || o_j_diff != o_0_diff ){
								cerr << "Big mistake.  Inversion Perversion!\n";
							}
							if( ( n_j_diff < 0 ) || o_j_diff < 0 ){
								cerr << "WTF?! Negative lcb length!\n";
							}
						}
					}
					if( lcb_list[ lcbI ].left_end[ seqI ] < 0 ){
						block_sum = breakpoint + 1;	// we're skipping over the new LCB so add it's length
						if( i == 1 ){
							right_lcb = lcb_list.size();
							// make sure that left_lcb is still to the left!
							if( left_lcb == lcbI )
								left_lcb = lcb_list.size();
						}
					}
					if( i == 0 && lcb_list[ lcbI ].left_end[ seqI ] > 0 )
						left_lcb = lcb_list.size();
					lcb_list.push_back( new_lcb );
					
					if( debugChecks() ){
						// scan through the LCB list for a sanity check
						for( seqJ = 0; seqJ < seq_count; seqJ++ ){
							vector< bool > visited( lcb_list.size(), false );
							vector< uint > last_lcbs( seq_count, 0 );
							uint vis_lcb = first_lcbs[ seqJ ];
							while( vis_lcb != (uint)-1 ){
								visited[ vis_lcb ] = true;
								if( lcb_list[ vis_lcb ].right_adjacency[ seqJ ] == (uint)-1 )
									last_lcbs[ seqJ ] = vis_lcb;
								else if( lcb_list[ lcb_list[ vis_lcb ].right_adjacency[ seqJ ] ].left_adjacency[ seqJ ] != vis_lcb )
									cerr << "Asymmetric link\n";
								vis_lcb = lcb_list[ vis_lcb ].right_adjacency[ seqJ ];
							}
							for( uint visI = 0; visI < visited.size(); visI++ ){
								if( visited[ visI ] == false )
									cerr << "Right traverse missed lcb: " << visI << endl;
							}

							visited = vector< bool >( lcb_list.size(), false );
							vis_lcb = last_lcbs[ seqJ ];
							while( vis_lcb != (uint)-1 ){
								visited[ vis_lcb ] = true;
								vis_lcb = lcb_list[ vis_lcb ].left_adjacency[ seqJ ];
							}
							for( uint visI = 0; visI < visited.size(); visI++ ){
								if( visited[ visI ] == false )
									cerr << "Left traverse missed lcb: " << visI << endl;
							}
						}
					}
				}
				// now do the same thing with the right end
				breakpoint = left_end + length;
			}
			if( debugChecks() )
				if( lcb_list[ left_lcb ].left_adjacency[ seqI ] == right_lcb )
					cerr << "kaboom\n";
			recursiveInvert( t.root, seqI, t, lcb_list, first_lcbs, left_lcb, right_lcb );
		}
	}

	// step B.3
	// break up alignment into LCBs while filtering gap columns
	
	// scan thru to delete any columns that are all gaps
	vector< string > evolved_seqs = alignment[ 0 ].sequences;
	vector< string > names = alignment[ 0 ].names;
	alignment.clear();
	for( lcbI = 0; lcbI < lcb_list.size(); lcbI++ ){
		gnAlignedSequences gnas;
		alignment.push_back( gnas );
		alignment[ lcbI ].sequences = vector< string >( seq_count );
		alignment[ lcbI ].names = names;
		vector< GappedIntervalSequence< SourceInterval >::iterator > cur_iters( intervals.size() );
		vector< gnSeqI > cur_starts( intervals.size() );
		gnSeqI lcb_len = absolut( lcb_list[ lcbI ].right_end[ 0 ] ) - absolut( lcb_list[ lcbI ].left_end[ 0 ] ) + 1;
		for( gnSeqI baseI = 0; baseI < lcb_len; baseI++ ){
			bool found_char = false;
			for( seqI = 0; seqI < seq_count; seqI++ ){
				cur_iters[ seqI ] = intervals[ seqI ].find( absolut( lcb_list[ lcbI ].left_end[ 0 ] ) + baseI - 1);
				cur_starts[ seqI ] = intervals[ seqI ].getStart( cur_iters[ seqI ] );
				if( !cur_iters[ seqI ]->isGap() )
					found_char = true;
			}
			// add the column
			if( found_char ){
				for( seqI = 0; seqI < seq_count; seqI++ ){
					if( cur_iters[ seqI ]->isGap() ){
						alignment[ lcbI ].sequences[ seqI ] += '-';
						continue;
					}
					gnSeqI source_base = cur_iters[ seqI ]->leftEnd();
					source_base += absolut( lcb_list[ lcbI ].left_end[ 0 ] ) - 1 + baseI - cur_starts[ seqI ];
					gnSeqC new_char;
					if( cur_iters[ seqI ]->sourceID() == TARGET_ID ){
						 new_char = evolved_seqs[ seqI ][ source_base ];
					}else{
						new_char = donor.sequences[ seqI ][ source_base ];
					}
//					if( new_char != 'A' && new_char != 'C' && new_char != 'G' && new_char != 'T' )
//						cerr << new_char;
					alignment[ lcbI ].sequences[ seqI ] += new_char;
				}
			}
		}
	}

	// Step B.4: Generate correctly rearranged, ungapped sequence
	// 			 Simultaneously update the coordinates in the LCB structs to
	//			 reflect the evolved sequence coordinates
	ungapped_seqs = vector< string > ( seq_count );
	const gnFilter *comp_filter = gnFilter::DNAComplementFilter();
	for( seqI = 0; seqI < seq_count; seqI++ ){
		uint cur_lcb = first_lcbs[ seqI ];
		// check for complete coverage
		if( debugChecks() && seqI == 0 ){
			while( cur_lcb != (uint)-1 ){
				uint next_lcb = lcb_list[ cur_lcb ].right_adjacency[ seqI ];
				if( next_lcb == (uint)-1 )
					break;
				int64 cc_left = lcb_list[ cur_lcb ].left_end[ seqI ];
				int64 cc_right = lcb_list[ cur_lcb ].right_end[ seqI ];
				int64 cn_left = lcb_list[ next_lcb ].left_end[ seqI ];
				int64 cn_right = lcb_list[ next_lcb ].right_end[ seqI ];

				if( absolut( lcb_list[ cur_lcb ].right_end[ seqI ] ) + 1 != absolut( lcb_list[ next_lcb ].left_end[ seqI ] ) )
					cerr << "Incorrect linkage or incomplete coverage\n";
				cur_lcb = next_lcb;
			}
			cur_lcb = first_lcbs[ seqI ];
		}
		while( cur_lcb != (uint)-1 ){
			string& gap_seq = alignment[ cur_lcb ].sequences[ seqI ];
			string nogap_seq;
			filterGaps( gap_seq, nogap_seq );
			bool revcomp = lcb_list[ cur_lcb ].left_end[ seqI ] < 0;
			if( revcomp )
				comp_filter->ReverseFilter( nogap_seq );
			lcb_list[ cur_lcb ].left_end[ seqI ] = ungapped_seqs[ seqI ].length() + 1;
			ungapped_seqs[ seqI ] += nogap_seq;
			lcb_list[ cur_lcb ].right_end[ seqI ] = ungapped_seqs[ seqI ].length();
			if( revcomp ){
				lcb_list[ cur_lcb ].left_end[ seqI ] *= -1;
				lcb_list[ cur_lcb ].right_end[ seqI ] *= -1;
			}
			cur_lcb = lcb_list[ cur_lcb ].right_adjacency[ seqI ];
		}
	}
	
}

class NotGap
{
public:
	bool operator()( const char& a )
	{
		return a != '-';
	}
};

void Alignment::writeAlignment( ostream& os, const PhyloTree<TreeNode>& t, boolean write_ancestors ) const {

	uint seq_count = alignment[ 0 ].sequences.size();
	uint seqI;
	uint alignI;

	if( debugCheckingLevel() > 1 ){
		// spit it out in a "Multi-CLUSTALW" format, nice and easy to read:
		for( alignI = 0; alignI < alignment.size(); alignI++ ){
			os << "> ";
			for( seqI = 0; seqI < alignment[ 0 ].sequences.size(); seqI++ ){
				if( seqI > 0 )
					os << '\t';
				os << lcb_list[ alignI ].left_end[ seqI ];
			}
			os << endl;
			alignment[ alignI ].outputClustalW( os );
			os << endl;
		}
	}
	// write it in extended Multi-FastA format (XMFA)
	for( alignI = 0; alignI < alignment.size(); alignI++ ){
		// check for empty LCBs--don't write them
		gnSeqI length_sum = 0;
		NotGap ng;
		for( seqI = 0; seqI < alignment[ alignI ].sequences.size(); seqI++ )
		{
			if( !write_ancestors && t[ seqI ].children.size() > 0 )
				continue;	// only write out the leaf nodes
			length_sum += std::count_if(alignment[ alignI ].sequences[ seqI ].begin(), alignment[ alignI ].sequences[ seqI ].end(), ng);
		}
		if( length_sum == 0 ){
			continue;
		}
		int seq_id = 0;
		for( seqI = 0; seqI < alignment[ 0 ].sequences.size(); seqI++ ){
			if( (!write_ancestors && t[ seqI ].children.size() > 0) || ungapped_seqs[ seqI ].size() == 0)
				continue;	// only write out the leaf nodes which have length > 0
			seq_id++;
			os << "> " << seq_id << ":";
			if( lcb_list[ alignI ].left_end[ seqI ] < 0 ){
				os << -lcb_list[ alignI ].right_end[ seqI ] << '-';
				os << -lcb_list[ alignI ].left_end[ seqI ] << " - ";
			}else{
				os << lcb_list[ alignI ].left_end[ seqI ] << "-";
				os << lcb_list[ alignI ].right_end[ seqI ] << " + ";
			}
			os << alignment[ 0 ].names[ seqI ] << endl;
			gnSeqI cur_pos = 0;
			for( ; cur_pos < alignment[ alignI ].sequences[ seqI ].length(); cur_pos += 80 ){
				gnSeqI cur_len = cur_pos + 80 < alignment[ alignI ].sequences[ seqI ].length() ? 80 : alignment[ alignI ].sequences[ seqI ].length() - cur_pos;
				os.write( alignment[ alignI ].sequences[ seqI ].data() + cur_pos, cur_len );
				os << endl;
			}
		}
		os << "=" << endl << endl;
	}

}

void Alignment::writeEvolvedSequences( ostream& os, const PhyloTree<TreeNode>& t, boolean write_ancestors ) const {
	// write out the ungapped evolved sequences
	uint seq_count = alignment[ 0 ].sequences.size();
	uint seqI;
	gnSequence out_seq;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		if( !write_ancestors && t[ seqI ].children.size() > 0 )
			continue;	// only write out the leaf nodes
		out_seq += ungapped_seqs[ seqI ];
		out_seq.setContigName( out_seq.contigListSize() - 1, alignment[ 0 ].names[ seqI ] );
	}
	gnFASSource::Write( out_seq, os, false, false );
}
