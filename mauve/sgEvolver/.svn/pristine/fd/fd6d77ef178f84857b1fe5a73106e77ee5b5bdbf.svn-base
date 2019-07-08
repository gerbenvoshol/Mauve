/*******************************************************************************
 * $Id: scoreAlignment.cpp,v 1.14 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libGenome/gnFilter.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Matrix.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/Aligner.h"

using namespace std;
using namespace genome;
using namespace mems;

// basic data structures

/** store a pair of aligned positions and the characters */
typedef struct aligned_coords_s {
	int64 pos1;
	int64 pos2;
	char char1;
	char char2;
} aligned_coords_t;


class AlignedCoordSeqIComparator {
public:
	boolean operator()(const aligned_coords_t& a, const aligned_coords_t& b) const{
		if( absolut(a.pos1) == absolut(b.pos1) )
			return absolut(a.pos2) < absolut(b.pos2);
		return absolut(a.pos1) < absolut(b.pos1);
	}
};


void constructCoordList( uint seqI, uint seqJ, IntervalList& iv_list, vector< aligned_coords_t >& coord_list, vector< gnSequence* >& seq_table ){

	//
	// pre-allocate the vector
	//
	gnSeqI ij_vec_size = 0;
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		ij_vec_size += iv_list[ivI].AlignmentLength();
	}
	coord_list = vector< aligned_coords_t >( ij_vec_size );

	//
	// fill in the vector with all aligned pairs
	//
	gnSeqI vecI = 0;	// current place in vector
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		GappedAlignment* aln;
		aln = dynamic_cast< GappedAlignment* >( iv_list[ ivI ].GetMatches()[0] );
		if( aln == NULL ){
			throw "Error:  expecting interval to contain a single GappedAlignment";
		}
		int64 pos1 = aln->Start( seqI );
		int64 pos2 = aln->Start( seqJ );
		// if rev. comp then we're starting at the other (right) side
		if( pos1 < 0 )
			pos1 -= aln->Length( seqI ) - 1;
		if( pos2 < 0 )
			pos2 -= aln->Length( seqJ ) - 1;

		const std::vector< std::string >& align_matrix = GetAlignment( *aln, seq_table );
		for( gnSeqI colI = 0; colI < aln->Length(); colI++ ){
			aligned_coords_t act;
			act.char1 = align_matrix[ seqI ][ colI ];
			act.char2 = align_matrix[ seqJ ][ colI ];
			act.pos1 = act.char1 == '-' ? 0 : pos1;
			act.pos2 = act.char2 == '-' ? 0 : pos2;
			coord_list[ vecI++ ] = act;
			if( act.char1 != '-' )
				pos1++;
			if( act.char2 != '-' )
				pos2++;
		}
	}

	//
	// sort the vector on aligned position
	//
	AlignedCoordSeqIComparator acsc;
	sort( coord_list.begin(), coord_list.end(), acsc );
}


const gnFilter* comp_filter = gnFilter::DNAComplementFilter();

void sanityCheck( string& aln_name, uint seqI, uint seqJ, aligned_coords_t& act, vector< string >& evolved_seqs ){
	if( act.pos1 != 0 ){
		char seqI_char = evolved_seqs[ seqI ][ absolut( act.pos1 ) - 1 ];
		if( act.pos1 < 0 )
			seqI_char = comp_filter->Filter( seqI_char );
		if( act.char1 != seqI_char ){
			cerr << "Error: " << aln_name << " character " << act.char1 << " instead of " << seqI_char << " at position " << act.pos1 << endl;
		}
	}

	if( act.pos2 != 0 ){
		char seqJ_char = evolved_seqs[ seqJ ][ absolut( act.pos2 ) - 1 ];
		if( act.pos2 < 0 )
			seqJ_char = comp_filter->Filter( seqJ_char );

		if( act.char2 != seqJ_char ){
			cerr << "Error: " << aln_name << " character " << act.char2 << " instead of " << seqJ_char << " at position " << act.pos2 << endl;
		}
	}
}

bool warn_missing = true;
int warn_count = 0;
const int WARN_MAX = 1000;

void compareAlignments( IntervalList& correct, IntervalList& calculated, vector< string >& evolved_seqs, vector< gnSequence* >& seq_table ){
	
	string cor_name = "correct";
	string calc_name = "calculated";
	gnSeqI sp_truepos = 0;
	gnSeqI sp_falsepos = 0;
	gnSeqI sp_trueneg = 0;
	gnSeqI sp_falseneg = 0;
	gnSeqI sp_extraneg = 0;
	gnSeqI wrong_strand = 0;

	// calculate total possible score
	gnSeqI sp_possible = 0;

	uint seqI = 0;
	uint seqJ = 0;

	// now for every entry in the correct alignment, look for a corresponding
	// entry in the calculated alignment.  Do sanity checking while we're at it
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			//
			// the correct coordinate matrix is expected to have one entry for
			// each character in seqI, plus an arbitrary number of "zero" entries
			// for other sequences (seqJ > seqI) that couldn't align to anything
			// in seqI
			//
			vector< aligned_coords_t > cor;
			constructCoordList( seqI, seqJ, correct, cor, seq_table );

			//
			// the calculated alignment is expected to account for every residue
			// at least once, either aligning it to another residue or a gap
			//
			vector< aligned_coords_t > calc;
			constructCoordList( seqI, seqJ, calculated, calc, seq_table );

			gnSeqI corI = 0;
			gnSeqI calcI = 0;

			// skip any gaps aligned to gaps
			while( corI < cor.size() && cor[ corI ].pos1 == 0 && cor[ corI ].pos2 == 0 )
				corI++;
			while( calcI < calc.size() && calc[ calcI ].pos1 == 0 && calc[ calcI ].pos2 == 0 )
				calcI++;

			//
			// get thru the unaligned region of seqI (where pos1 is zero)
			//
			while( corI < cor.size() && calcI < calc.size() && 
				cor[ corI ].pos1 == 0 && calc[ calcI ].pos1 == 0 )
			{
				// sanity check the "correct" aligned characters...
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );

				while( calcI < calc.size() && calc[ calcI ].pos1 == 0 &&
					absolut(calc[ calcI ].pos2) < absolut(cor[ corI ].pos2) )
				{
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					// these are all false positives, catch them later
					calcI++;
				}

				uint match_count = 0;
				while( calcI < calc.size() && calc[ calcI ].pos1 == 0 &&
					absolut(calc[ calcI ].pos2) == absolut(cor[ corI ].pos2) ){
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					match_count++;
					calcI++;
				}
				sp_trueneg += match_count > 0 ? 1 : 0;
				sp_extraneg += match_count > 0 ? match_count - 1 : 0;
				// catch false positives later
				corI++;
			}

			// scan thru any remaining false negatives
			// that were in the calculated alignment
			while( calcI < calc.size() && calc[ calcI ].pos1 == 0 ){
				// sanity check the aligned characters...
				sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
				calcI++;
				sp_falseneg++;
			}
			while( corI < cor.size() && cor[ corI ].pos1 == 0 ){
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );
				// this is a false positive...
				corI++;
			}
			
			//
			// now look at the rest of seqI
			//
			for( ; corI < cor.size(); corI++ ){
				// sanity check the "correct" aligned characters
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );

				// if this is an aligned pair it gets counted towards
				// the total possible points in the sum of pairs score
				if( cor[ corI ].pos2 != 0 )
					sp_possible++;

				// make sure calc is up to where cor is...
				while( calcI < calc.size() && (absolut(calc[ calcI ].pos1) < absolut(cor[ corI ].pos1)) ){
					// sanity check the aligned characters
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );

					// cor skipped more than one char and is missing something...
					cerr << "correct alignment is missing characters\n";
					calcI++;
				}

				// make sure cor is up to where calc is...
				if( calcI == calc.size() || (absolut(calc[ calcI ].pos1) > absolut(cor[ corI ].pos1)) ){
					// calc is missing something!
					// assume it's a gap
					if( warn_missing )
					{
						cerr << "Warning, calculated alignment missing prediction for seqI " << seqI << " position " << absolut(cor[ corI ].pos1) << ".  Assuming gap prediction\n";
						warn_count++;
						if( warn_count > WARN_MAX )
						{
							cerr << "Too many warnings... no longer reporting\n";
							warn_missing = false;
						}
					}
					if( cor[ corI ].pos2 == 0 )
						sp_trueneg++;
					else
						sp_falseneg++;
					continue;
				}

				// if all is well, we should have found the corresponding entries in
				// cor and calc
				int cur_count = 0;
				// examing all predictions for this residue in the calculated alignment...
				while( calcI < calc.size() && (absolut(calc[ calcI ].pos1) == absolut(cor[ corI ].pos1)) ){
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					
					// check the strands that were aligned
					bool cor_parity = ((cor[ corI ].pos1 > 0 && cor[ corI ].pos2 > 0) ||
										(cor[ corI ].pos1 < 0 && cor[ corI ].pos2 < 0)) ?
										true : false;

					bool calc_parity = ((calc[ calcI ].pos1 > 0 && calc[ calcI ].pos2 > 0) ||
										(calc[ calcI ].pos1 < 0 && calc[ calcI ].pos2 < 0)) ?
										true : false;

					// if they match then increment the sp_score
					if( absolut(calc[ calcI ].pos2) == absolut(cor[ corI ].pos2) ){
						if( calc[ calcI ].pos2 == 0 )
							sp_trueneg++;	// correct alignment to a gap
						else if( cor_parity == calc_parity )
							sp_truepos++;	// correct positions on correct strand!!
						else
						{
							sp_falsepos++;	// correct positions, wrong strand
							wrong_strand++;
						}
					}else{
						if( calc[ calcI ].pos2 == 0 )
							sp_falseneg++;	// incorrect alignment to a gap
						else
							sp_falsepos++;	// incorrectly aligned

					}
					calcI++;
				}
			}
		}
	}
	
	
	double sp_accuracy = (double)sp_truepos / (double)sp_possible;
	cout << "Sum of pairs accuracy: " << sp_accuracy << endl;
	cout << "Sum of pairs error rate: " << (double)sp_falsepos / (double)sp_truepos << endl;
	cout << "The error rate gives the number of incorrect orthology predictions per correct prediction." << endl;
	cout << "Sum of pairs positive predictive value: " << (double)sp_truepos / (double)(sp_truepos + sp_falsepos) << endl;
	cout << "trueneg: " << sp_trueneg << "  falseneg: " << sp_falseneg << "  extraneg: " << sp_extraneg << endl;
	cout << "wrong_strand: " << wrong_strand << endl;
}








template< typename FwdIt >
bool findFirstCorrectlyAlignedPair( 
	const FwdIt& cor_first,
	const FwdIt& cor_last,
	const FwdIt& calc_first,
	const FwdIt& calc_last,
	FwdIt& cor_iter,
	FwdIt& calc_iter,
	int dir
	)
{
	while( dir == -1 || dir == 1 && calc_iter != calc_last && cor_iter != cor_last )
	{
		if( genome::absolut(calc_iter->pos1) < genome::absolut(cor_iter->pos1) )
		{
			if( dir == 1 )
				++calc_iter;
			else if( cor_iter == cor_last )
				return false;
			else
				--cor_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos1) > genome::absolut(cor_iter->pos1) )
		{
			if( dir == 1 )
				++cor_iter;
			else if( calc_iter == calc_last )
				return false;
			else
				--calc_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos2) != genome::absolut(cor_iter->pos2) ||
			calc_iter->pos2 == 0)
		{
			if( calc_iter == calc_last || cor_iter == cor_last )
				return false;
			calc_iter += dir;
			cor_iter += dir;
			continue;
		}
		bool cor_a = cor_iter->pos1 < 0;
		bool cor_b = cor_iter->pos2 < 0;
		bool cor_parity = cor_a == cor_b;
		bool calc_a = calc_iter->pos1 < 0;
		bool calc_b = calc_iter->pos2 < 0;
		bool calc_parity = calc_a == calc_b;
		if( cor_parity != calc_parity )
		{
			// calc has the opposite strand aligned
			++calc_iter;
			++cor_iter;
			continue;
		}
		break;	// these are aligned correctly
	}
	if( dir == 1 && (cor_iter == cor_last || calc_iter == calc_last) )
	{
		calc_iter = calc_last;
		return false;
	}
}

size_t findLcb( vector< LCB >& calc_pair_adj, int64 pos )
{
	size_t calc_ivI = 0;
	for( ; calc_ivI < calc_pair_adj.size(); ++calc_ivI )
	{
		if( calc_pair_adj[calc_ivI].left_end[0] == NO_MATCH )
			continue;
		if( calc_pair_adj[calc_ivI].left_end[0] <= genome::absolut( pos ) &&
			genome::absolut( pos ) < calc_pair_adj[calc_ivI].right_end[0] )
			break;
	}
	return calc_ivI;
}


void markAllCorrectlyAlignedLcbs( 
	const vector< aligned_coords_t >::iterator& cor_first,
	const vector< aligned_coords_t >::iterator& cor_last,
	const vector< aligned_coords_t >::iterator& calc_first,
	const vector< aligned_coords_t >::iterator& calc_last,
	vector< LCB >& calc_pair_adj,
	boost::dynamic_bitset<>& found_lcbs
	)
{

	vector< aligned_coords_t >::iterator cor_iter = cor_first;
	vector< aligned_coords_t >::iterator calc_iter = calc_first;
	while( calc_iter != calc_last && cor_iter != cor_last )
	{
		if( genome::absolut(calc_iter->pos1) < genome::absolut(cor_iter->pos1) )
		{
			++calc_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos1) > genome::absolut(cor_iter->pos1) )
		{
			++cor_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos2) != genome::absolut(cor_iter->pos2)  ||
			calc_iter->pos2 == 0)
		{
			++calc_iter;
			++cor_iter;
			continue;
		}
		bool cor_a = cor_iter->pos1 < 0;
		bool cor_b = cor_iter->pos2 < 0;
		bool cor_parity = cor_a == cor_b;
		bool calc_a = calc_iter->pos1 < 0;
		bool calc_b = calc_iter->pos2 < 0;
		bool calc_parity = calc_a == calc_b;
		if( cor_parity != calc_parity )
		{
			// calc has the opposite strand aligned
			++calc_iter;
			++cor_iter;
			continue;
		}
		// these are aligned correctly

		size_t calc_ivI = findLcb( calc_pair_adj, calc_iter->pos1 );
		found_lcbs[calc_ivI] = true;
		// move calc_iter to the next LCB
		aligned_coords_t next_act;
		next_act.pos1 = calc_pair_adj[calc_ivI].right_end[0] + 1;
		AlignedCoordSeqIComparator acsi;
		calc_iter = std::lower_bound( calc_first, calc_last, next_act, acsi );
	}
}


void computeLCBaccuracy( IntervalList& correct, IntervalList& calculated, vector< string >& evolved_seqs, vector< gnSequence* >& seq_table )
{	
	string cor_name = "correct";
	string calc_name = "calculated";
	double bp_sp_truepos = 0;
	double bp_sp_falsepos = 0;
	double bp_sp_falseneg = 0;
	uint seqI = 0;
	uint seqJ = 0;


	// find all pairwise LCBs
	Matrix< vector< LCB > > cor_pairwise_adjs(seq_table.size(),seq_table.size());
	Matrix< vector< LCB > > calc_pairwise_adjs(seq_table.size(),seq_table.size());
	for( seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ )
		{
			vector< size_t > projection( 2 );
			projection[0] = seqI;
			projection[1] = seqJ;

			// construct pairwise Interval projections
			vector< AbstractMatch* > cor_mpa_list;
			vector< AbstractMatch* > calc_mpa_list;
			for( size_t corI = 0; corI < correct.size(); corI++ )
			{
				if( correct[corI].LeftEnd(seqI) == NO_MATCH || correct[corI].LeftEnd(seqJ) == NO_MATCH )
					continue;
				MatchProjectionAdapter mpa_tmp( &correct[corI], projection );
				vector< bitset_t > mpa_aln;
				mpa_tmp.GetAlignment(mpa_aln);
				bitset_t intersect = mpa_aln[0] & mpa_aln[1];
				if(!intersect.any())
					continue;	// nothing was aligned among seqI and seqJ in this interval, so skip it!
				cor_mpa_list.push_back( mpa_tmp.Copy() );
				if( cor_mpa_list.back()->Orientation(0) == AbstractMatch::reverse )
					cor_mpa_list.back()->Invert();
			}
			for( size_t calcI = 0; calcI < calculated.size(); calcI++ )
			{
				if( calculated[calcI].LeftEnd(seqI) == NO_MATCH || calculated[calcI].LeftEnd(seqJ) == NO_MATCH )
					continue;
				MatchProjectionAdapter mpa_tmp( &calculated[calcI], projection );
				vector< bitset_t > mpa_aln;
				mpa_tmp.GetAlignment(mpa_aln);
				bitset_t intersect = mpa_aln[0] & mpa_aln[1];
				if(!intersect.any())
					continue;	// nothing was aligned among seqI and seqJ in this interval, so skip it!
				calc_mpa_list.push_back( mpa_tmp.Copy() );
				if( calc_mpa_list.back()->Orientation(0) == AbstractMatch::reverse )
					calc_mpa_list.back()->Invert();
			}
			vector< vector< AbstractMatch* > > LCB_list;
			vector< gnSeqI > breakpoints;
			IdentifyBreakpoints( cor_mpa_list, breakpoints );
			ComputeLCBs_v2( cor_mpa_list, breakpoints, LCB_list );
			vector< double > lcb_scores( LCB_list.size(), 0 );
			computeLCBAdjacencies_v3( LCB_list, lcb_scores, cor_pairwise_adjs(seqI,seqJ) );

			breakpoints.clear();
			LCB_list.clear();
			IdentifyBreakpoints( calc_mpa_list, breakpoints );
			ComputeLCBs_v2( calc_mpa_list, breakpoints, LCB_list );
			lcb_scores = vector< double >( LCB_list.size(), 0 );
			computeLCBAdjacencies_v3( LCB_list, lcb_scores, calc_pairwise_adjs(seqI,seqJ) );
		}
	}

	// calculate total possible score
	double bp_sp_possible = 0;
	double bp_dist_sum = 0;
	double bp_dist_count = 0;
	vector< double > bp_dist;

	Matrix< vector<double> > left_dists( seq_table.size(), seq_table.size() );
	Matrix< vector<double> > right_dists( seq_table.size(), seq_table.size() );

	Matrix< boost::dynamic_bitset<> > found_lcbs( seq_table.size(), seq_table.size() );
	for( seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ )
		{
			boost::dynamic_bitset<> tmp_bs ( calc_pairwise_adjs(seqI,seqJ).size(), false );	// tracks whether calculated LCBs were found in the correct alignment
			found_lcbs(seqI,seqJ) = tmp_bs;
			vector<double> asdf( calc_pairwise_adjs(seqI,seqJ).size(), (std::numeric_limits<double>::max)() );
			left_dists(seqI, seqJ) = asdf;
			left_dists(seqJ, seqI) = asdf;
			right_dists(seqI, seqJ) = asdf;
			right_dists(seqJ, seqI) = asdf;
		}
	}

	// now for every entry in the correct alignment, look for a corresponding
	// entry in the calculated alignment.  Do sanity checking while we're at it
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			//
			// the correct coordinate matrix is expected to have one entry for
			// each character in seqI, plus an arbitrary number of "zero" entries
			// for other sequences (seqJ > seqI) that couldn't align to anything
			// in seqI
			//
			vector< aligned_coords_t > cor;
			constructCoordList( seqI, seqJ, correct, cor, seq_table );

			//
			// the calculated alignment is expected to account for every residue
			// at least once, either aligning it to another residue or a gap
			//
			vector< aligned_coords_t > calc;
			constructCoordList( seqI, seqJ, calculated, calc, seq_table );

			gnSeqI corI = 0;
			gnSeqI calcI = 0;

			vector< LCB >& cor_pair_adj = cor_pairwise_adjs(seqI,seqJ);
			vector< LCB >& calc_pair_adj = calc_pairwise_adjs(seqI,seqJ);
			for( size_t cor_ivI = 0; cor_ivI < cor_pair_adj.size(); cor_ivI++ )
			{
				bp_sp_possible++;

				int64 lend_seqI = cor_pair_adj[cor_ivI].left_end[0];
				int64 rend_seqI = cor_pair_adj[cor_ivI].right_end[0]-1;

				aligned_coords_t act_left;
				aligned_coords_t act_right;
				act_left.pos1 = lend_seqI;
				act_left.pos2 = (std::numeric_limits<int64>::min)();
				act_right.pos1 = rend_seqI;
				act_right.pos2 = (std::numeric_limits<int64>::max)();
				AlignedCoordSeqIComparator acsi;
				vector< aligned_coords_t >::iterator cor_first;
				vector< aligned_coords_t >::iterator cor_last;
				cor_first = std::lower_bound( cor.begin(), cor.end(), act_left, acsi );
				cor_last = std::upper_bound( cor.begin(), cor.end(), act_right, acsi );

				vector< aligned_coords_t >::iterator calc_first;
				vector< aligned_coords_t >::iterator calc_last;
				calc_first = std::lower_bound( calc.begin(), calc.end(), act_left, acsi );
				calc_last = std::upper_bound( calc.begin(), calc.end(), act_right, acsi );

				// find the left-most LCB that has a correctly aligned character
				vector< aligned_coords_t >::iterator calc_iter = calc_first;
				vector< aligned_coords_t >::iterator cor_iter = cor_first;
				findFirstCorrectlyAlignedPair( cor_first, cor_last, calc_first, calc_last, cor_iter, calc_iter, 1 );
				if( calc_iter == calc_last )
				{
					bp_sp_falseneg++;	// didn't align anything correctly in this region
					continue;
				}
				bp_sp_truepos++;	// found this LCB!

				markAllCorrectlyAlignedLcbs( cor_first, cor_last, calc_first, calc_last, calc_pair_adj, found_lcbs(seqI,seqJ) );

				// determine which LCB in the calculated ivs we're inside of
				size_t calc_ivI = findLcb( calc_pair_adj, calc_iter->pos1 );
				bp_dist.push_back( (double)lend_seqI - (double)calc_pair_adj[calc_ivI].left_end[0] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;
				if( genome::absolut( bp_dist.back() ) < genome::absolut( left_dists(seqI,seqJ)[calc_ivI] ) )
					left_dists(seqI,seqJ)[calc_ivI] = bp_dist.back();

				bool parity = (cor_pair_adj[cor_ivI].left_end[1] < 0) == (calc_pair_adj[calc_ivI].left_end[1] < 0);
				if( parity )
					bp_dist.push_back( (double)cor_pair_adj[cor_ivI].left_end[1] - (double)calc_pair_adj[calc_ivI].left_end[1] );
				else
				{
					bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[1] + (double)cor_pair_adj[cor_ivI].right_end[1] );
					if( calc_pair_adj[calc_ivI].left_end[1] < 0 )
						bp_dist.back() *= -1;
				}
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;

				if( genome::absolut( bp_dist.back() ) < genome::absolut( left_dists(seqJ,seqI)[calc_ivI] ) )
					left_dists(seqJ,seqI)[calc_ivI] = bp_dist.back();

				// now do the same for the other side
				// determine how far away the end of the LCB prediction is
				vector< aligned_coords_t >::iterator rcor_first;
				vector< aligned_coords_t >::iterator rcor_last;
				rcor_first = cor_last - 1;
				rcor_last = cor_first;

				vector< aligned_coords_t >::iterator rcalc_first;
				vector< aligned_coords_t >::iterator rcalc_last;
				rcalc_first = calc_last - 1;
				rcalc_last = calc_first;

				// find the left-most LCB that has a correctly aligned character
				vector< aligned_coords_t >::iterator rcalc_iter = rcalc_first;
				vector< aligned_coords_t >::iterator rcor_iter = rcor_first;
				bool found = findFirstCorrectlyAlignedPair( rcor_first, rcor_last, rcalc_first, rcalc_last, rcor_iter, rcalc_iter, -1 );

				// determine which LCB in the calculated ivs we're inside of
				calc_ivI = findLcb( calc_pair_adj, rcalc_iter->pos1 );
				bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[0] - (double)cor_pair_adj[cor_ivI].right_end[0] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;
				if( genome::absolut( bp_dist.back() ) < genome::absolut( right_dists(seqI,seqJ)[calc_ivI] ) )
					right_dists(seqI,seqJ)[calc_ivI] = bp_dist.back();

				parity = (cor_pair_adj[cor_ivI].left_end[1] < 0) == (calc_pair_adj[calc_ivI].left_end[1] < 0);
				if( !parity )
				{
					bp_dist.push_back( (double)cor_pair_adj[cor_ivI].left_end[1] + (double)calc_pair_adj[calc_ivI].left_end[1] );
					if( calc_pair_adj[calc_ivI].left_end[1] > 0 )
						bp_dist.back() *= -1;
				}else
					bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[1] - (double)cor_pair_adj[cor_ivI].right_end[1] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;

				if( genome::absolut( bp_dist.back() ) < genome::absolut( right_dists(seqJ,seqI)[calc_ivI] ) )
					right_dists(seqJ,seqI)[calc_ivI] = bp_dist.back();
			}
		}
	}

	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			// compute falsepos as the number of lcbs that weren't found in the correct aln
			bp_sp_falsepos += found_lcbs(seqI, seqJ).size() - found_lcbs(seqI, seqJ).count();
		}
	}

	double bp_sp_sensitivity = bp_sp_truepos / bp_sp_possible;
	double bp_sp_ppv = bp_sp_truepos / (bp_sp_truepos + bp_sp_falsepos);
	double bp_sp_mean_dist = bp_dist_sum / bp_dist_count;
	double bp_sp_stddev_dist = 0;
	for( size_t i = 0; i < bp_dist.size(); ++i )
		bp_sp_stddev_dist += (bp_dist[i] - bp_sp_mean_dist) * (bp_dist[i] - bp_sp_mean_dist);
	bp_sp_stddev_dist /= bp_dist.size()-1;
	bp_sp_stddev_dist = sqrt( bp_sp_stddev_dist );

	double bp_sp_mean_2 = 0;
	double bp_sp_count_2 = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
				{
					bp_sp_mean_2 += left_dists(seqI,seqJ)[lcbI];
					bp_sp_count_2++;
				}
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
				{
					bp_sp_mean_2 += right_dists(seqI,seqJ)[lcbI];
					bp_sp_count_2++;
				}
			}
		}
	}
	bp_sp_mean_2 /= bp_sp_count_2;

	double bp_sp_stddev_2 = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					bp_sp_stddev_2 += (left_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2) * (left_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2);
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					bp_sp_stddev_2 += (right_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2) * (right_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2);
			}
		}
	}
	bp_sp_stddev_2 /= (bp_sp_count_2 - 1);
	bp_sp_stddev_2 = sqrt(bp_sp_stddev_2);

	cout << "Sum of pairs LCB sensitivity: " << bp_sp_sensitivity << std::endl;
	cout << "Sum of pairs LCB positive predictive value: " << bp_sp_ppv << std::endl;
//	cout << "Mean distance between predicted breakpoint and true breakpoint: " << bp_sp_mean_dist << std::endl;
//	cout << "Standard deviation: " << bp_sp_stddev_dist << std::endl;

//	cout << "\nSecond bp accuracy method:\n";
	cout << "Mean distance between predicted breakpoint and true breakpoint: " << bp_sp_mean_2 << std::endl;
	cout << "Standard deviation: " << bp_sp_stddev_2 << std::endl;

	// calculate quartile statistics
	vector< gnSeqI > all_dists;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					all_dists.push_back( (gnSeqI)genome::absolut( left_dists(seqI,seqJ)[lcbI] ));
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					all_dists.push_back( (gnSeqI)genome::absolut( right_dists(seqI,seqJ)[lcbI] ) );
			}
		}
	}
	std::sort( all_dists.begin(), all_dists.end() );
	if( all_dists.size() > 0 )
	{
		cout << "Absolute distance, min: " << all_dists.front() << endl;
		cout << "Absolute distance, first quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.25) ] << endl;
		cout << "Absolute distance, second quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.5) ] << endl;
		cout << "Absolute distance, third quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.75) ] << endl;
		cout << "Absolute distance, max: " << all_dists.back() << endl;
	}
}

/** records the aligned segments flanking an indel */
typedef struct indel_s
{
	aligned_coords_t left_block_left;
	aligned_coords_t left_block_right;
	aligned_coords_t right_block_left;
	aligned_coords_t right_block_right;
} indel_t;

typedef struct indel_map_s
{
	vector< size_t > s1_l;	/**< left block mapping for sequence 1 */
	vector< size_t > s1_r;	/**< left block mapping for sequence 1 */
	vector< size_t > s2_l;	/**< right block mapping for sequence 2 */
	vector< size_t > s2_r;	/**< right block mapping for sequence 2 */
} indel_map_t;

const size_t NO_INDEL = (std::numeric_limits<size_t>::max)();

void printIndel( indel_t i )
{
	cerr << "i.left_block_left.pos1 " << i.left_block_left.pos1 << endl;
	cerr << "i.left_block_right.pos1 " << i.left_block_right.pos1 << endl;
	cerr << "--\n";
	cerr << "i.right_block_left.pos1 " << i.right_block_left.pos1 << endl;
	cerr << "i.right_block_right.pos1 " << i.right_block_right.pos1 << endl;
	cerr << "**\n";
	cerr << "i.left_block_left.pos2 " << i.left_block_left.pos2 << endl;
	cerr << "i.left_block_right.pos2 " << i.left_block_right.pos2 << endl;
	cerr << "i.right_block_left.pos2 " << i.right_block_left.pos2 << endl;
	cerr << "i.right_block_right.pos2 " << i.right_block_right.pos2 << endl;
}

bool checkCollinear( const aligned_coords_t& a, const aligned_coords_t& b )
{
	if( (a.pos2 > 0 && b.pos2 < 0) || (a.pos2 < 0 && b.pos2 > 0) )
		return false;
	if( (a.pos1 >= b.pos1 && a.pos2 >= b.pos2 ) ||
		(a.pos1 <= b.pos1 && a.pos2 <= b.pos2 ) )
		return true;
	return false;
}

void constructIndelList( const vector< aligned_coords_t >& coords, vector< indel_t >& indels, indel_map_t& imap,
						size_t seq1_size, size_t seq2_size)
{
	indels.clear();

	// initialize sequence position to indel ID maps
	imap.s1_l.clear();
	imap.s1_l.resize( seq1_size+1, NO_INDEL );
	imap.s1_r.clear();
	imap.s1_r.resize( seq1_size+1, NO_INDEL );
	imap.s2_l.clear();
	imap.s2_l.resize( seq2_size+1, NO_INDEL );
	imap.s2_r.clear();
	imap.s2_r.resize( seq2_size+1, NO_INDEL );

	// skip any leading gaps
	size_t cI = 0;
	while( cI < coords.size() && (coords[ cI ].pos1 == 0 || coords[ cI ].pos2 == 0) )
		cI++;
	if(cI >= coords.size())
		return;	// nothing aligned, so no indels in this pair

	indel_t cur_ind;
	cur_ind.right_block_left = coords[cI];
	size_t prevI = cI;
	bool first = true;
	for( cI++; cI < coords.size(); cI++ )
	{
		// skip gaps
		while( cI < coords.size() && (coords[ cI ].pos1 == 0 || coords[ cI ].pos2 == 0) )
			cI++;
		if(cI >= coords.size() || 
			(coords[cI].pos1 - coords[prevI].pos1 == 1 && coords[cI].pos2 - coords[prevI].pos2 == 1))
		{
			if( cI < coords.size() )
				prevI = cI;
			continue;	// still in an aligned block
		}

		cur_ind.right_block_right = coords[prevI];
		if(!first)
		{
			if(!checkCollinear(cur_ind.left_block_left, cur_ind.left_block_right) ||
				!checkCollinear(cur_ind.left_block_right, cur_ind.right_block_left) ||
				!checkCollinear(cur_ind.right_block_left, cur_ind.right_block_right) )
			{
				printIndel(cur_ind);
				cerr << "ohshit\n";
			}
			indels.push_back( cur_ind );
		}
		first = false;
		if(!checkCollinear(coords[cI], coords[prevI]))
		{
			// hit an LCB boundary
			prevI = cI;
			first = true;
			cur_ind.right_block_left = coords[cI];	// update this so left_block_left gets set properly later
			continue;
		}
		// hit a gap
		cur_ind.left_block_left = cur_ind.right_block_left;
		cur_ind.left_block_right = coords[prevI];
		cur_ind.right_block_left = coords[cI];
		prevI = cI;
	}

	// get the last one
	cur_ind.right_block_right = coords[prevI];
	if(!first)
		indels.push_back( cur_ind );

	// create sequence position to indel ID maps
	for( size_t i = 0; i < indels.size(); i++ )
	{
		size_t s = min( absolut(indels[i].left_block_left.pos1), absolut(indels[i].left_block_right.pos1) );
		size_t e = max( absolut(indels[i].left_block_left.pos1), absolut(indels[i].left_block_right.pos1) );
		for( size_t j = s; j <= e; j++ )
			imap.s1_l[j] = i;
		s = min( absolut(indels[i].right_block_left.pos1), absolut(indels[i].right_block_right.pos1) );
		e = max( absolut(indels[i].right_block_left.pos1), absolut(indels[i].right_block_right.pos1) );
		for( size_t j = s; j <= e; j++ )
			imap.s1_r[j] = i;

		s = min( absolut(indels[i].left_block_left.pos2), absolut(indels[i].left_block_right.pos2) );
		e = max( absolut(indels[i].left_block_left.pos2), absolut(indels[i].left_block_right.pos2) );
		for( size_t j = s; j <= e; j++ )
			imap.s2_l[j] = i;
		s = min( absolut(indels[i].right_block_left.pos2), absolut(indels[i].right_block_right.pos2) );
		e = max( absolut(indels[i].right_block_left.pos2), absolut(indels[i].right_block_right.pos2) );
		for( size_t j = s; j <= e; j++ )
			imap.s2_r[j] = i;
	}
}

// sets seq1 to be the reference genome -- always positive values
void homogenizeParity( vector< aligned_coords_t >& coords )
{
	for( size_t i = 0; i < coords.size(); i++ )
		if( coords[i].pos1 < 0 )
		{
			coords[i].pos1 *= -1;
			coords[i].pos2 *= -1;
		}
		AlignedCoordSeqIComparator acsc;
		sort( coords.begin(), coords.end(), acsc );
}

void computeIndelAccuracy( IntervalList& correct, IntervalList& calculated, vector< string >& evolved_seqs, vector< gnSequence* >& seq_table )
{	
	// how to do this?
	// reduce a pairwise alignment into ungapped blocks
	double indel_sp_truepos = 0;
	double indel_sp_perfectpos = 0;
	double indel_sp_falsepos = 0;
	double indel_sp_falseneg = 0;
	vector< pair< int, int > > indel_bounds;	// left and right boundary scores for all indels that were classified TP
	vector< pair< int, int > > gap_sizes;		// seq_i and seq_j gap sizes for all indels that were classified TP
	bitset_t perfect_prediction;	// true if there were no intervening indels predicted	
	uint seqI = 0;
	uint seqJ = 0;


	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			//
			// the correct coordinate matrix is expected to have one entry for
			// each character in seqI, plus an arbitrary number of "zero" entries
			// for other sequences (seqJ > seqI) that couldn't align to anything
			// in seqI
			//
			vector< aligned_coords_t > cor;
			constructCoordList( seqI, seqJ, correct, cor, seq_table );
			homogenizeParity( cor );

			//
			// the calculated alignment is expected to account for every residue
			// at least once, either aligning it to another residue or a gap
			//
			vector< aligned_coords_t > calc;
			constructCoordList( seqI, seqJ, calculated, calc, seq_table );
			homogenizeParity( calc );

			vector< indel_t > cor_indels;
			indel_map_t cor_map;
			constructIndelList( cor, cor_indels, cor_map, evolved_seqs[seqI].size(), evolved_seqs[seqJ].size() );

			vector< indel_t > calc_indels;
			indel_map_t calc_map;
			constructIndelList( calc, calc_indels, calc_map, evolved_seqs[seqI].size(), evolved_seqs[seqJ].size() );

			// for each indel in the true alignment, check whether it satisfies the
			// criteria for a TP in the predicted alignment:
			// it must have at least one correctly aligned nucleotide pair in the left
			// and right flanking blocks, and one nucleotide correctly aligned to a gap
			// within the indel.
			size_t cur_indel_tp = 0;
			for( size_t corI = 0; corI < cor_indels.size(); corI++ )
			{
				// find the nearest corresponding indel in calc_indels
				size_t calcI_left = 0;
				size_t calcI_right = 0;
				int64 lI;
				for( lI = cor_indels[corI].left_block_right.pos1; lI >= cor_indels[corI].left_block_left.pos1; lI-- )
				{
					calcI_left = calc_map.s1_l[absolut(lI)];
					if( calcI_left == NO_INDEL)
						continue;	// nothing to see here
					else{
						// check whether the aligned positions have the same generalized offset
						// i.e. if they're on the same diagonal.  if so, then we've got a hit
						int64 cor_diag = cor_indels[corI].left_block_left.pos1 - cor_indels[corI].left_block_left.pos2;
						int64 calc_diag = calc_indels[calcI_left].left_block_left.pos1 - calc_indels[calcI_left].left_block_left.pos2;
						if( cor_diag == calc_diag )
							break;
					}
				}
				// did we find a match to the left diagonal?
				if( lI < cor_indels[corI].left_block_left.pos1 )
					continue;

				// find the nearest indel in calc_indels corresponding to the right-side block
				int64 rI;
				for( rI = cor_indels[corI].right_block_left.pos1; rI <= cor_indels[corI].right_block_right.pos1; rI++ )
				{
					calcI_right = calc_map.s1_r[absolut(rI)];
					if( calcI_right == NO_INDEL)
						continue;	// nothing to see here
					else{
						// check whether the aligned positions have the same generalized offset
						// i.e. if they're on the same diagonal.  if so, then we've got a hit
						int64 cor_diag = cor_indels[corI].right_block_left.pos1 - cor_indels[corI].right_block_left.pos2;
						int64 calc_diag = calc_indels[calcI_right].right_block_left.pos1 - calc_indels[calcI_right].right_block_left.pos2;
						if( cor_diag == calc_diag )
							break;
					}
				}
				// did we find a match to the right diagonal?
				if( rI > cor_indels[corI].right_block_right.pos1 )
					continue;

				// did we correctly align at least one site to a gap?
				bool good = 
					(calc_indels[calcI_left].left_block_right.pos1 + 1 < cor_indels[corI].right_block_left.pos1 &&
					calc_indels[calcI_right].right_block_left.pos1 - 1 > cor_indels[corI].left_block_right.pos1 ) ||
					(calc_indels[calcI_left].left_block_right.pos2 + 1 < cor_indels[corI].right_block_left.pos2 &&
					calc_indels[calcI_right].right_block_left.pos2 - 1 > cor_indels[corI].left_block_right.pos2);
				if(good)
				{
					cur_indel_tp++;
					int b1l = cor_indels[corI].left_block_right.pos1 - calc_indels[calcI_left].left_block_right.pos1;
					int b1r = calc_indels[calcI_right].right_block_left.pos1 - cor_indels[corI].right_block_left.pos1;
					indel_bounds.push_back(make_pair( b1l, b1r ) );

					int g1size = cor_indels[corI].right_block_left.pos1 - cor_indels[corI].left_block_right.pos1 - 1;
					int g2size = absolut(cor_indels[corI].left_block_right.pos2 - cor_indels[corI].right_block_left.pos2) - 1;
					gap_sizes.push_back( make_pair( g1size, g2size ) );

					// the indel was perfectly predicted if no intervening
					// indels were predicted.  record a perfect indel
					perfect_prediction.push_back( calcI_left == calcI_right );
				}
			}
			indel_sp_falseneg += cor_indels.size() - cur_indel_tp;
			indel_sp_falsepos += calc_indels.size() - cur_indel_tp;
			indel_sp_truepos += cur_indel_tp;
		}
	}

	cout << "Indel SP truepos: " << indel_sp_truepos << endl;
	cout << "Indel SP falsepos: " << indel_sp_falsepos << endl;
	cout << "Indel SP falseneg: " << indel_sp_falseneg << endl;
	double sens = 1;
	if( indel_sp_truepos+indel_sp_falseneg > 0 )
		sens = (double)indel_sp_truepos/(double)(indel_sp_truepos+indel_sp_falseneg);
	cout << "Indel SP sensitivity: " << sens << endl;
	double ppv = 1;
	if( indel_sp_truepos+indel_sp_falsepos > 0 )
		ppv = (double)indel_sp_truepos/(double)(indel_sp_truepos+indel_sp_falsepos);
	cout << "Indel SP PPV: " << ppv << endl;

	// write out the indel boundary/size statistics
	ofstream indel_bounds_out( "indel_bound_error_by_gap_size.txt" );
	for( size_t i = 0; i < indel_bounds.size(); i++ )
		indel_bounds_out << indel_bounds[i].first << '\t' << indel_bounds[i].second << '\t' << gap_sizes[i].first << '\t' << gap_sizes[i].second << '\t' << (perfect_prediction[i] ? "p" : "n") << endl;
	indel_bounds_out.close();

	vector< size_t > ib2( indel_bounds.size() * 2 );
	size_t j = 0;
	for( size_t i = 0; i < indel_bounds.size(); i++ )
	{
		ib2[j++] = indel_bounds[i].first;
		ib2[j++] = indel_bounds[i].second;
	}
	std::sort( ib2.begin(), ib2.end() );
	int min = 0;
	int q1 = 0;
	int med = 0;
	int q3 = 0;
	int max = 0;
	double mean = 0;
	double var = 0;
	if(ib2.size() > 0)
	{
		min = ib2.front();
		q1 = ib2[ (size_t)(ib2.size() * 0.25) ];
		med = ib2[ (size_t)(ib2.size() * 0.5) ];
		q3 = ib2[ (size_t)(ib2.size() * 0.75) ];
		max = ib2.back();
		for( size_t i = 0; i < ib2.size(); i++ )
			mean += ib2[i];
		mean /= (double)(ib2.size());
		for( size_t i = 0; i < ib2.size(); i++ )
			var += (mean - (double)ib2[i])*(mean - (double)ib2[i]);
		var /= (double)(ib2.size());
	}
	cout << "Indel boundary prediction summary statistics:\n";
	cout << "min: " << min << endl;
	cout << "q1: " << q1 << endl;
	cout << "median: " << med << endl;
	cout << "q3: " << q3 << endl;
	cout << "max: " << max << endl;
	cout << "mean: " << mean << endl;
	cout << "variance: " << var << endl;
}



/**
 * program to score alignments
 * reads in a "correct" alignment and a calculated alignment
 * scores the calculated alignment based on the correct one
 */
int main( int argc, char* argv[] ){
	
	if( argc < 3 ){
		cout << "scoreAlignment <correct alignment> <calculated alignment> <evolved sequence file> [--disable-lcb-scoring]\n";
		cout << "Use --disable-lcb-scoring when scoring an alignment that may align a given nucleotide more than once.\n";
		return -1;
	}

	boolean score_lcbs = true;
	boolean debug_mismatches = false;	/**< turns on code to debug mismatches in evolved and aligned base pairs */
	string correct_fname = argv[ 1 ];
	string calculated_fname = argv[ 2 ];
	string evolved_fname;
	if( argc > 3 ){
		debug_mismatches = true;
		evolved_fname = argv[ 3 ];
	}

	ifstream correct_in;
	correct_in.open( correct_fname.c_str() );
	if( !correct_in.is_open() ){
		cerr << "Error opening " << correct_fname << endl;
		return -1;
	}
	ifstream calculated_in;
	calculated_in.open( calculated_fname.c_str() );
	if( !calculated_in.is_open() ){
		cerr << "Error opening " << calculated_fname << endl;
		return -1;
	}
	if( argc > 4 )
	{
		string lcb_score_arg = argv[4];
		if( lcb_score_arg == "--disable-lcb-scoring" )
			score_lcbs = false;
	}
//try{
	IntervalList correct_ivs;
	IntervalList calculated_ivs;
	correct_ivs.ReadStandardAlignment( correct_in );
	correct_in.close();
	calculated_ivs.ReadStandardAlignment( calculated_in );
	calculated_in.close();

	if( calculated_ivs.size() == 0 )
	{
		cerr << "WARNING!  No alignment found in input file.  Assuming gap predictions!\n";
		warn_missing = false;
	}

	gnSequence empty_seq;
	vector< gnSequence* > seq_table( correct_ivs[0].SeqCount(), &empty_seq );
	uint seq_count = seq_table.size();
	const gnFilter* comp_filter = gnFilter::DNAComplementFilter();
	
	gnSequence evolved_gnseqs;
	vector< string > evolved_seqs( seq_count );
	if( debug_mismatches ){
		evolved_gnseqs.LoadSource( evolved_fname );
		for( uint i = 0; i < seq_count; i++ ){
			evolved_seqs[ i ] = evolved_gnseqs.contig( i ).ToString();
		}
	}
	
	compareAlignments( correct_ivs, calculated_ivs, evolved_seqs, seq_table );
	if( score_lcbs )
		computeLCBaccuracy( correct_ivs, calculated_ivs, evolved_seqs, seq_table );
	computeIndelAccuracy( correct_ivs, calculated_ivs, evolved_seqs, seq_table );


/*	
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}catch( char const* c ){
	cerr << c << endl;
}catch(...){
	cerr << "Unhandled exception" << endl;
}
*/
}


