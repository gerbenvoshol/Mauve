#include "libMems/IntervalList.h"
#include "libMems/CompactGappedAlignment.h"
#include "libGenome/gnFASSource.h"
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;
using namespace genome;
using namespace mems;

// don't bother writing out evolved sequence contigs shorter than MINIMUM_LENGTH
const int MINIMUM_LENGTH = 100; 

int main( int argc, char* argv[] )
{
	if( argc < 4 )
	{
		cerr << "Usage: " << "breakSimulatedGenomeOnAncestralContigs" << " <ancestral FastA> <evolved sequences> <evolved xmfa file>\n";
		return -1;
	}

	string ancestral_in_name = argv[1];
	string evolved_in_name = argv[2];
	string xmfa_in_name = argv[3];
	vector< string > seq_filename;
	map< string, uint > seq_name_map;

	// read in the alignment
	ifstream xmfa_in;
	xmfa_in.open( xmfa_in_name.c_str() );
	if( !xmfa_in.is_open() ){
		cerr << "Error opening " << xmfa_in_name << endl;
		return -1;
	}
	
	IntervalList evolved_ivs;
	evolved_ivs.ReadStandardAlignmentCompact( xmfa_in );
	xmfa_in.close();
	
	// read in the ancestral FastA
	gnSequence ancestor;
	gnSequence evolved_seqs;
	ancestor.LoadSource(ancestral_in_name);
	evolved_seqs.LoadSource(evolved_in_name);

	// identify the contig breaks
	gnSeqI pos = 0;
	vector< gnSeqI > breaks;
	for(int i = 0; i < ancestor.contigListLength(); i++){
		pos += ancestor.contig(i).length();
		breaks.push_back(pos);
	}

	int seq_count = evolved_seqs.contigListSize() - 1;	
	vector< vector< gnSeqI > > new_breakpoints(seq_count);
	// find the positions in the evolved sequences which correspond to each ancestral contig break
	for(int posI = 0; posI < breaks.size(); posI++){
		vector<gnSeqI> pos;
		vector<bool> column;
		// which alignment block is the current breakpoint in?
		int blockI = 0;
		for(; blockI < evolved_ivs.size(); blockI++){
			if(evolved_ivs[blockI].LeftEnd(seq_count) <= breaks[posI] && 
			   breaks[posI] <= evolved_ivs[blockI].RightEnd(seq_count)){
				break;
			}
		}
		if(breaks[posI] >= evolved_seqs.contig(seq_count).length())
			break;	// we have exceeded the length of the ancestral simulated sequence
		CompactGappedAlignment<>* cga = dynamic_cast< CompactGappedAlignment<>* >(evolved_ivs[blockI].GetMatches()[0]);
		gnSeqI col = cga->SeqPosToColumn( evolved_seqs.contigListSize() - 1, breaks[posI] );
		cga->GetColumn( col, pos, column );
		for( int seqI = 0; seqI < evolved_seqs.contigListSize() - 1; seqI++){
			new_breakpoints[seqI].push_back(pos[seqI]);
		}
	}

	// walk through each evolved sequence, breaking it on the ancestral contig boundaries
	vector< gnSequence > broken_seqs(seq_count);
	for( int seqI = 0; seqI < evolved_seqs.contigListSize() - 1; seqI++){
		new_breakpoints[seqI].push_back(1);
		new_breakpoints[seqI].push_back(evolved_seqs.contigLength(seqI));
		sort(new_breakpoints[seqI].begin(), new_breakpoints[seqI].end());
		for(int breakI=1; breakI < new_breakpoints[seqI].size(); breakI++){
			gnSeqI len = new_breakpoints[seqI][breakI] - new_breakpoints[seqI][breakI-1];
			if(len < MINIMUM_LENGTH) continue;
			broken_seqs[seqI] += evolved_seqs.contig(seqI).subseq(new_breakpoints[seqI][breakI-1],len);
			string contig_name = "contig_";
			contig_name += boost::lexical_cast<string>(breakI);
			broken_seqs[seqI].setContigName(broken_seqs[seqI].contigListSize()-1, contig_name);
		}
		string fasta_fname = evolved_seqs.contigName(seqI) + ".fasta";
		ofstream fasta_out(fasta_fname.c_str());
		for(int i=0; i<broken_seqs[seqI].contigListSize(); i++){
			fasta_out << ">" + broken_seqs[seqI].contigName(i) << endl;
			fasta_out << broken_seqs[seqI].contig(i).ToString() << endl;
		}
		fasta_out.close();
	}
	
}


