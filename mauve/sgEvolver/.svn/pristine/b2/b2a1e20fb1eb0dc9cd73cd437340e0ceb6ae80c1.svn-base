#include "libMems/IntervalList.h"
#include <string>
#include <fstream>
#include <sstream>

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc < 4 )
	{
		cerr << "Usage: " << "maf2xmfa" << " <maf input file> <xmfa output file> <sequence 1 name> ... <sequence N name>\n";
		return -1;
	}
	int seq_count = argc - 3;

	string maf_in_name = argv[1];
	string xmfa_out_name = argv[2];
	vector< string > seq_filename;
	map< string, uint > seq_name_map;
	for( int i = 3; i < argc; ++i )
	{
		seq_filename.push_back( argv[i] );
		seq_name_map[seq_filename.back()] = i - 3;
	}
	ifstream maf_in( maf_in_name.c_str() );
	if( !maf_in.is_open() )
	{
		cerr << "Error opening \"" << maf_in_name << "\"\n";
		return -2;
	}

	ofstream xmfa_out( xmfa_out_name.c_str() );
	if( !xmfa_out.is_open() )
	{
		cerr << "Error opening \"" << xmfa_out_name << "\"\n";
		return -2;
	}

	cout << "Reading " << maf_in_name << endl;
	// construct an IntervalList from the maf file
	IntervalList iv_list;
	iv_list.seq_filename = seq_filename;
	LoadSequences( iv_list, &cout );

	int read_state = 0;
	string cur_line;
	GappedAlignment gal( seq_count, 0 );
	vector< string > aln_table( seq_count );
	gnSeqI cur_aln_len = 0;
	uint line_count = 0;
	while( getline( maf_in, cur_line ) )
	{
		line_count++;
		if( cur_line == "" )
		{
			if( cur_aln_len > 0 )
			{
				// finish the block
				for( size_t i = 0; i < seq_count; ++i )
				{
					if( aln_table[i] == "" )
						aln_table[i].resize( cur_aln_len, '-' );
				}
				gal.SetAlignment( aln_table );
				vector<AbstractMatch*> asdf( 1, gal.Copy() );
				Interval iv( asdf.begin(), asdf.end() );
				iv_list.push_back( iv );
				// clean up for the next block
				for( size_t i = 0; i < seq_count; ++i )
				{
					aln_table[i].clear();
					gal.SetLeftEnd(i,0);
				}
				cur_aln_len = 0;
			}
			continue;
		}
		if( cur_line[0] == 'a' )
		{
			// score line, do nothing
			continue;
		}
		if( cur_line[0] == 's' )
		{
			// sequence line
			stringstream line_str( cur_line );
			string seq_name;
			gnSeqI lend;
			gnSeqI length;
			string strand;
			gnSeqI total_seqlen;
			string seq;
			line_str >> seq_name;
			line_str >> seq_name;
			line_str >> lend;
			line_str >> length;
			line_str >> strand;
			line_str >> total_seqlen;
			line_str >> seq;
			cur_aln_len = seq.size();
			uint seq_id = seq_name_map[seq_name];
			if(strand == "-"){
				lend = total_seqlen - lend - length + 1;
			}else{
				lend += 1;
			}
			cout << "seq_id: " << seq_id << "\tcur_aln_len: " << cur_aln_len << "\tlend: " << lend + 1 << "\tstrand: " << strand << "\tseq_name: " << seq_name << " glen " << total_seqlen << endl;
			gal.SetLeftEnd( seq_id, lend );
			gal.SetLength( length, seq_id );
			gal.SetOrientation( seq_id, strand == "+" ? AbstractMatch::forward : AbstractMatch::reverse );
			aln_table[seq_id] = seq;
			continue;
		}
	}
	cout << "Read " << line_count << " lines\n";
	cout << "Read " << iv_list.size() << " LCBs\n";
	cout << "Writing output to " << xmfa_out_name << endl;
	try{
		iv_list.WriteStandardAlignment( xmfa_out );
	}catch(gnException& gne){
		cerr << "Error code " << gne.GetCode().GetInt() << " message " << gne.GetMessage() << endl;
	}
}


