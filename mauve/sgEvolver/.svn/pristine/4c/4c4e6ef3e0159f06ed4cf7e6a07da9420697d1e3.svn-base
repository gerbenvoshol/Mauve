#include "libGenome/gnSequence.h"
#include "libGenome/gnFASSource.h"
#include <algorithm>

using namespace std;
using namespace genome;

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cerr << "Usage: " << argv[0] << " <mfa name>\n";
		cerr << "One file will be written for each sequence entry in the MFA file\n";
		return -1;
	}
	gnSequence mfa_seq;
	mfa_seq.LoadSource( argv[1] );
	for( uint i = 0; i < mfa_seq.contigListSize(); ++i )
	{
		string ofname = mfa_seq.contigName(i);
		remove( ofname.begin(), ofname.end(), ';' );
		remove( ofname.begin(), ofname.end(), ' ' );
		ofstream ofile( ofname.c_str() );
		if( !ofile.is_open() )
		{
			cerr << "Error opening \"" << ofname << "\"\n";
			return -1;
		}
		gnSequence contig_seq = mfa_seq.contig(i);
		gnFASSource::Write( contig_seq, ofile, false, false );
	}
	return 0;
}

