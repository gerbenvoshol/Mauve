#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gnGapSource.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnSourceSpec.h"
#include <cstring>

using namespace std;
using namespace genome;

gnGapSource::gnGapSource()
{
	m_spec = new gnGenomeSpec();
	gnFragmentSpec* fragspec = new gnFragmentSpec();
	gnSourceSpec* sspec = new gnSourceSpec(this);
	m_spec->AddSpec(fragspec);
	fragspec->AddSpec(sspec);
}

gnGapSource* gnGapSource::Clone() const {
	return new gnGapSource( *this );
}

gnGapSource::gnGapSource( const gnGapSource& gnss ) :
m_spec( gnss.m_spec )
{}

gnGapSource& gnGapSource::operator=( const gnGapSource& gnss ) {
	m_spec = gnss.m_spec;
	return *this;
}

void gnGapSource::Close() {}

void gnGapSource::Open( string openString ){};
void gnGapSource::Open(){};
string gnGapSource::GetOpenString() const{ return ""; };
uint32 gnGapSource::GetContigListLength() const{ return 1; }
boolean gnGapSource::HasContig( const string& name ) const{ return true; }
uint32 gnGapSource::GetContigID( const string& name ) const{ return 0; }
string gnGapSource::GetContigName( const uint32 i ) const{ return ""; }

gnSeqI gnGapSource::GetContigSeqLength( const uint32 i ) const {
	return GNSEQI_END / 2;
}


const gnFilter* gnGapSource::GetFilter() const {
	return NULL;
}

void gnGapSource::SetFilter( gnFilter* filter ) {}

boolean gnGapSource::Read( const uint64 pos, char* buf, gnSeqI& bufLen) {
	memset( buf, '-', bufLen );
	return true;
}

boolean gnGapSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ) {
	memset( buf, '-', bufLen );
	return true;
}

gnGenomeSpec *gnGapSource::GetSpec() const{
	return m_spec->Clone();
}

