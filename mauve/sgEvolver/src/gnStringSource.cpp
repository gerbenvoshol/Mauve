#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gnStringSource.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnSourceSpec.h"
#include <cstring>

using namespace std;
using namespace genome;

gnStringSource::gnStringSource( string& source ) :
m_pFilter( NULL ),
source_str( source )
{
	m_spec = new gnGenomeSpec();
	gnFragmentSpec* fragspec = new gnFragmentSpec();
	gnSourceSpec* sspec = new gnSourceSpec(this);
	m_spec->AddSpec(fragspec);
	fragspec->AddSpec(sspec);
}

gnStringSource* gnStringSource::Clone() const {
	return new gnStringSource( *this );
}

gnStringSource::gnStringSource( const gnStringSource& gnss ) :
m_pFilter( NULL ),
source_str( gnss.source_str ),
m_spec( gnss.m_spec )
{}

gnStringSource& gnStringSource::operator=( const gnStringSource& gnss ) {
	source_str = gnss.source_str;
	m_spec = gnss.m_spec;
	m_pFilter = gnss.m_pFilter;
	return *this;
}

/**
 * Clears the source string.
 */
void gnStringSource::Close() {
	source_str.clear();
}

void gnStringSource::Open( string openString ){};
void gnStringSource::Open(){};
string gnStringSource::GetOpenString() const{ return ""; };
uint32 gnStringSource::GetContigListLength() const{ return 1; }
boolean gnStringSource::HasContig( const string& name ) const{ return true; }
uint32 gnStringSource::GetContigID( const string& name ) const{ return 0; }
string gnStringSource::GetContigName( const uint32 i ) const{ return ""; }

/**
 * Get the total number of base pairs in the specified contig.
 * @param i The index of the contig or ALL_CONTIGS.
 * @return The length in base pairs of the specified contig.
 */
gnSeqI gnStringSource::GetContigSeqLength( const uint32 i ) const {
	return source_str.size();
}


const gnFilter* gnStringSource::GetFilter() const {
	return m_pFilter;
}

void gnStringSource::SetFilter( gnFilter* filter ) {
	m_pFilter = filter;
}

boolean gnStringSource::Read( const uint64 pos, char* buf, gnSeqI& bufLen) {
	if( pos >= source_str.size() )
		return false;
	bufLen = pos + bufLen < source_str.size() ? bufLen : source_str.size() - pos;
	strncpy( buf, source_str.data() + pos, bufLen );
	return true;
}

boolean gnStringSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ) {
	if( m_pFilter == NULL )
		return Read( start, buf, bufLen );
	if( start >= source_str.size() )
		return false;
	gnSeqI bufI = 0;
	for( gnSeqI charI = start; charI < source_str.size(); charI++ ){
		if( m_pFilter->IsValid( source_str[ charI ] ) )
			buf[ bufI++ ] = source_str[ charI ];
		if( bufI == bufLen )
			break;
	}
	bufLen = bufI;
	return true;
}

gnGenomeSpec *gnStringSource::GetSpec() const{
	return m_spec->Clone();
}

