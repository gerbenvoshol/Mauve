#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnGapSource_h_
#define _gnGapSource_h_

#include "libGenome/gnBaseSource.h"
#include "libGenome/gnGenomeSpec.h"

#include <string>
#include "libGenome/gnFilter.h"

/**
 * 
 */
class GNDLLEXPORT gnGapSource : public genome::gnBaseSource
{
public:
	gnGapSource();
	~gnGapSource(){}
	gnGapSource( const gnGapSource& gnss );
	gnGapSource& operator=( const gnGapSource& gnss );
	virtual gnGapSource* Clone() const;

	/**
	 * Not implemented
	 */
	virtual void Open( std::string openString );
	/**
	 * Not implemented.
	 */
	virtual void Open();
	/**
	 * Not implemented.
	 */
	virtual void Close();
	/**
	 * Not implemented
	 */
	virtual std::string GetOpenString() const;

	/**
	 * This always returns 1 because the gap source only has a single contig.
	 * @return 1.
	 */
	virtual uint32 GetContigListLength() const;
	/**
	 * Always returns true
	 */
	virtual boolean HasContig( const std::string& name ) const;
	/**
	 * Always returns 0
	 */
	virtual uint32 GetContigID( const std::string& name ) const;
	/**
	 * Always returns an empty string.
	 */
	virtual std::string GetContigName( const uint32 i ) const;
	/**
	 * Always returns GNSEQI_END, the maximum length of a gap sequence.
	 * @param i The index of the contig or ALL_CONTIGS.
	 * @return GNSEQI_END
	 */
	virtual gnSeqI GetContigSeqLength( const uint32 i ) const;

	/**
	 * gnGapSource doesn't filter anything so it returns NULL.
	 * @return NULL.
	 */
	virtual const genome::gnFilter* GetFilter() const;
	/**
	 * Not implemented
	 */
	virtual void SetFilter( genome::gnFilter* filter );
	
	/**
	 * Returns a string of gap characters of the requested length
	 */
	virtual boolean Read( const uint64 pos, char* buf, gnSeqI& bufLen);
	/**
	 * Returns a string of gap characters of the requested length
	 */
	virtual boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=genome::ALL_CONTIGS );
	virtual genome::gnGenomeSpec *GetSpec() const;
protected:
	genome::gnGenomeSpec* m_spec;
private:
};// class gnGapSource

#endif
	// _gnGapSource_h_
