#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnStringSource_h_
#define _gnStringSource_h_

#include "libGenome/gnBaseSource.h"
#include "libGenome/gnGenomeSpec.h"

#include <string>
#include "libGenome/gnFilter.h"

/**
 * 
 */
class GNDLLEXPORT gnStringSource : public genome::gnBaseSource
{
public:
	gnStringSource( std::string& source );
	~gnStringSource(){}
	gnStringSource( const gnStringSource& gnss );
	gnStringSource& operator=( const gnStringSource& gnss );
	virtual gnStringSource* Clone() const;

	/**
	 * Not implemented
	 */
	virtual void Open( std::string openString );
	/**
	 * Not implemented.
	 */
	virtual void Open();
	/**
	 * Clears the source string.
	 */
	virtual void Close();
	/**
	 * Not implemented
	 */
	virtual std::string GetOpenString() const;

	/**
	 * This always returns 1 because string sources only have a single contig.
	 * @return The number of contigs in this source.
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
	 * Get the total number of base pairs in the specified contig.
	 * @param i The index of the contig or ALL_CONTIGS.
	 * @return The length in base pairs of the specified contig.
	 */
	virtual gnSeqI GetContigSeqLength( const uint32 i ) const;

	/**
	 * Get the filter currently being used to filter unwanted characters out of read sequences.
	 * @return A pointer to the gnFilter currently in use.
	 */
	virtual const genome::gnFilter* GetFilter() const;
	/**
	 * Set the filter that will be used to filter unwanted characters out of the sequence data.
	 * @param filter The filter to remove unwanted characters from the sequence.
	 * @throws NullPointer is thrown if the specified filter pointer is null.
	 */
	virtual void SetFilter( genome::gnFilter* filter );

	virtual boolean Read( const uint64 pos, char* buf, gnSeqI& bufLen);
	virtual boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=genome::ALL_CONTIGS );
	virtual genome::gnGenomeSpec *GetSpec() const;
protected:
	const genome::gnFilter* m_pFilter;
	std::string source_str;
	genome::gnGenomeSpec* m_spec;
private:
};// class gnStringSource

#endif
	// _gnStringSource_h_
