// CommonCS/LibCommon/Src/CDXTest.cpp
// Contains: Test function for CDX class library
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

// This file contains a simple test function which reads a CDX stream that's already
// in memory, and turns around and writes it out to a stdio file.

// Don't warn about truncation of decorated names
#pragma warning (disable:4786)

#include "CDXStdObjects.h"
#include "CDXIO.h"
#include "cs_assert.h"

#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

#define kCDX_HeaderLength 28

void CDXTest(const string &outfile, const UINT8 *p, size_t s);
void CDXTest(const string &outfile, const string &);
void CDXTest(ostream &dest, istream &src);

// CountFragments
//
// This is an example of something you might do with a CDXObject.

class FragmentCounter
{
public:
	int m_numFrags;
	map<INT16, INT16, std::less<INT16> > m_elementCount;
	FragmentCounter() : m_numFrags(0) {}
	void operator()(CDXObject *);
	void operator()(const CDXObjectsByTag::value_type &i) { operator()( i.second ); }
};

void FragmentCounter::operator()( CDXObject *obj )
{
	ASSERT(obj != NULL);
	CDXFragment *frag = dynamic_cast<CDXFragment *> ( obj );
	if (frag == 0)
	{
		// Not a fragment.  Count the fragments in every object contained in this one.
		CDXObjectsRange allObjects = obj->ContainedObjects();
		for_each(allObjects.begin(), allObjects.end(), *this);
	}
	else
	{
		++m_numFrags;
		frag->CountElements(m_elementCount);
	}
}

// CDXTest
// Convert the blob of memory to a CDXObject, then write that out to a file
// Throws a variety of exceptions, including runtime_error and invalid_cdx_error,
// if something goes wrong.
void CDXTest(const string &outfile, const UINT8 *p, size_t s)
{
	ASSERT(p != NULL);
	string str((const char *)p, s);
	ofstream os(outfile.c_str(), ios::out | ios::trunc | ios::binary);
	istringstream is(str, ios::in | ios::binary);
	CDXTest(os, is);
}

void CDXTest(const string &outfile, const string &s)
{
	ofstream os(outfile.c_str(), ios::out | ios::trunc);
	istringstream is(s, ios::in | ios::binary);
	CDXTest(os, is);
}

void CDXTest(ostream &dest, istream &src)
{
	try
	{		
		char hdrBuff[kCDX_HeaderLength];
		src.read(hdrBuff, kCDX_HeaderLength);
		
		// Make a factory; just use the standard classes in CDXStdObjects.cpp
		CDXObjectFactory factory;
		CDXAddStandardFactories(factory);

		// Read one object and everything it contains
		CDXistream is(src);
		unique_ptr<CDXObject> obj ( factory.ReadOneObject(is) );
		// Count the fragments and elements
		FragmentCounter fragCounter;
		fragCounter(obj.get());

		// Print out the result of the count
		printf("Fragments: %d", fragCounter.m_numFrags);
		
		for (map<INT16, INT16, std::less<INT16> >::iterator ifrag = fragCounter.m_elementCount.begin();
			 ifrag != fragCounter.m_elementCount.end();
			 ++ifrag)
			printf("Element %d occurs %d times", (*ifrag).first, (*ifrag).second);

		// Open a file, and dump the objects.  At the moment, CDX files created by
		// ChemDraw don't include an object header for the document object.  Since
		// WriteTo will write this, we subtract 6 from the header length so that
		// the 6-byte document object header doesn't screw up ChemDraw.
		
		//dest.write(hdrBuff, kCDX_HeaderLength-6);
		//(*obj).WriteTo(CDXostream(dest));
		dest << kCDXML_HeaderString;
		XMLDataSink ds(dest);
		(*obj).XMLWrite(ds);
	} // try
	catch (const exception &ex)
	{
		string what = ex.what();
	}
}
