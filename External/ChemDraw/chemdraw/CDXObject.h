// CommonCS/LibCommon/Hdr/CDXObject.h
// Contains: Program-independent class library for managing CDX objects
//           This file describes the interface for the generic CDXObject and CDXAttribute classes
// Copyright (c) 1986-2011, PerkinElmer, Inc., All Rights Reserved

// BSD 3-Clause License
// 
// Copyright (c) 1986-2025, CambridgeSoft Corp, Revvity Inc and others.
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#pragma once

#if TARGET_OS_WIN32
// Don't warn about truncation of long symbols
#pragma warning (disable:4786)
// Don't warn about unknown pragmas
#pragma warning (disable:4068)
#endif

//#include "CoreChemistryAPI.h"
#include "CDXIO.h"
#include "cs_swapBytes.h"
#include "Base64Converter.h"
#include <cassert>
#include <functional>
#include <map>
#include <set>
using std::map;
using std::set;
using std::pair;

using std::ostringstream;
using std::istringstream;
using std::stringstream;
using std::ostream;
using std::ofstream;
using std::streamsize;
using std::exception;
using std::nothrow_t;
using std::bad_alloc;
using std::swap;
using std::make_pair;
using std::fill_n;
using std::ios;
using std::multimap;
using std::find;
using std::streambuf;
using std::endl;
using std::hex;
using std::vector;
using std::ostream;
using std::min;
using std::max;
using std::string;

#define CDX_INHERIT_SRC(m_val)	m_val(src.m_val)

// Classes defined in this file:
class CDXObjectFactory;			// Factory class for creating CDXObjects
class CDXObject;				// The abstract base class for objects
class CDXAttribute;				// Base class for complex or unrecognized attributes
class invalid_cdx_error;		// exception thrown when CDX has invalid format of some kind

// Forward defines
class CDXFontTable;

/*
+===========================================================================================+
|																							|
|		CdxObjectsByTag																		|
|																							|
+-------------------------------------------------------------------------------------------+
|	CDX Objects can be referenced by type or by ID.  An object may contain multiple objects	|
|	of the same type, but non-zero IDs should be unique.									|
|	A CDXObjectsByTag object "owns" its objects, that is, they are deleted when the			|
|	CDXObjectsByTag is destructed.															|
+===========================================================================================+
*/
class CDXObjectsByTag : public std::multimap<CDXTag, CDXObject *, std::less<CDXTag> >
{
public:
    CDXObjectsByTag() {}
	CDXObjectsByTag::iterator insertObj(CDXObject *obj); // Insert this object; its tag is obtained from the object
	void WriteTo(CDXDataSink &) const;
	void XMLWrite(XMLDataSink &) const;
	virtual ~CDXObjectsByTag();
};

inline CDXObject * GetObject(CDXObjectsByTag::const_iterator i) { return (*i).second; }


/*
+===========================================================================================+
|																							|
|		CDXObjectsRange																		|
|																							|
+===========================================================================================+
*/
class CORE_CHEMISTRY_API CDXObjectsRange
{
	std::pair<CDXObjectsByTag::const_iterator, CDXObjectsByTag::const_iterator> m_range;
public:
	CDXObjectsRange()	{}
	CDXObjectsRange(const std::pair<CDXObjectsByTag::const_iterator, CDXObjectsByTag::const_iterator> &r) : m_range(r) {}
#if !defined _MSC_VER  ||  _MSC_VER >= 1200
	// Sun Studio compiler has trouble initializing with m_range(r) but accepts its components as (r.first,r.second)
	// -- problem seems to occur only when initialzing const_iterator with iterator
	CDXObjectsRange(const std::pair<CDXObjectsByTag::iterator, CDXObjectsByTag::iterator> r) : m_range(r.first, r.second) {}
#else // VC5 users
	CDXObjectsRange(const std::pair<CDXObjectsByTag::iterator, CDXObjectsByTag::iterator> &r)
			: m_range(std::pair<CDXObjectsByTag::const_iterator, CDXObjectsByTag::const_iterator>(r.first,r.second))	{}
#endif
	CDXObjectsRange(CDXObjectsByTag::const_iterator b, CDXObjectsByTag::const_iterator e) : m_range(b, e) {}
	CDXObjectsByTag::const_iterator begin()	const	{ return m_range.first; }
	CDXObjectsByTag::const_iterator end()	const	{ return m_range.second; }
	int		size() const;
	bool	empty() const	{ return m_range.first == m_range.second; }
};



/*
+===========================================================================================+
|																							|
|		CDXObjectsByID																		|
|																							|
+===========================================================================================+
*/
class CDXObjectsByID : public std::map<CDXObjectID, CDXObject *, std::less<CDXObjectID> >
{
public:
    CDXObjectsByID() {}
	void insertObj(CDXObject *); // Insert this object if it has a non-zero ID; its ID is obtained from the object
};

inline CDXObject * GetObject(CDXObjectsByID::const_iterator i) { return (*i).second; }

/*
+===========================================================================================+
|																							|
|										CdxWalker											|
|																							|
|																							|
+-------------------------------------------------------------------------------------------+
|	Walks all objects in the CDX tree rooted at the supplied node.  The search is depth-	|
|	first.  Penned by HEH, 9/29/05.															|
+===========================================================================================+
*/
class CdxWalker
{
public:
				CdxWalker (const CDXObject& top)	// typically a CDXDocument::
                                : m_top (top)
                                , m_index (0)
                                {}
	virtual		~CdxWalker()		{}

	bool		Next (const CDXObject*& pObj);

private:
	void		Build (const CDXObject* pObj);

	const CDXObject&				m_top;
	std::vector<const CDXObject*>	m_contents;
	int							m_index;	// 0 signifies m_contents hasn't been built.
};

/*
+===========================================================================================+
|																							|
|								CdxGetChildrenIgnoringGroups								|
|																							|
+-------------------------------------------------------------------------------------------+
|	Used to enumerate the children of a given node ignoring Group objects.  That is, if		|
|	a child is a Group, then its children are substituted in place of the Group.			|
|	If a CDXTag is supplied, then like ContainedObjects(), the returned objects will be		|
|	limited to those of the given tag type.  Penned by HEH, 3/17/99.						|
+===========================================================================================+
*/
class CORE_CHEMISTRY_API CdxGetChildrenIgnoringGroups
{
public:
						CdxGetChildrenIgnoringGroups (const CDXObject *pObj, CDXTag tag = 0);
						~CdxGetChildrenIgnoringGroups()	{}
	const CDXObject*	NextChild();
	bool				NextChild (const CDXObject* &pObj)	{ pObj = NextChild(); return (pObj != NULL); }

protected:
	void							MoveIteratorToInterestingItem();
	const CDXObject *				m_curObj;
	CDXTag							m_tag;
	std::vector<const CDXObject*>	m_levelsToTry;	// stack
	CDXObjectsByTag::const_iterator	m_curIt;	// When NextChild() is invoked, guaranteed to point to something interesting, or we're done
};

/*
+===========================================================================================+
|																							|
|		CDXAttributes																		|
|																							|
+-------------------------------------------------------------------------------------------+
|	Although one might have multiple attributes of the same type in a given object,			|
|	I think we're going to permit only a single attribute of a given type, so use a map,	|
|	not a multimap.  A CDXAttributes object "owns" its attributes, that is, they are		|
|	deleted when the CDXAttributes is destructed.											|
+===========================================================================================+
*/
class CDXAttributes : public std::map<CDXTag, CDXAttribute *, std::less<CDXTag> > {
public:
	std::pair<iterator,bool> Insert(CDXAttribute *); // Insert this attribute; its tag is obtained from the object
	bool Remove (CDXAttribute* pAttr);
	void WriteTo(CDXDataSink &) const;
	void XMLWrite(XMLDataSink &) const;
	virtual ~CDXAttributes();
};

inline CDXAttribute * GetAttribute(CDXAttributes::const_iterator i) { return (*i).second; }

// ****************************
// ** class CDXObjectFactory **
// ****************************
//
// This is the class that knows which CDXObject subclass to create
// when a particular object tag is seen.  Generally, only one instance
// of this class ever exists.  You might want to subclass this and
// add suppliers for your CDXObject subclasses in the constructor.

class CORE_CHEMISTRY_API CDXObjectFactory {
public:
	// A pointer to a function that constructs a particular subclass of CDXObject
	typedef CDXObject * (*CDXObjectSupplierProc)(CDXTag tag, CDXObjectID id);

private:
	// This stores the factory 
	typedef std::map<CDXTag, CDXObjectSupplierProc, std::less<CDXTag> > CDXSupplierMap;
	CDXSupplierMap m_supplierMap;

public:
	CDXObjectFactory();
	virtual ~CDXObjectFactory();

	// Add a factory function for a tag represented by a specialized subclass
	void AddFactory(CDXTag t, CDXObjectSupplierProc f);

	// Here's the main entry point for reading a CDX file
	// It reads one complete object and everything it contains.
	CDXObject *ReadOneObject(CDXDataSource &);

	// This is used to allocate the appropriate subclass of CDXObject,
	// once we've read the tag and thus know that some kind of object is coming.
	CDXObject *AllocateObject(CDXObjectID, CDXTag tag);
};


// ************************
// ** class CDXAttribute **
// ************************
//
// This is the abstract base class for complex or unrecognized CDX attributes.
// Simple attributes are stored in specific fields of the object.

class CORE_CHEMISTRY_API CDXAttribute {
protected:
	INT8	*m_data;
	size_t	m_size;
	CDXTag	m_tag;
public:
	// Standard constructors, copy operator, destructor
	CDXAttribute(CDXTag tag, INT8 *data = 0, size_t siz = 0);
#if 0 // If we're going to use these, they must swap bytes when appropriate
	CDXAttribute(CDXTag tag, UINT32 datum) : m_tag(tag), m_size(4), m_data(0) {}
	CDXAttribute(CDXTag tag, UINT16 datum) : m_tag(tag), m_size(2), m_data(0) {}
	CDXAttribute(CDXTag tag, UINT8  datum) : m_tag(tag), m_size(1), m_data(0) {}
#endif
	CDXAttribute(const CDXAttribute &a);
	CDXAttribute & operator=(const CDXAttribute &a);
	~CDXAttribute() { delete [] m_data; }

	CDXTag      GetTag	()			const { return m_tag; }
	size_t		GetSize	()			const { return m_size; }
	const void*	GetData	()			const { return m_data; }
//#if 0 // If we're going to use these, they must swap bytes when appropriate
	INT8		GetData_INT8	()	const	{ assert (GetSize() == 1);  return *( INT8 *) m_data; }
	UINT8		GetData_UINT8	()	const	{ assert (GetSize() == 1);  return *(UINT8 *) m_data; }
	INT16		GetData_INT16	()	const	{ assert (GetSize() == 2);  return SwapBytes (*( INT16*)m_data); }
	UINT16		GetData_UINT16	()	const	{ assert (GetSize() == 2);  return SwapBytes (*(UINT16*)m_data); }
	INT32		GetData_INT32	()	const	{ assert (GetSize() == 4);  return SwapBytes (*( INT32*)m_data); }
	UINT32		GetData_UINT32	()	const	{ assert (GetSize() == 4);  return SwapBytes (*(UINT32*)m_data); }
//#endif
	void		SetData	(INT8 *data, size_t siz) { delete [] m_data; m_data = data; m_size = siz; }

	// Write the contents of the object
	void WriteTo(CDXDataSink &) const;

	// Read the length and data
	void ReadFrom(CDXDataSource &);

	// Write the contents in XML form
	void XMLWrite(XMLDataSink &sink) const;
};

class invalid_cdx_error : public std::exception {
	CDXObjectID m_objectID;
	CDXTag m_objectTag;
	CDXTag m_attribTag;
	std::string m_what;

public:
	invalid_cdx_error(CDXObjectID objID_arg, CDXTag objTag_arg, CDXTag attribTag_arg = 0);
	virtual ~invalid_cdx_error() throw() {}	// because of linux gcc expects same throw characteristiscs as in the base class
	virtual const char* what() const throw() {
		return m_what.c_str();
	}
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//	Use the ff. macro whenever a derived class of CDXObject is declared.  Then implement
//  the copy ctor.  Clone() can be used to copy a CDXObject when its exact type is
//	unknown.  We're not going to allow operator =, so make it private.
#define DECLARE_CDX_OBJECT(CDXType) \
private: \
CDXType& operator=(const CDXType &a); /* not implemented */ \
public: \
CDXType(const CDXType &a); \
virtual	~CDXType(); \
virtual CDXObject*	Clone()	const	{ return NEW CDXType (*this); }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// *********************
// ** class CDXObject **
// *********************
//
// This is the abstract base class for all CDX objects.
// It can read and write itself, and contains methods for
// traversing the object tree.

class CORE_CHEMISTRY_API CDXObject
{
private:

	CDXObjectID 					m_objectID;			// This object's ID
	CDXObjectsByTag*				m_contents;			// All the objects contained in this object
	mutable CDXObjectsByID*			m_contentsByID;		// Those objects contained in this object with non-zero IDs
	CDXAttributes*					m_attributes;		// Attributes not understood by the subclass
	CDXObject*						m_parent;

	static const CDXObjectsByTag	kEmptyContents;
	static const CDXObjectsByID		kEmptyContentsByID;

	friend class CDXObjectFactory;
	friend class XMLParser;

	// Read the contents of this object
	// This is only called from ~CDXObjectFactory::ReadOneObject, right after the
	// object is constructed.
	void ReadFrom(CDXDataSource &, CDXObjectFactory &);

	// Skip the contents of this object
	static void SkipFrom(CDXDataSource &);

	// Read and store the data for one attribute
	void ReadOneAttribute(CDXDataSource &, CDXTag);

	// Skip the data for one attribute
	static void SkipOneAttribute(CDXDataSource &, CDXTag);

	// Contents deletion is not ordinarily needed, as contents are self-destructing.
	void DeleteContents();

	CDXObject& operator=(const CDXObject &a); // assignment is not implemented

protected:
	// Store an attribute
	// This may be overridden by subclasses to store certain attributes in a more efficient
	// subclass-specific way than the m_attributes list.  If the tag represents an attribute
	// that is not stored in a subclass-specific way, then the override should call
	// CDXObject::StoreAttribute to store it in a generic way.
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);

	// Write attributes
	// This may be overridden by subclasses to write attributes which were stored in a more
	// efficient subclass-specific way than the m_attributes list.  An override should call
	// CDXObject::WriteAttributesTo to write any generic attributes that were not stored in 
	// this way.
	virtual void WriteAttributesTo(CDXDataSink &) const;
	
	// FinishReading
	// This may be overridden by subclasses to do anything that needs to be done when
	// the end of the object is encountered.  Some subclasses may want to do validity
	// checking at this point, or perhaps interpret some attributes that can't be
	// interpreted until all the other attributes are read.
	virtual void FinishReading() {}

	// concrete subclasses must implement this to define the XML element name
	virtual std::string XMLObjectName() const = 0;

	// Write attributes
	// This may be overridden by subclasses to write attributes which were stored in a more
	// efficient subclass-specific way than the m_attributes list.  An override should call
	// CDXObject::XMLWriteAttributes to write any generic attributes that were not stored in 
	// this way.
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	// Override this to return true if you need to write content between the start-tag
	// and end-tag (i.e. <node id="1">some content</node>).  If you don't, and there
	// are no contained objects, we'll use an empty-element tag (i.e. <node id="1" />).
	virtual bool XMLNeedToWriteContent() const;

	// Subclasses can override this if they want to write content themselves.
	// Note that contained objects are written automatically - subclasses need
	//    not write them as long as the contained objects are stored in the
	//    normal way (i.e. in m_contents).
	// The default implementation writes nothing.
	virtual void XMLWriteContent(XMLDataSink &) const;

	// Subclasses can override this if they want to do any cleanup after their content is written.
	// This function is called after all properties and subobjects are written, but before the
	// final close tag is emitted.
	// The default implementation does nothing.
	virtual void XMLWriteContentDone(XMLDataSink &) const;

    std::vector<CDXObjectID> ReadObjectIDList(CDXDataSource& dataSource, size_t size);

public:
	// ctor and dtor
	CDXObject(CDXTag tag, CDXObjectID id);
	CDXObject(const CDXObject &a);
	virtual	~CDXObject();
	
	// concrete subclasses must implement this virtual constructor
	virtual CDXObject*	Clone()	const = 0;

	// IsValid() returns True if the object appears sound.  It is not recursive.
	// Although not fully implemented, it is useful in tracking down certain
	// memory corruption problems.
	virtual bool	IsValid() const	{	return true; }

	CDXObjectID				GetObjectID() const		{ return m_objectID; }
	virtual CDXTag			GetTag     () const = 0;	// Returns the object type, like typeid() does.
	CDXObject*				GetParent  ()			{ return m_parent; }
	const CDXObject*		GetParent  () const		{ return m_parent; }
	const CDXAttributes*	GetAttributes()	const	{ return m_attributes; }
	const CDXAttribute*		GetAttribute (CDXDatumID id)	const
								{	if (m_attributes == NULL) return NULL;
										CDXAttributes::const_iterator itAt = m_attributes->find (id);
									if (itAt == m_attributes->end())	return NULL;
									return itAt->second;
								}
	// It is not permissible to change an object's ID if it has a parent referencing that ID.
	void					SetObjectID (CDXObjectID id)	{ assert (GetParent() == NULL); m_objectID = id; }

	// Remove all attributes with the given tag. Return true if something was deleted.
	// (Currently only one instance of a given tag is permitted.)
    bool RemoveAttribute(CDXDatumID id);

	// Write the contents of the object
	void WriteTo(CDXDataSink &) const;

	// Get the raw pointer to the vector of all contained objects (might be NULL!)
	const CDXObjectsByTag *GetContents() const { return m_contents; }

	// Get a range for all contained objects
    CDXObjectsRange ContainedObjects() const;

	// Get a range for contained objects of a given type
    CDXObjectsRange ContainedObjects(CDXDatumID typ) const;

	// Look through attributes for given type
	bool HasAttribute(CDXDatumID id) const
				{ if (m_attributes == NULL) return false; return m_attributes->find(id) != m_attributes->end(); }

	// Add an attribute.  The CDXObject is responsible for disposing of the memory.
    void AddAttribute(CDXTag tag, void *pData, size_t size);
    void AddAttribute(CDXAttribute *pAttrib);

	void GenerateContentsByID() const;

	CDXObject *FindByID(CDXObjectID id) const
				{ GenerateContentsByID();
				if (m_contentsByID == NULL) return NULL;
				CDXObjectsByID::const_iterator i = m_contentsByID->find(id);
				return (i == m_contentsByID->end())? NULL : GetObject(i); }

	CDXObject *FindByIDRecursive(CDXObjectID id) const;

	// Add a newly created child to this object's m_contents and m_contentsByID
	// The child is then owned by this object.  Note that ownership transfers
	// when you call AddChild.  If we fail, we delete it before throwing an exception.
	// Otherwise, it will be deleted in our destructor.
	void AddChild(CDXObject *child);

	// Deletes the given child and removes its reference from the parent's corral
	// of objects.
	void RemoveChild(CDXObject *child, bool andDeleteIt = true);

	// Transfers the given child object from the current parent (this) to a different parent.
	void TransferChildTo(CDXObject *child, CDXObject *destParent);

	// Transfers all the children to a different parent.
	void TransferChildrenTo (CDXObject *destParent);

	friend void CFW51_CDX (const CDXObject *pObj, CDXDataSink &pSink);
	
	// Add some text to this object
    virtual void AppendText(const CDXString& inString) {}
	
	// Figure out the bounds of this object and all its children
	virtual CDXRectangle *BoundsIncludingChildren() const;

	// This writes the xml form of this object.  It uses other (virtual) methods
	// to do the object-specific parts, and calls itself recursively for contained
	// objects.  XMLWrite is the main entry point to write an XML stream.
	void XMLWrite(XMLDataSink &) const;

	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void XMLStoreCharacterData(const std::string &);
    
    // These methods give subclasses an opportunity to prepare text for writing or finalize it
    // after it has been read, using the font table to provide context.
    virtual void FinalizeTextAfterRead(const CDXFontTable& inFontTable)     {}
    virtual void PrepareTextForWrite(const CDXFontTable& inFontTable)       {}
};

class CDXUnknownObject : public CDXObject
{
// This is used when we run into an object of a type which we do not recognize.
public:
	CDXUnknownObject(CDXTag tag, CDXObjectID id) : CDXObject(kCDXObj_UnknownObject, id) { m_unknownTag = tag; }
	virtual CDXObject*	Clone()	const;
	virtual std::string XMLObjectName() const;
	virtual CDXTag	GetTag() const	{ return kCDXObj_UnknownObject; }	// not m_unknownTag!
	// GetTag() is intended to return the true class of the object, like typeid() does.

	CDXTag		GetUnknownsTag() const	{ return m_unknownTag; }

private:
	CDXTag	m_unknownTag;
};

// Used so that sorted lists (like sets) are in objectID order
template<class T>
struct CDXless
{
	inline bool operator()(const T& __X, const T& __Y) const
		{return (__X->GetObjectID() < __Y->GetObjectID()); }
};

CORE_CHEMISTRY_API const string &CDXMLAttributeName(CDXTag);
CORE_CHEMISTRY_API CDXTag CDXMLAttributeID(const char *);

// these functions below with the second parameter an enumerated type are here because of 
// CSBR-49939 CDXToCDXML on Solaris fails exporting arrows

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondOrder  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondDisplay  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondDoublePosition  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondTopology  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondReactionParticipation  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBondCIPType  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXSideType  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXConstraintType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXGeometricFeature  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXGraphicType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBracketType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXRectangleType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXOvalType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXSymbolType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXLineType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXArrowType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXOrbitalType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXFrameType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBracketUsage  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXPolymerRepeatPattern  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXPolymerFlipType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXArrowHeadType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXArrowHeadPosition  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXFillType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXNoGoType  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXNodeType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXRadical  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXRingBondCount  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXUnsaturation  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXReactionStereo  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXTranslation  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXAbundance  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXExternalConnectionType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXAtomGeometry  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXAtomCIPType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXEnhancedStereoType  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXObjectTagType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXPositioningType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXPageDefinition  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXSpectrumYType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXSpectrumXType  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXSpectrumClass  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXTextJustification  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXAminoAcidTermini  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, CDXBioShapeType  v);

CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, double v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, long   v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, int    v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, short  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, char   v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, unsigned long  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, unsigned int   v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, unsigned short v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, unsigned char  v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, bool v);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const std::string &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXPoint2D &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXPoint3D &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXFloat64Point3D &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXFloat64Matrix3D &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const Point32 &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXRectangle &);
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const CDXFontStyle &);

// These translate strings to and from an XML friendly escaped format
CORE_CHEMISTRY_API void XMLPut(XMLDataSink &sink, const std::string &s);
CORE_CHEMISTRY_API void XMLTranslateEscapedString(const std::string& in, std::string& out);
CORE_CHEMISTRY_API void CDXPut(CDXDataSink& sink, const std::string& s);


template <class T>
void XMLPut(XMLDataSink &sink, const std::vector<T> &v)
{
	for (typename std::vector<T>::const_iterator i = v.begin();  i != v.end();  ++i)
	{
		if (i != v.begin()) sink.os << " ";
		XMLPut(sink, *i);
	}
}

template <class T>
void CDXMLPutAttribute(XMLDataSink &sink, CDXTag tag, const T &val)
{
	sink.os << " " << CDXMLAttributeName(tag) << "=\"";
	XMLPut(sink, val);
	sink.os << "\"" << GetTextEOL();
}

void CDXMLPutAttribute(XMLDataSink &sink, CDXTag tag, bool v = true);
void CDXMLPutAttribute(XMLDataSink &sink, vector<CDXChemProp> props);
void CDXMLPutAttributeForObjectIDList(XMLDataSink& sink, CDXTag tag, const std::vector<CDXObjectID>& basisObjects);

std::string HexToBinary(const std::string &s);
std::string BinaryToHex(const std::string &s);

class CDXFontTable;
void CDXPushFontTable(const CDXFontTable *fontTable);
void CDXPopFontTable(const CDXFontTable *fontTable);
