// CDXObject.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2009 CambridgeSoft Corp., All Rights Reserved

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

// This file contains the implementation for the generic CDXObject,
// CDXAttribute, and CDXObjectFactory classes.  It does not know anything
// about the various kinds of objects and attributes.
//
// It relies on the CDXIO module for low level data transfer.


#if TARGET_OS_WIN32
#pragma warning( disable : 4786 )	// identifier was truncated to 'number' characters in the debug information
#endif

#include "CDXStdObjects.h"
#include "CDXObject.h"
#include "cs_specialchardefs.h"
#include "cs_charUtils.h"
#include "CDXMLNames.h"
#include "XMLParser.h"
#include "cs_assert.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "cs_lock.h"
#include "CDXUnicode.h"
#include "UTF8Iterator.h"

static cs::CriticalSection s_lock;

using namespace std;

// Added by Glysade
namespace cs {
    const char* kSymbolDegree   = u8"\u00B0";  // °
    const char* kSymbolEllipsis = u8"\u2026";  // …
    const char* kSymbolEnDash   = u8"\u2013";  // –
    const char* kSymbolEmDash   = u8"\u2014";  // —
    const char* kSymbolPlusMinus= u8"\u00B1";  // ±
    const char* kSymbolBullet   = u8"\u2022";  // •
    const char* kSymbolCenterDot= u8"\u00B7";  // ·
    const char* kSymbolReg      = u8"\u00AE";  // ®
    const char* kSymbolCopyright= u8"\u00A9";  // ©
    const char* kSymbolAngstrom = u8"\u212B";  // Å
    const char* kSymbolMicro    = u8"\u00B5";  // µ
    const char* kSymbolCent     = u8"\u00A2";  // ¢
    const char* kSymbolPound    = u8"\u00A3";  // £
}
// End added by glysade
/**
 * Fixes invalid utf-8 strings by transcoding illegal one-byte utf-8 encodings to one of the special characters known
 * to ChemDraw, or removing them if no such conversion can be made in a cross-platform way.
 *
 * @param inInternalUTF8String The string to be made safe
 *
 * @return The transformed string
 */
string MakeStringSafe(const string& inInternalUTF8String)
{
    std::string utf8;
    // Reserve space for original string plus a few expansions (1 byte to 2 or 3 bytes)
    utf8.reserve(inInternalUTF8String.size() + 16);

    // Handle our special-case characters. Note that we can't just search and replace
    // as we might replace bytes that are part of a valid UTF-8 stream. Instead we
    // build up a copy of the string using an iterator.
    for (UTF8Iterator i(inInternalUTF8String); !i.AtEnd(); ++i)
    {
        const size_t utf8Length = i.GetUTF8Length();
        bool copyOriginal = true;

        if (utf8Length == 1)
        {
            // If we have a 1 byte UTF-8 length but one of our special multi-byte characters, make sure the
            // multi-byte version goes into the output.

            const auto extraWideChar = i.GetCharacter();
            
            if (extraWideChar >= 127)
            {
                // This is an illegal utf-8 char so we try to fix it or skip it

                copyOriginal = false;
                switch (extraWideChar)
                {
                case kUnicodeCodePointEnDash:
                    utf8 += cs::kSymbolEnDash;
                    break;

                case kUnicodeCodePointEmDash:
                    utf8 += cs::kSymbolEmDash;
                    break;

                case kUnicodeCodePointBullet:
                    utf8 += cs::kSymbolBullet;
                    break;

                case kUnicodeCodePointCenterDot:
                    utf8 += cs::kSymbolCenterDot;
                    break;

                case kUnicodeCodePointDegree:
                    utf8 += cs::kSymbolDegree;
                    break;

                case kUnicodeCodePointEllipsis:
                    utf8 += cs::kSymbolEllipsis;
                    break;

                case kUnicodeCodePointAngstrom:
                    utf8 += cs::kSymbolAngstrom;
                    break;

                case kUnicodeCodePointCent:
                    utf8 += cs::kSymbolCent;
                    break;

                case kUnicodeCodePointPound:
                    utf8 += cs::kSymbolPound;
                    break;

                case kUnicodeCodePointCopyright:
                    utf8 += cs::kSymbolCopyright;
                    break;

                case kUnicodeCodePointReg:
                    utf8 += cs::kSymbolReg;
                    break;

                case kUnicodeCodePointPlusMinus:
                    utf8 += cs::kSymbolPlusMinus;
                    break;

                case kUnicodeCodePointMicro:
                    utf8 += cs::kSymbolMicro;
                    break;

                case kUnicodeCodePointReplacement:
                default:
                    // Reject/excise this character as unknown, invalid, or unstable between Mac and Windows since transcoding result in differences
                    break;
                }
            }
        }

        if (copyOriginal)
        {
            // Copy the original data over to the new string
            utf8 += inInternalUTF8String.substr(i.GetByteIndex(), utf8Length);
        }
    }

    return utf8;
}


namespace
{
    /**
     * Determine if a utf8 character is legal in XML.
     * XML 1.0 disallows all ascii control characters except for \t \r and \n.
     * XML 1.1 allows for all but NUL if escaped using the &#x; construction but this does not work with our XMLParser.
     *
     * @param  wchar the utf8 character to be tested
     *
     * @return whether character is allowed in output
     */
    bool IsLegalXML(char32_t wchar)
    {
        const char32_t tab = 9;
        const char32_t lf = 10;
        const char32_t cr = 13;
        const char32_t space = 32;

        const bool illegal = (wchar < space) && (wchar != tab) && (wchar != lf) && (wchar != cr);

        return !illegal;
    }

    /**
     * Removes any illegal characters from string. Although XML supports utf8, some ascii characters are forbidden.
     * See IsLegalXML() for details. This will also truncate the string at the first null.
     *
     * @param utf8Source the input string
     *
     * @return the input string with any illegal characters removed
     */
    string RemoveIllegalXMLChars(const string& utf8Source)
    {
      return utf8Source;
      /* XXX We don't have IsLegalXML nor UTF8Iterator
        string result;
        result.reserve(utf8Source.size());

        UTF8Iterator iter(utf8Source);
        for (; !iter.AtEnd(); iter.Next())
        {
            // Filter out illegal control characters
            if (IsLegalXML(iter.GetCharacter()))
            {
                result += utf8Source.substr(iter.GetByteIndex(), iter.GetUTF8Length());
            }
        }

        return result;
       */
    }
}

// ************************
// ** class CDXAttribute **
// ************************

// class CDXAttribute
// Copy constructor

CDXAttribute::CDXAttribute(const CDXAttribute &a)
 : m_tag(a.m_tag), m_size(a.m_size), m_data(NULL)
{
	if (a.m_data != NULL)
	{
		m_data = NEW INT8[m_size];
		memcpy(m_data, a.m_data, m_size);
	}
}

CDXAttribute::CDXAttribute (CDXTag tag, INT8 *data /* = NULL */, size_t siz /* = 0 */)
	:	m_tag(tag)
	,	m_size(siz)
	,	m_data(data)	//AMT-Note: If specified, this parameter must be allocated with: new <type>[]
{
}

// class CDXAttribute
// Copy operator

CDXAttribute &
CDXAttribute::operator= (const CDXAttribute &a)
{
	m_tag = a.m_tag;
	m_size = a.m_size;

	if (a.m_data != NULL)
	{
		INT8 *newData = NEW INT8[a.m_size];
		INT8 *saveData = m_data;
		m_data = newData;
		delete [] saveData;
	}
	else
	{
		delete [] m_data;
		m_data = NULL;
	}

	return *this;
}

// CDXAttribute::WriteTo
// Write this attribute to the data sink

void CDXAttribute::WriteTo(CDXDataSink &sink) const
{
	// LOGPRINTF(("CDXAttribute::WriteTo tag %4.4x size %4d", m_tag, m_size));

	// Zero is an invalid tag value
	ASSERT(m_tag != 0);
	if (m_tag == 0)
		return;

	// First, the tag which describes what kind of attribute this is.
	sink.PutTag( m_tag );
	
	// Next comes the size.  If this is greater than 0xFFFE, then we
	// write the special value, 0xFFFF, which signals that the real length
	// follows.
	if (m_size <= 0xFFFE)
		sink.Put( UINT16(m_size) );
	else
	{
		sink.Put( UINT16(0xFFFF) );
		sink.Put( UINT32(m_size) );
	}

	// Finally, the data.
	if (m_data != NULL)
		sink.Put((UINT8 *)m_data, m_size);
}

inline char ToHex(unsigned char c)
{
	return c >= 10 ? (c - 10 + 'A') : (c + '0');
}

void CDXAttribute::XMLWrite(XMLDataSink &sink) const
{
	// This is used for tags which we don't understand.
	// This should never occur in XML, but we want to know about them, so write something.

	sink.os << " " << CDXMLAttributeName(m_tag) << "=\"|x|";
	unsigned char *p;
	size_t i;
	for (p = (unsigned char *)GetData(), i = GetSize();  i > 0;  --i, ++p)
		sink.os << ToHex(*p >> 4) << ToHex(*p & 0xF);
	sink.os << "\"" << GetTextEOL();
}

// ****************************
// ** class CDXObjectFactory **
// ****************************

// class CDXObjectFactory
// constructor

CDXObjectFactory::CDXObjectFactory()
{
}

// class CDXObjectFactory
// destructor

CDXObjectFactory::~CDXObjectFactory()
{
}

// CDXObjectFactory::ReadOneObject
// Read an object, and the objects and attributes contained in it, from the data source.
//
// This is the main entry point that outsiders will call to read
// a CDX object.  If the first thing in the data stream is an
// attribute tag, we pretend there was a document object there.

CDXObject *CDXObjectFactory::ReadOneObject(CDXDataSource &src)
{
	unique_ptr<CDXObject> obj;

	CDXTag tag = src.GetTag();
	if (!CDXIsObjectTag(tag))
	{
		// Pretend there was a document object
		unique_ptr<CDXObject> w(AllocateObject(0, kCDXObj_Document));
		obj = std::move(w);
		
		// Take care of that first attribute
		(*obj).ReadOneAttribute(src, tag);
	}
	else
	{
		// If the stream started with an object tag, then use the internal
		// method to allocate the right kind of object.
		unique_ptr<CDXObject> w( AllocateObject(src.GetObjectID(), tag) );
		obj = std::move(w);
	}

	(*obj).ReadFrom(src, *this);
    (*obj).FinishReading();
    
	return obj.release();
}

// CDXObjectFactory::AddFactory
//
// Add a factory function for a tag represented by a specialized subclass.
// If the key is repeated, new factory will replace the old one.

void CDXObjectFactory::AddFactory(CDXTag t, CDXObjectSupplierProc f)
{
	CDXSupplierMap::iterator pos = m_supplierMap.find(t);
	if (pos != m_supplierMap.end())
		m_supplierMap.erase(pos);

	m_supplierMap.insert( CDXSupplierMap::value_type (t, f) ); 
}

// CDXObjectFactory::AllocateObject
//
// This is used internally to allocate the appropriate subclass of CDXObject, once we've read the
// tag and thus know that some kind of object is coming.

CDXObject *CDXObjectFactory::AllocateObject(CDXObjectID id, CDXTag tag)
{
	CDXObject *obj;

	// Check to see if a specialized factory function for this tag has been registered.
	CDXSupplierMap::iterator i = m_supplierMap.find(tag);
	if (i != m_supplierMap.end())
	{
		// There's a specialized subclass for this tag.  Use its factory to create it.
		obj = (* (*i).second) (tag, id);
	}
	else	// user-defined object, or a std object we haven't provided a class for
	{
		obj = new CDXUnknownObject(tag, id);
	}

	return obj;
}

// *********************
// ** class CDXObject **
// *********************

CDXObject::CDXObject (const CDXObject &src)	// Clone the object
	:	m_objectID	(src.m_objectID)
	,	m_contents(NULL)
	,	m_contentsByID(NULL)
//	,	m_tag		(src.GetTag())
	,	m_parent	(NULL)
	,	m_attributes(NULL)
{
	if (src.m_contents != NULL)
	{
		m_contents = new CDXObjectsByTag;
		// Clone the contained objects
		for (CDXObjectsByTag::const_iterator  itTag = src.m_contents->begin();  itTag != src.m_contents->end();  ++itTag)
		{
			const CDXObject*	pChild = itTag->second;
			AddChild (pChild->Clone());	// recursive
		}
	}
	// Clone the contained attributes
	if (src.m_attributes != NULL)
	{
		for (CDXAttributes::const_iterator itAttr = src.m_attributes->begin();  itAttr != src.m_attributes->end();  ++itAttr)
		{
			const CDXAttribute	*pAttr = itAttr->second;
			AddAttribute (new CDXAttribute (*pAttr));
		}
	}
}

bool CDXObject::RemoveAttribute(CDXDatumID id) 
{ 
    if (m_attributes == NULL)
    {
        return false;
    }
    
    return m_attributes->Remove(const_cast<CDXAttribute*>(GetAttribute(id)));
}

void CDXObject::AddAttribute(CDXTag tag, void *pData, size_t size)
{
    if (m_attributes == NULL)
    {
        m_attributes = new CDXAttributes;
    }
    
    m_attributes->Insert(new CDXAttribute(tag, (INT8 *)pData, size));
}
void CDXObject::AddAttribute(CDXAttribute *pAttrib)
{
    if (m_attributes == NULL)
    {
        m_attributes = new CDXAttributes;
    }
    
    m_attributes->Insert(pAttrib);
}

// CDXObject destructor
// The contained objects and attributes are self-deleting.

CDXObject::~CDXObject()
{
	DeleteAndNull(m_contents);
	DeleteAndNull(m_contentsByID);
	DeleteAndNull(m_attributes);
}


const CDXObjectsByTag	CDXObject::kEmptyContents = CDXObjectsByTag();
const CDXObjectsByID	CDXObject::kEmptyContentsByID = CDXObjectsByID();


// Get a range for all contained objects
CDXObjectsRange CDXObject::ContainedObjects() const
{
    if (m_contents != NULL)
    {
        return CDXObjectsRange(m_contents->begin(), m_contents->end());
    }
    return CDXObjectsRange(kEmptyContents.begin(), kEmptyContents.end());
}

// Get a range for contained objects of a given type
CDXObjectsRange CDXObject::ContainedObjects(CDXDatumID typ) const
{
    if (m_contents != NULL)
    {
        return m_contents->equal_range(typ);
    }
    return CDXObjectsRange(kEmptyContents.begin(), kEmptyContents.end());
}

// CDXObject::StoreAttribute
//
// This stores the attributes in the m_attributes vector.

void CDXObject::StoreAttribute(CDXDataSource &src, CDXTag t, size_t s)
{
	/*
	+---------------------------------------------------------------------------------------+
	|	If this attribute is not user-defined, either it is bad, or we (the CDX reader) are	|
	|	broken.  All attributes which are not user-defined must be handled by the subclass,	|
	|	so XML comes out right.																|
	+---------------------------------------------------------------------------------------+
	*/

	// Silence this ASSERT so CF databases don't have lots of ASSERTs popping
	// up during debug. Should be turned back on once v11 is released.
//	ASSERT((t & kCDXTag_UserDefined) != 0);
	unique_ptr<CDXAttribute> a;
	
	if (s == 0) {
		unique_ptr<CDXAttribute> w( NEW CDXAttribute( t ) );
		a = std::move(w);
	}
	else
	{
		// Store an uninterpreted attribute.
		// Make a new, empty CDXAttribute
		unique_ptr<CDXAttribute> w( NEW CDXAttribute(t) );
		a = std::move(w);
		// Allocate some storage for the value, pass ownership to the attribute
		char *p = new char[s];
		a.get()->SetData((INT8 *)p, s);
		// Read the data
		src.GetBytes(p, s);
	}

	// Store the attribute in our attribute vector
	AddAttribute (a.get());
	a.release();
}

// CDXObject::WriteAttributesTo
//
// This writes any attributes which are not stored in m_attributes.  The base class has
// none; this is to be overridden by subclasses if they have any.

void CDXObject::WriteAttributesTo(CDXDataSink &sink) const
{
	if (m_attributes != NULL)
		m_attributes->WriteTo(sink);			// The contained attributes
}

// CDXObject::XMLWriteAttributes
//
// This writes any attributes which are not stored in m_attributes.  The base class has
// none; this is to be overridden by subclasses if they have any.

void CDXObject::XMLWriteAttributes(XMLDataSink &sink) const
{
	if (m_attributes != NULL)
		m_attributes->XMLWrite(sink);			// The contained attributes
}

string CDXObject::XMLObjectName() const
{
	ostringstream os;
	os << "object" << GetTag();
	return os.str();
}

// CDXObject::ReadFrom
//
// Read the next object from a CDXDataSource.  This figures out what kind of object is
// coming, and then reads the entire object, including its attributes and contained objects.

void CDXObject::ReadFrom(CDXDataSource &src, CDXObjectFactory &factory)
{
	CDXTag tag;
	while ((tag = src.GetTag()) != kCDXProp_EndObject)
	{
		if (src.IgnoreObject(tag, this))
			SkipFrom(src);
		else if (!CDXIsObjectTag(tag))
			ReadOneAttribute(src, tag);
		else
		{
			// Allocate the right kind of CDXObject
			unique_ptr<CDXObject> obj(factory.AllocateObject(src.GetObjectID(), tag));

			// Read the contents of the object
			(*obj).ReadFrom(src, factory);

			(*obj).FinishReading();

			AddChild(obj.release()); // It's now owned by our m_contents
		}
	}
}

// CDXObject::SkipFrom
//
// Skip the next object from a CDXDataSource.  This figures out what kind of object is
// coming, and then skips the entire object, including its attributes and contained objects.

void CDXObject::SkipFrom(CDXDataSource &src)
{
	src.GetObjectID();

	CDXTag tag;
	while ((tag = src.GetTag()) != kCDXProp_EndObject)
	{
		if (!CDXIsObjectTag(tag))
			SkipOneAttribute(src, tag);
		else
			SkipFrom(src);
	}
}

// CDXObject::ReadOneAttribute
//
// Read the rest of an attribute, and store it in this object.  This is used
// internally; subclasses should override StoreAttribute, not ReadOneAttribute.

void CDXObject::ReadOneAttribute(CDXDataSource &src, CDXTag tag)
{
	ASSERT (!CDXIsObjectTag(tag));

	// Read the length of the attribute
	UINT32 len = src.GetUINT16();
#ifdef _HEH_
	ASSERT (len < 10000);
#endif

	// An attribute larger than 0xFFFE stores a "length" of 0xFFFF followed by
	// the real length as a UINT32
	if (len == 0xFFFF)
		len = src.GetUINT32();

	// Data corruption will frequently show up here.  It's better to trap
	// it now than to endure an access violation down the road.
	if (len > 0x10000000)
		throw std::runtime_error ("CDX data format error");

	if (src.IgnoreTag(tag, this))
		src.SkipBytes(len);
	else
		StoreAttribute(src, tag, len);	
}

// CDXObject::SkipOneAttribute
//
// Skip the rest of an attribute.

void CDXObject::SkipOneAttribute(CDXDataSource &src, CDXTag tag)
{
	ASSERT (!CDXIsObjectTag(tag));

	// Read the length of the attribute
	UINT32 len = src.GetUINT16();
#ifdef _HEH_
	ASSERT (len < 10000);
#endif

	// An attribute larger than 0xFFFE stores a "length" of 0xFFFF followed by
	// the real length as a UINT32
	if (len == 0xFFFF)
		len = src.GetUINT32();

	// Data corruption will frequently show up here.  It's better to trap
	// it now than to endure an access violation down the road.
	if (len > 0x10000000)
		throw std::runtime_error ("CDX data format error");

	src.SkipBytes(len);
}

// CDXObject::WriteTo
//
// Write out an entire object, its attributes, and contained objects.
// This should not be overridden.  Subclasses can override WriteAttributesTo
// in order to store attributes that are handled by the subclass.

void CDXObject::WriteTo(CDXDataSink &sink) const
{
	// LOGPRINTF(("CDXObject::WriteTo tag %4.4x ID %4d; %d attribs, %d subobjects",
	//			GetTag(), m_objectID, m_attributes.size(), m_contents.size()));

	// An object consists of the tag and ID, the attributes,
	// the contained objects, and an EndObject tag

	// Normally the tag is the object's intrinsic type.  There is one exception:
	// an object we did not recognize.
	CDXTag	tag = GetTag();
	if (tag == kCDXObj_UnknownObject)
		tag = ((CDXUnknownObject*)(this))->GetUnknownsTag();
	sink.PutTag(GetTag());
	sink.PutObjectID(m_objectID);
	WriteAttributesTo(sink);			// Attributes handled by the subclass
	if (m_contents != NULL)
		m_contents->WriteTo(sink);			// The contained objects
	sink.PutTag(kCDXProp_EndObject);
}

// CDXObject::XMLWrite
//
// This writes the xml form of this object.  It uses other (virtual) methods
// to do the object-specific parts, and calls itself recursively for contained
// objects.  XMLWrite is the main entry point to write an XML stream.
void CDXObject::XMLWrite(XMLDataSink &sink) const
{
	// First write the opening tag and the attributes.
	sink.os << "<" << XMLObjectName() << GetTextEOL();
	// The id is the only totally generic tag
	if (m_objectID != 0)
		sink.os << " " << kCDXML_id << "=\"" << m_objectID << "\"" << GetTextEOL();

	// This is overridden by subclasses to write the object-specific data
	XMLWriteAttributes(sink);
	
	// If there's any 
	if (!XMLNeedToWriteContent() && (m_contents == NULL || m_contents->empty()))
		sink.os << "/>";
	else
	{
		sink.os << ">";
		XMLWriteContent(sink);
		if (m_contents != NULL)
			m_contents->XMLWrite(sink);			// The contained objects
		XMLWriteContentDone(sink);
		sink.os << "</" << XMLObjectName() << ">";
	}
}

// CDXObject::XMLNeedToWriteContent
//
// Override this to return true if you need to write content between the start-tag
// and end-tag (i.e. <node id="1">some content</node>).  If you don't, and there
// are no contained objects, we'll use an empty-element tag (i.e. <node id="1" />).
bool CDXObject::XMLNeedToWriteContent() const
{
	return false; // Unless we're overridden, an empty-element tag is OK
}

// CDXObject::XMLWriteContent
//
// Subclasses can override this if they want to write content themselves.
// Note that contained objects are written automatically - subclasses need
//    not write them as long as the contained objects are stored in the
//    normal way (i.e. in m_contents).
// The default implementation writes nothing.
void CDXObject::XMLWriteContent(XMLDataSink &) const
{
}

// CDXObject::XMLWriteContentDone
//
// Subclasses can override this if they want to do any cleanup after their content is written.
// This function is called after all properties and subobjects are written, but before the
// final close tag is emitted.
// The default implementation does nothing.
void CDXObject::XMLWriteContentDone(XMLDataSink &) const
{
}

inline UINT8 fromhex(char c)
{
	return (c >= 'a') ? (c - 'a' + 10) : ((c >= 'A') ? (c - 'A' + 10) : (c - '0'));
}

inline char tohex(UINT8 c)
{
	return (c > 9) ? (c + 'A' - 10) : (c + '0');
}

string HexToBinary(const string &s)
{
	if ((s.length() & 1) != 0) // odd length is an error
		throw std::runtime_error("hex strings must have an even number of characters");
	string result(s.length() / 2, '\0');	// Initialize with the right number of nulls
	string::iterator ri = result.begin();
	for (string::const_iterator si = s.begin();  si != s.end();  ++ri, si += 2)
		*ri = (fromhex(si[0]) << 4) | fromhex(si[1]);
	return result;
}

string BinaryToHex(const string &s)
{
	string result(s.length() * 2, '\0');	// Initialize with the right number of nulls
	string::iterator ri = result.begin();
	for (string::const_iterator si = s.begin();  si != s.end();  ri += 2, si += 1)
	{
		ri[0] = tohex(UINT8(si[0]) >> 4);
		ri[1] = tohex(UINT8(si[0]) & 0xF);
	}
	return result;
}

void CDXObject::XMLStoreAttribute(CDXTag t, const string &s)
{
	const string *result;
	string binstr;
	if (s.substr(0,3) == "|x|")
	{
		binstr = HexToBinary(s.substr(3, string::npos));
		result = &binstr;
	}
	else
	{
//		ASSERT (false); // We probably don't do the right thing with non-hex attributes
		result = &s;
	}
	UINT8 *p = new UINT8[result->size()];
#pragma warning (disable:4996)	// "'std::copy': Function call with parameters that may be unsafe."
	std::copy(result->begin(), result->end(), p);
	AddAttribute (new CDXAttribute(t, (INT8 *)p, result->size()));
}

void CDXObject::XMLStoreCharacterData(const string &)
{
}

// CDXObject constructor
//
// Construct the object, and add it to its parent

CDXObject::CDXObject(CDXTag , CDXObjectID id)
	: m_objectID(id)
	, m_contents(NULL)
	, m_contentsByID(NULL)
//	, m_tag(tag)
	, m_attributes(NULL)
	, m_parent(NULL)
{
}
//-----------------------------------------
void CDXObject::AddChild(CDXObject* child)
{
    if (!child)
    {
        ASSERT(false);
        return;
    }

    // Child's affiliation with a different parent should have been broken by now.
    ASSERT(child->m_parent == NULL);

    // If it has an ID, add it to the ObjectsByID map
    if (m_contentsByID && child->GetObjectID())
    {
        m_contentsByID->insertObj(child);
    }

    // Add it to the parent's ObjectsByTag multimap
    if (!m_contents)
    {
        m_contents = new CDXObjectsByTag;
    }

    m_contents->insertObj(child);
    child->m_parent = this;
}

//-----------------------------------------
void CDXObject::RemoveChild(CDXObject *pChild, bool andDeleteIt /* = true */)
{
	ASSERT (pChild->GetParent() == this);
	if (m_contents != NULL)
	{
		ASSERT (m_contents->find (pChild->GetTag()) != m_contents->end());
		pair< CDXObjectsByTag::iterator, CDXObjectsByTag::iterator> p = m_contents->equal_range(pChild->GetTag());
		for (CDXObjectsByTag::iterator it = p.first; it != p.second; ++it)
			if (it->second == pChild)
			{
				m_contents->erase(it);
				break;
			}
	}
	if (m_contentsByID != NULL)
	{
		CDXObjectsByID::iterator objsIt = m_contentsByID->find (pChild->GetObjectID());
		if (objsIt != m_contentsByID->end())
			m_contentsByID->erase (objsIt);
	}
	if (andDeleteIt)
		delete pChild;
	else
		pChild->m_parent = NULL;
}
//-----------------------------------------
void CDXObject::TransferChildTo(CDXObject *pChild, CDXObject *destParent)
{
	RemoveChild (pChild, false);
	destParent->AddChild (pChild);
}
//-----------------------------------------
void CDXObject::TransferChildrenTo(CDXObject *destParent)
{
	if (m_contents == NULL)
		return;

	while (!m_contents->empty())
	{
		CDXObject*	pChild = m_contents->begin()->second;
		TransferChildTo (pChild, destParent);
	}
}
//-----------------------------------------
void CDXObject::GenerateContentsByID() const
{
	if (m_contentsByID == NULL && m_contents != NULL)
	{
		if (m_contentsByID == NULL)
			m_contentsByID = new CDXObjectsByID;
		for (CDXObjectsByTag::const_iterator it = m_contents->begin(); it != m_contents->end(); ++it)
			m_contentsByID->insertObj( it->second );
	}
}
//-----------------------------------------
CDXObject *CDXObject::FindByIDRecursive(CDXObjectID id) const
{
	CDXObject *retVal = FindByID(id);
	if (retVal != NULL)
		return retVal;

	if (m_contents != NULL)
		for (CDXObjectsByTag::const_iterator it = m_contents->begin(); retVal == NULL && it != m_contents->end(); ++it)
			retVal = it->second->FindByIDRecursive(id);

	return retVal;
}

vector<CDXObjectID> CDXObject::ReadObjectIDList(CDXDataSource& dataSource, size_t size)
{
    vector<CDXObjectID> objectIDList;

    if ((size % sizeof(CDXObjectID)) != 0)
    {
        throw std::runtime_error("Object list attribute value is invalid");
    }

    size_t numItems = size / sizeof(CDXObjectID);
    while (numItems > 0)
    {
        numItems--;
        const CDXObjectID objectID = dataSource.GetUINT32();
        if (objectID < 0)
        {
            continue;
        }

        objectIDList.push_back(objectID);
    }

    return objectIDList;
}

// ***************************
// ** class CDXObjectsRange **
// ***************************
int CDXObjectsRange::size() const
{
	const CDXObjectsRange&	range = *this;
#ifdef __solaris
	// Alternative specific to Solaris: 	int siz; distance (range.begin(), range.end(), siz); return siz;
	CDXObjectsByTag::const_iterator itBegin = range.begin();
	CDXObjectsByTag::const_iterator itEnd = range.end();
	int n = 0;
	for (; itBegin!=itEnd; ++itBegin)
		++n;
	return n;
#else
	return (int)distance (range.begin(), range.end());
#endif
}


// ***************************
// ** class CDXObjectsByTag **
// ***************************

// CDXObjectsByTag is a collection of CDXObject's.
//
// CDX Objects can be referenced by type or by ID.  An object may contain multiple objects
// of the same type, but non-zero IDs should be unique.
// A CDXObjectsByTag object "owns" its objects, that is, they are deleted when
// the CDXObjectsByTag is destructed.

// CDXObjectsByTag::insert
//
// Insert an object into a CDXObjectsByTag map.

CDXObjectsByTag::iterator
CDXObjectsByTag::insertObj(CDXObject *obj)
{
//	compiler complains about this override when using std namespace
//	return multimap<CDXTag, CDXObject *, std::less<CDXTag> >::insert( value_type(obj->GetTag(), obj) FILELINE );

	return insert( value_type(obj->GetTag(), obj) );
}

// CDXObjectsByTag::WriteTo
//
// Write all of the objects in the collection

void CDXObjectsByTag::WriteTo(CDXDataSink &sink) const
{
	for (const_iterator i = begin();  i != end();  ++i)
		GetObject(i)->WriteTo(sink);
}

void CDXObjectsByTag::XMLWrite(XMLDataSink &sink) const
{
	for (const_iterator i = begin();  i != end();  ++i)
		GetObject(i)->XMLWrite(sink);
}

// CDXObjectsByTag::~CDXObjectsByTag
//
// delete all of the objects in the collection

CDXObjectsByTag::~CDXObjectsByTag()
{
	for (const_iterator i = begin();  i != end();  ++i)
		delete GetObject(i);
}

// **************************
// ** class CDXObjectsByID **
// **************************

// CDXObjectsByID::insert
//
// Insert an object into a CDXObjectsByID map.

void
CDXObjectsByID::insertObj(CDXObject *obj)
{
	CDXObjectID objID = obj->GetObjectID();

	if (objID != 0)
	{
		/*pair<iterator,bool> result = */insert( value_type(objID, obj) );
//		if (!result.second)
//			throw std::runtime_error("Multiple objects with same ID");
	}
}



// ****************************************
// *********** class CdxWalker ************
// ****************************************

bool CdxWalker::Next (const CDXObject*& pObj)
{
	if (!m_index)
		Build (&m_top);
	if (m_index < m_contents.size())
	{
		pObj = m_contents [m_index++];
		return true;
	}
	m_index = 0;	// so they can restart
	return false;
}
//-----------------------------------------------------------
void CdxWalker::Build (const CDXObject* pObj)
{
	m_contents.push_back (pObj);
	CDXObjectsRange	children = pObj->ContainedObjects();
	for (CDXObjectsByTag::const_iterator it = children.begin();  it != children.end();  it++)
		Build (it->second);
}


// ****************************************
// ** class CdxGetChildrenIgnoringGroups **
// ****************************************

//----------------------------------------------------------
CdxGetChildrenIgnoringGroups::CdxGetChildrenIgnoringGroups (const CDXObject *pObj, CDXTag tag /* = 0 */)
	: m_curObj(NULL)
	, m_tag	(tag)
{
	m_levelsToTry.push_back(pObj);
	MoveIteratorToInterestingItem();
}
//----------------------------------------------------------
void CdxGetChildrenIgnoringGroups::MoveIteratorToInterestingItem()
{
	do
	{
		if (m_curObj == NULL)
		{
			if (m_levelsToTry.empty())
				return;

			m_curObj = m_levelsToTry.back();
			m_levelsToTry.pop_back();

			ASSERT(m_curObj != NULL);
			if (m_curObj == NULL)
				continue;

			if (m_curObj->GetContents() == NULL)
			{
				m_curObj = NULL;
				continue;
			}

			m_curIt = m_curObj->GetContents()->begin();
		}

		const CDXObjectsByTag *curChildren = m_curObj->GetContents();
		if (curChildren != NULL)
		{
			while (m_curIt != curChildren->end())
			{
				const CDXDatumID id = (CDXDatumID)m_curIt->first;
				if (id == m_tag)
					return;	// Could include fragment, and perhaps even group or page.
				if (id == kCDXObj_Group  ||  id == kCDXObj_Fragment  ||  id == kCDXObj_Page)
					m_levelsToTry.push_back (m_curIt->second);
				else if (m_tag == kCDXProp_EndObject)
					return;
				++m_curIt;
			}
		}

		m_curObj = NULL;
	} while (!m_levelsToTry.empty());
}
//----------------------------------------------------------
const CDXObject*  CdxGetChildrenIgnoringGroups::NextChild()
{
	if (m_curObj == NULL)
		return NULL;

	const CDXObjectsByTag *curChildren = m_curObj->GetContents();
	if (m_curIt == curChildren->end())
		return NULL;

	const CDXObject	*pResult = (m_curIt++)->second;
	MoveIteratorToInterestingItem();
	return pResult;
}
//----------------------------------------------------------



// Figure out the bounds of this object and all its children
CDXRectangle *CDXObject::BoundsIncludingChildren() const
{
	CDXRectangle *r = NULL;
	CDXObjectsRange subObjs = ContainedObjects();
	for (CDXObjectsByTag::const_iterator  subObj = subObjs.begin();  subObj != subObjs.end();  subObj++)
	{
		unique_ptr<CDXRectangle> subBounds(subObj->second->BoundsIncludingChildren());
		if (subBounds.get() != NULL)
		{
			if (r == NULL)
			{
				r = subBounds.release();
			}
			else
			{
				r->left = ((r->left < subBounds->left) ? r->left : subBounds->left);
				r->top = ((r->top < subBounds->top) ? r->top : subBounds->top);
				r->right = ((r->right > subBounds->right) ? r->right : subBounds->right);
				r->bottom = ((r->bottom > subBounds->bottom) ? r->bottom : subBounds->bottom);
			}
		}
	}
	return r;
}

// *************************
// ** class CDXUnknownObject **
// *************************
//
// This is the object which is created when we run into a type which we don't know about
CDXObject* CDXUnknownObject::Clone() const
{
	return new CDXUnknownObject(*this);
}

std::string CDXUnknownObject::XMLObjectName() const
{
	return kCDXML_unknown;
}


// *************************
// ** class CDXAttributes **
// *************************

// CDXAttributes is a collection of CDXAttribute's.
//
// Although one might have multiple attributes of the same type in a given object,
// I think we're going to permit only a single attribute of a given type, so use
// a map, not a multimap.
// A CDXAttributes object "owns" its attributes, that is, they are deleted when
// the CDXAttributes is destructed.

// CDXAttributes::insert
//
// Insert an attribute into a CDXAttributes map.

pair<CDXAttributes::iterator,bool>
CDXAttributes::Insert(CDXAttribute *a)
{
	Remove (a);
	pair<CDXAttributes::iterator,bool>	result = insert( value_type(a->GetTag(), a) );
	return result;
}

// CDXAttributes::Remove
//
// Remove an attribute into a CDXAttributes map.

bool CDXAttributes::Remove (CDXAttribute* pAttr)
{
	ASSERT (pAttr != NULL);
	if (!pAttr)
		return false;
	iterator it = find (pAttr->GetTag());
	if (it == end())
		return false;
	delete GetAttribute(it);
	erase (it);
	return true;
}

// CDXAttributes::WriteTo
//
// Write all of the attributes in the collection

void CDXAttributes::WriteTo(CDXDataSink &sink) const
{
	for (const_iterator i = begin();  i != end();  ++i)
		GetAttribute(i)->WriteTo(sink);
}

// CDXAttributes::XMLWrite
//
// Write all of the attributes in the collection

void CDXAttributes::XMLWrite(XMLDataSink &sink) const
{
	for (const_iterator i = begin();  i != end();  ++i)
		GetAttribute(i)->XMLWrite(sink);
}

// CDXAttributes::~CDXAttributes
//
// Delete all of the attributes in the collection

CDXAttributes::~CDXAttributes()
{
	if (empty())
		return;
	for (const_iterator i = begin();  i != end();  ++i)
		delete GetAttribute(i);
	clear();
}

// invalid_cdx_error constructor
//
// This is the exception thrown when a problem is found with the format of the CDX stream

invalid_cdx_error::invalid_cdx_error(CDXObjectID objID_arg, CDXTag objTag_arg, CDXTag attribTag_arg)
		: m_objectID(objID_arg),
		  m_objectTag(objTag_arg),
		  m_attribTag(attribTag_arg)
{
	ostringstream s;
	s << "Invalid CDX format in object ID "
	  << m_objectID
	  << " type "
	  << m_objectTag;
	if (m_attribTag != 0)
		s << " attribute "
		  << m_attribTag;
	m_what = s.str();
}


map<CDXTag, string>											sCDXMLAttributeNames;
map<const char *, CDXTag, XMLCharPtrLess<const char *> >	sCDXMLAttributeIDs;

static void InitAttributeMapsIfNeeded(void)
{    
	if (!sCDXMLAttributeNames.empty())
		return;
    PROTECT_GLOBAL_AND_STATIC_DATA(s_lock);
	
	sCDXMLAttributeNames[kCDXProp_CreationUserName] = "CreationUserName";
	sCDXMLAttributeNames[kCDXProp_CreationDate] = "CreationDate";
	sCDXMLAttributeNames[kCDXProp_CreationProgram] = "CreationProgram";
	sCDXMLAttributeNames[kCDXProp_ModificationUserName] = "ModificationUserName";
	sCDXMLAttributeNames[kCDXProp_ModificationDate] = "ModificationDate";
	sCDXMLAttributeNames[kCDXProp_ModificationProgram] = "ModificationProgram";
	sCDXMLAttributeNames[kCDXProp_Name] = "Name";
	sCDXMLAttributeNames[kCDXProp_Comment] = "Comment";
	sCDXMLAttributeNames[kCDXProp_CartridgeData] = "CartridgeData";
	sCDXMLAttributeNames[kCDXProp_ZOrder] = "Z";
	sCDXMLAttributeNames[kCDXProp_RegistryNumber] = "RegistryNumber";
	sCDXMLAttributeNames[kCDXProp_RegistryAuthority] = "RegistryAuthority";
	sCDXMLAttributeNames[kCDXProp_RepresentsProperty] = "RepresentsProperty";
	sCDXMLAttributeNames[kCDXProp_FontTable] = "FontTable";
	sCDXMLAttributeNames[kCDXProp_2DPosition] = "p";
	sCDXMLAttributeNames[kCDXProp_3DPosition] = "xyz";
	sCDXMLAttributeNames[kCDXProp_2DExtent] = "extent";
	sCDXMLAttributeNames[kCDXProp_3DExtent] = "extent3D";
	sCDXMLAttributeNames[kCDXProp_BoundingBox] = "BoundingBox";
	sCDXMLAttributeNames[kCDXProp_FixInplaceExtent]  = "FixInPlaceExtent";
	sCDXMLAttributeNames[kCDXProp_FixInplaceGap]  = "FixInPlaceGap";
	sCDXMLAttributeNames[kCDXProp_Side]  = "Side";
	sCDXMLAttributeNames[kCDXProp_BoundsInParent] = "BoundsInParent";
	sCDXMLAttributeNames[kCDXProp_3DHead] = "Head3D";
	sCDXMLAttributeNames[kCDXProp_3DTail] = "Tail3D";
	sCDXMLAttributeNames[kCDXProp_TopLeft] = "TopLeft";
	sCDXMLAttributeNames[kCDXProp_TopRight] = "TopRight";
	sCDXMLAttributeNames[kCDXProp_BottomRight] = "BottomRight";
	sCDXMLAttributeNames[kCDXProp_BottomLeft] = "BottomLeft";
	sCDXMLAttributeNames[kCDXProp_3DCenter] = "Center3D";
	sCDXMLAttributeNames[kCDXProp_3DMajorAxisEnd] = "MajorAxisEnd3D";
	sCDXMLAttributeNames[kCDXProp_3DMinorAxisEnd] = "MinorAxisEnd3D";
	sCDXMLAttributeNames[kCDXProp_RotationAngle] = "RotationAngle";
	sCDXMLAttributeNames[kCDXProp_ColorTable] = "ColorTable";
	sCDXMLAttributeNames[kCDXProp_ForegroundColor] = "color";
	sCDXMLAttributeNames[kCDXProp_BackgroundColor] = "bgcolor";
	sCDXMLAttributeNames[kCDXProp_ForegroundAlpha] = "alpha";
	sCDXMLAttributeNames[kCDXProp_BackgroundAlpha] = "bgalpha";
    sCDXMLAttributeNames[kCDXProp_HighlightColor] = "highlightColor";
	sCDXMLAttributeNames[kCDXProp_Node_Type] = "NodeType";
	sCDXMLAttributeNames[kCDXProp_Node_LabelDisplay] = "LabelDisplay";
	sCDXMLAttributeNames[kCDXProp_Node_Element] = "Element";
	sCDXMLAttributeNames[kCDXProp_Node_NeedsClean] = "NeedsClean";
	sCDXMLAttributeNames[kCDXProp_Atom_ElementList] = "ElementList";
	sCDXMLAttributeNames[kCDXProp_Atom_GenericList] = "GenericList";
	sCDXMLAttributeNames[kCDXProp_Atom_Isotope] = "Isotope";
	sCDXMLAttributeNames[kCDXProp_Atom_Charge] = "Charge";
	sCDXMLAttributeNames[kCDXProp_Atom_Radical] = "Radical";
	sCDXMLAttributeNames[kCDXProp_Atom_Formula] = "NodeFormula";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictFreeSites] = "FreeSites";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictImplicitHydrogens] = "ImplicitHydrogens";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictRingBondCount] = "RingBondCount";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictUnsaturatedBonds] = "UnsaturatedBonds";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictRxnChange] = "RxnChange";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictRxnStereo] = "RxnStereo";
	sCDXMLAttributeNames[kCDXProp_Atom_Translation] = "Translation";
    sCDXMLAttributeNames[kCDXProp_Atom_IsotopicAbundance] = "IsotopicAbundance";
    sCDXMLAttributeNames[kCDXProp_Atom_ExternalConnectionType] = "ExternalConnectionType";
    sCDXMLAttributeNames[kCDXProp_Atom_ExternalConnectionNum] = "ExternalConnectionNum";
	sCDXMLAttributeNames[kCDXProp_Atom_AbnormalValence] = "AbnormalValence";
	sCDXMLAttributeNames[kCDXProp_Atom_NumHydrogens] = "NumHydrogens";
	sCDXMLAttributeNames[kCDXProp_Atom_HDot] = "HDot";
	sCDXMLAttributeNames[kCDXProp_Atom_HDash] = "HDash";
	sCDXMLAttributeNames[kCDXProp_Atom_Geometry] = "Geometry";
	sCDXMLAttributeNames[kCDXProp_Atom_BondOrdering] = "BondOrdering";
	sCDXMLAttributeNames[kCDXProp_Node_Attachments] = "Attachments";
    sCDXMLAttributeNames[kCDXProp_Node_HydrogenBondAttachmentAtoms] = "HydrogenBondAttachmentAtoms";
    sCDXMLAttributeNames[kCDXProp_Node_HydrogenBonds] = "HydrogenBonds";
	sCDXMLAttributeNames[kCDXProp_Atom_GenericNickname] = "GenericNickname";
	sCDXMLAttributeNames[kCDXProp_Atom_AtomNumber] = "AtomNumber";
	sCDXMLAttributeNames[kCDXProp_Atom_ResidueID] = "ResidueID";
	sCDXMLAttributeNames[kCDXProp_Atom_AltGroupID] = "AltGroupID";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictSubstituentsUpTo] = "SubstituentsUpTo";
	sCDXMLAttributeNames[kCDXProp_Atom_RestrictSubstituentsExactly] = "SubstituentsExactly";
	sCDXMLAttributeNames[kCDXProp_Atom_CIPStereochemistry] = "AS";
	sCDXMLAttributeNames[kCDXProp_Atom_LinkCountLow] = "LinkCountLow";
	sCDXMLAttributeNames[kCDXProp_Atom_LinkCountHigh] = "LinkCountHigh";
	sCDXMLAttributeNames[kCDXProp_Atom_EnhancedStereoType] = "EnhancedStereoType";
	sCDXMLAttributeNames[kCDXProp_Atom_EnhancedStereoGroupNum] = "EnhancedStereoGroupNum";
	sCDXMLAttributeNames[kCDXProp_Mole_Racemic] = "Racemic";
	sCDXMLAttributeNames[kCDXProp_Mole_Absolute] = "Absolute";
	sCDXMLAttributeNames[kCDXProp_Mole_Relative] = "Relative";
	sCDXMLAttributeNames[kCDXProp_Mole_Weight] = "Weight";
	sCDXMLAttributeNames[kCDXProp_Mole_Formula] = "Formula";
	sCDXMLAttributeNames[kCDXProp_Frag_ConnectionOrder] = "ConnectionOrder";
	sCDXMLAttributeNames[kCDXProp_Frag_SequenceType] = "SequenceType";
	sCDXMLAttributeNames[kCDXProp_Frag_IsFromGuidedStereo] = "StereoGuided";
	sCDXMLAttributeNames[kCDXProp_Bond_Order] = "Order";
	sCDXMLAttributeNames[kCDXProp_Bond_Connectivity] = "Connectivity";
    sCDXMLAttributeNames[kCDXProp_Bond_Connectivity_Routed] = "Routed";
	sCDXMLAttributeNames[kCDXProp_Bond_Display] = "Display";
	sCDXMLAttributeNames[kCDXProp_Bond_Display2] = "Display2";
	sCDXMLAttributeNames[kCDXProp_Bond_DoublePosition] = "DoublePosition";
	sCDXMLAttributeNames[kCDXProp_Bond_Begin] = "B";
	sCDXMLAttributeNames[kCDXProp_Bond_End] = "E";
	sCDXMLAttributeNames[kCDXProp_Bond_RestrictTopology] = "Topology";
	sCDXMLAttributeNames[kCDXProp_Bond_RestrictRxnParticipation] = "RxnParticipation";
	sCDXMLAttributeNames[kCDXProp_Bond_BeginAttach] = "BeginAttach";
	sCDXMLAttributeNames[kCDXProp_Bond_EndAttach] = "EndAttach";
	sCDXMLAttributeNames[kCDXProp_Bond_CIPStereochemistry] = "BS";
	sCDXMLAttributeNames[kCDXProp_Bond_BondOrdering] = "BondCircularOrdering";
	sCDXMLAttributeNames[kCDXProp_Bond_CrossingBonds] = "CrossingBonds";
    sCDXMLAttributeNames[kCDXProp_Bond_BeginExternalNum] = "BeginExternalNum";
    sCDXMLAttributeNames[kCDXProp_Bond_EndExternalNum] = "EndExternalNum";
	sCDXMLAttributeNames[kCDXProp_Text] = "Text";
	sCDXMLAttributeNames[kCDXProp_Justification] = "Justification";
	sCDXMLAttributeNames[kCDXProp_LabelJustification] = "LabelJustification";
	sCDXMLAttributeNames[kCDXProp_CaptionJustification] = "CaptionJustification";
	sCDXMLAttributeNames[kCDXProp_AminoAcidTermini] = "AminoAcidTermini";
	sCDXMLAttributeNames[kCDXProp_ShowSequenceTermini] = "ShowSequenceTermini";
	sCDXMLAttributeNames[kCDXProp_ShowSequenceBonds] = "ShowSequenceBonds";
    sCDXMLAttributeNames[kCDXProp_ShowSequenceUnlinkedBranches] = "ShowSequenceUnlinkedBranches";
    sCDXMLAttributeNames[kCDXProp_MonomerRenderingStyle] = "MonomerRenderingStyle";
	sCDXMLAttributeNames[kCDXProp_ResidueWrapCount] = "ResidueWrapCount";
	sCDXMLAttributeNames[kCDXProp_ResidueBlockCount] = "ResidueBlockCount";
	sCDXMLAttributeNames[kCDXProp_LineHeight] = "LineHeight";
	sCDXMLAttributeNames[kCDXProp_LabelLineHeight] = "LabelLineHeight";
	sCDXMLAttributeNames[kCDXProp_CaptionLineHeight] = "CaptionLineHeight";
	sCDXMLAttributeNames[kCDXProp_WordWrapWidth] = "WordWrapWidth";
	sCDXMLAttributeNames[kCDXProp_LineStarts] = "LineStarts";
	sCDXMLAttributeNames[kCDXProp_LabelAlignment] = "LabelAlignment";
	sCDXMLAttributeNames[kCDXProp_InterpretChemically] = "InterpretChemically";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowTerminalCarbonLabels] = "ShowTerminalCarbonLabels";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowNonTerminalCarbonLabels] = "ShowNonTerminalCarbonLabels";
	sCDXMLAttributeNames[kCDXProp_Atom_HideImplicitHydrogens] = "HideImplicitHydrogens";
	sCDXMLAttributeNames[kCDXProp_Positioning] = "PositioningType";
	sCDXMLAttributeNames[kCDXProp_PositioningAngle] = "PositioningAngle";
	sCDXMLAttributeNames[kCDXProp_PositioningOffset] = "PositioningOffset";
	sCDXMLAttributeNames[kCDXProp_MacPrintInfo] = "MacPrintInfo";
	sCDXMLAttributeNames[kCDXProp_WinPrintInfo] = "WinPrintInfo";
	sCDXMLAttributeNames[kCDXProp_PrintMargins] = "PrintMargins";
	sCDXMLAttributeNames[kCDXProp_ChainAngle] = "ChainAngle";
	sCDXMLAttributeNames[kCDXProp_BondSpacing] = "BondSpacing";
	sCDXMLAttributeNames[kCDXProp_BondSpacingAbs] = "BondSpacingAbs";
	sCDXMLAttributeNames[kCDXProp_BondSpacingType] = "BondSpacingType";
	sCDXMLAttributeNames[kCDXProp_BondLength] = "BondLength";
	sCDXMLAttributeNames[kCDXProp_BoldWidth] = "BoldWidth";
	sCDXMLAttributeNames[kCDXProp_LineWidth] = "LineWidth";
	sCDXMLAttributeNames[kCDXProp_MarginWidth] = "MarginWidth";
	sCDXMLAttributeNames[kCDXProp_HashSpacing] = "HashSpacing";
	sCDXMLAttributeNames[kCDXProp_LabelStyleFont] = "LabelFont";
	sCDXMLAttributeNames[kCDXProp_CaptionStyleFont] = "CaptionFont";
	sCDXMLAttributeNames[kCDXProp_LabelStyleSize] = "LabelSize";
	sCDXMLAttributeNames[kCDXProp_CaptionStyleSize] = "CaptionSize";
	sCDXMLAttributeNames[kCDXProp_LabelStyleFace] = "LabelFace";
	sCDXMLAttributeNames[kCDXProp_CaptionStyleFace] = "CaptionFace";
	sCDXMLAttributeNames[kCDXProp_LabelStyleColor] = "LabelColor";
	sCDXMLAttributeNames[kCDXProp_CaptionStyleColor] = "CaptionColor";
	sCDXMLAttributeNames[kCDXProp_CaptionJustification] = "CaptionJustification";
	sCDXMLAttributeNames[kCDXProp_FractionalWidths] = "FractionalWidths";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowQuery] = "ShowAtomQuery";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowStereo] = "ShowAtomStereo";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowAtomNumber] = "ShowAtomNumber";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowResidueID] = "ShowResidueID";
	sCDXMLAttributeNames[kCDXProp_Atom_ShowEnhancedStereo] = "ShowAtomEnhancedStereo";
    sCDXMLAttributeNames[kCDXProp_Atom_ShowAtomID] = "ShowAtomID";
    sCDXMLAttributeNames[kCDXProp_Atom_AtomID] = "AtomID";
    sCDXMLAttributeNames[kCDXProp_Bond_ShowQuery] = "ShowBondQuery";
	sCDXMLAttributeNames[kCDXProp_Bond_ShowRxn] = "ShowBondRxn";
	sCDXMLAttributeNames[kCDXProp_Bond_ShowStereo] = "ShowBondStereo";
	sCDXMLAttributeNames[kCDXProp_Magnification] = "Magnification";
	sCDXMLAttributeNames[kCDXProp_WidthPages] = "WidthPages";
	sCDXMLAttributeNames[kCDXProp_HeightPages] = "HeightPages";
	sCDXMLAttributeNames[kCDXProp_DrawingSpaceType] = "DrawingSpace";
	sCDXMLAttributeNames[kCDXProp_Width] = "Width";
	sCDXMLAttributeNames[kCDXProp_Height] = "Height";
	sCDXMLAttributeNames[kCDXProp_PageOverlap] = "PageOverlap";
	sCDXMLAttributeNames[kCDXProp_Header] = "Header";
	sCDXMLAttributeNames[kCDXProp_HeaderPosition] = "HeaderPosition";
	sCDXMLAttributeNames[kCDXProp_Footer] = "Footer";
	sCDXMLAttributeNames[kCDXProp_FooterPosition] = "FooterPosition";
	sCDXMLAttributeNames[kCDXProp_PrintTrimMarks] = "PrintTrimMarks";
	sCDXMLAttributeNames[kCDXProp_Window_IsZoomed] = "WindowIsZoomed";
	sCDXMLAttributeNames[kCDXProp_Window_Position] = "WindowPosition";
	sCDXMLAttributeNames[kCDXProp_Window_Size] = "WindowSize";
	sCDXMLAttributeNames[kCDXProp_Graphic_Type] = "GraphicType";
	sCDXMLAttributeNames[kCDXProp_Line_Type] = "LineType";
	sCDXMLAttributeNames[kCDXProp_Arrow_Type] = "ArrowType";
	sCDXMLAttributeNames[kCDXProp_Rectangle_Type] = "RectangleType";
	sCDXMLAttributeNames[kCDXProp_Oval_Type] = "OvalType";
	sCDXMLAttributeNames[kCDXProp_Orbital_Type] = "OrbitalType";
	sCDXMLAttributeNames[kCDXProp_Bracket_Type] = "BracketType";
	sCDXMLAttributeNames[kCDXProp_Symbol_Type] = "SymbolType";
	sCDXMLAttributeNames[kCDXProp_Curve_Type] = "CurveType";
	sCDXMLAttributeNames[kCDXProp_Frame_Type] = "FrameType";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_Size] = "HeadSize";
	sCDXMLAttributeNames[kCDXProp_Arc_AngularSize] = "AngularSize";
	sCDXMLAttributeNames[kCDXProp_CornerRadius] = "CornerRadius";
	sCDXMLAttributeNames[kCDXProp_Bracket_LipSize] = "LipSize";
	sCDXMLAttributeNames[kCDXProp_Bracket_Usage] = "BracketUsage";
	sCDXMLAttributeNames[kCDXProp_Polymer_RepeatPattern] = "PolymerRepeatPattern";
	sCDXMLAttributeNames[kCDXProp_Polymer_FlipType] = "PolymerFlipType";
	sCDXMLAttributeNames[kCDXProp_Curve_Points] = "CurvePoints";
	sCDXMLAttributeNames[kCDXProp_Curve_Points3D] = "CurvePoints3D";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_Type] = "ArrowheadType";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_Head] = "ArrowheadHead";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_Tail] = "ArrowheadTail";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_CenterSize] = "ArrowheadCenterSize";
	sCDXMLAttributeNames[kCDXProp_Arrowhead_Width] = "ArrowheadWidth";
	sCDXMLAttributeNames[kCDXProp_ShadowSize] = "ShadowSize";
	sCDXMLAttributeNames[kCDXProp_Arrow_ShaftSpacing] = "ArrowShaftSpacing";
	sCDXMLAttributeNames[kCDXProp_Arrow_EquilibriumRatio] = "ArrowEquilibriumRatio";
	sCDXMLAttributeNames[kCDXProp_Arrow_SourceID] = "ArrowSource";
	sCDXMLAttributeNames[kCDXProp_Arrow_TargetID] = "ArrowTarget";
	sCDXMLAttributeNames[kCDXProp_Fill_Type] = "FillType";
	sCDXMLAttributeNames[kCDXProp_Curve_Spacing] = "CurveSpacing";
	sCDXMLAttributeNames[kCDXProp_Closed] = "Closed";
	sCDXMLAttributeNames[kCDXProp_Arrow_NoGo] = "NoGo";
	sCDXMLAttributeNames[kCDXProp_Arrow_Dipole] = "Dipole";
	sCDXMLAttributeNames[kCDXProp_Picture_Edition] = "Edition";
	sCDXMLAttributeNames[kCDXProp_Picture_EditionAlias] = "EditionAlias";
	sCDXMLAttributeNames[kCDXProp_MacPICT] = "MacPICT";
	sCDXMLAttributeNames[kCDXProp_WindowsMetafile] = "WindowsMetafile";
	sCDXMLAttributeNames[kCDXProp_EnhancedMetafile] = "EnhancedMetafile";
	sCDXMLAttributeNames[kCDXProp_OLEObject] = "OLEObject";
	sCDXMLAttributeNames[kCDXProp_Compressed_MacPICT] = "CompressedMacPICT";
	sCDXMLAttributeNames[kCDXProp_Compressed_WindowsMetafile] = "CompressedWindowsMetafile";
	sCDXMLAttributeNames[kCDXProp_Compressed_EnhancedMetafile] = "CompressedEnhancedMetafile";
	sCDXMLAttributeNames[kCDXProp_Compressed_OLEObject] = "CompressedOLEObject";
	sCDXMLAttributeNames[kCDXProp_Uncompressed_MacPICT_Size] = "UncompressedMacPICTSize";
	sCDXMLAttributeNames[kCDXProp_Uncompressed_WindowsMetafile_Size] = "UncompressedWindowsMetafileSize";
	sCDXMLAttributeNames[kCDXProp_Uncompressed_EnhancedMetafile_Size] = "UncompressedEnhancedMetafileSize";
	sCDXMLAttributeNames[kCDXProp_Uncompressed_OLEObject_Size] = "UncompressedOLEObjectSize";
	sCDXMLAttributeNames[kCDXProp_GIF] = "GIF";
	sCDXMLAttributeNames[kCDXProp_TIFF] = "TIFF";
	sCDXMLAttributeNames[kCDXProp_PNG] = "PNG";
	sCDXMLAttributeNames[kCDXProp_PDF] = "PDF";
	sCDXMLAttributeNames[kCDXProp_JPEG] = "JPEG";
	sCDXMLAttributeNames[kCDXProp_BMP] = "BMP";
	sCDXMLAttributeNames[kCDXProp_Spectrum_XSpacing] = "XSpacing";
	sCDXMLAttributeNames[kCDXProp_Spectrum_XLow] = "XLow";
	sCDXMLAttributeNames[kCDXProp_Spectrum_XType] = "XType";
	sCDXMLAttributeNames[kCDXProp_Spectrum_YType] = "YType";
	sCDXMLAttributeNames[kCDXProp_Spectrum_XAxisLabel] = "XAxisLabel";
	sCDXMLAttributeNames[kCDXProp_Spectrum_YAxisLabel] = "YAxisLabel";
	sCDXMLAttributeNames[kCDXProp_Spectrum_Class] = "Class";
	sCDXMLAttributeNames[kCDXProp_Spectrum_YLow] = "YLow";
	sCDXMLAttributeNames[kCDXProp_Spectrum_YScale] = "YScale";
	sCDXMLAttributeNames[kCDXProp_TLC_OriginFraction] = "OriginFraction";
	sCDXMLAttributeNames[kCDXProp_TLC_SolventFrontFraction] = "SolventFrontFraction";
	sCDXMLAttributeNames[kCDXProp_TLC_ShowOrigin] = "ShowOrigin";
	sCDXMLAttributeNames[kCDXProp_TLC_ShowSolventFront] = "ShowSolventFront";
	sCDXMLAttributeNames[kCDXProp_ShowBorders] = "ShowBorders";
	sCDXMLAttributeNames[kCDXProp_TLC_ShowSideTicks] = "ShowSideTicks";
	sCDXMLAttributeNames[kCDXProp_TLC_Rf] = "Rf";
	sCDXMLAttributeNames[kCDXProp_TLC_Tail] = "Tail";
	sCDXMLAttributeNames[kCDXProp_TLC_ShowRf] = "ShowRf";
	sCDXMLAttributeNames[kCDXProp_GEP_ShowScale] = "ShowScale";
	sCDXMLAttributeNames[kCDXProp_GEP_LaneLabelsAngle] = "LabelsAngle";
	sCDXMLAttributeNames[kCDXProp_GEP_ScaleUnit] = "UnitID";
	sCDXMLAttributeNames[kCDXProp_GEP_StartRange] = "StartRange";
	sCDXMLAttributeNames[kCDXProp_GEP_EndRange] = "EndRange";
	sCDXMLAttributeNames[kCDXProp_GEP_Value] = "BandValue";
	sCDXMLAttributeNames[kCDXProp_GEP_ShowValue] = "ShowValue";
	sCDXMLAttributeNames[kCDXProp_GEP_AxisWidth] = "AxisWidth";
	sCDXMLAttributeNames[kCDXProp_NamedAlternativeGroup_TextFrame] = "TextFrame";
	sCDXMLAttributeNames[kCDXProp_NamedAlternativeGroup_GroupFrame] = "GroupFrame";
	sCDXMLAttributeNames[kCDXProp_NamedAlternativeGroup_Valence] = "Valence";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Reactants] = "ReactionStepReactants";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Products] = "ReactionStepProducts";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Plusses] = "ReactionStepPlusses";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Arrows] = "ReactionStepArrows";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_ObjectsAboveArrow] = "ReactionStepObjectsAboveArrow";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_ObjectsBelowArrow] = "ReactionStepObjectsBelowArrow";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Atom_Map] = "ReactionStepAtomMap";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Atom_Map_Manual] = "ReactionStepAtomMapManual";
	sCDXMLAttributeNames[kCDXProp_ReactionStep_Atom_Map_Auto] = "ReactionStepAtomMapAuto";
    sCDXMLAttributeNames[kCDXProp_RxnAutonumber_Conditions] = "RxnAutonumberConditions";
    sCDXMLAttributeNames[kCDXProp_RxnAutonumber_Style] = "RxnAutonumberStyle";
    sCDXMLAttributeNames[kCDXProp_RxnAutonumber_Start] = "RxnAutonumberStart";
    sCDXMLAttributeNames[kCDXProp_RxnAutonumber_Format] = "RxnAutonumberFormat";
	sCDXMLAttributeNames[kCDXProp_Template_PaneHeight] = "PaneHeight";
	sCDXMLAttributeNames[kCDXProp_Template_NumRows] = "NumRows";
	sCDXMLAttributeNames[kCDXProp_Template_NumColumns] = "NumColumns";
	sCDXMLAttributeNames[kCDXProp_Group_Integral] = "Integral";
	sCDXMLAttributeNames[kCDXProp_ObjectTag_Type] = "TagType";
	sCDXMLAttributeNames[kCDXProp_ObjectTag_Tracking] = "Tracking";
	sCDXMLAttributeNames[kCDXProp_ObjectTag_Persistent] = "Persistent";
	sCDXMLAttributeNames[kCDXProp_ObjectTag_Value] = "Value";
	sCDXMLAttributeNames[kCDXProp_SplitterPositions] = "SplitterPositions";
	sCDXMLAttributeNames[kCDXProp_PageDefinition] = "PageDefinition";
	sCDXMLAttributeNames[kCDXProp_IgnoreWarnings] = "IgnoreWarnings";
	sCDXMLAttributeNames[kCDXProp_ChemicalWarning] = "Warning";
	sCDXMLAttributeNames[kCDXProp_Visible] = "Visible";
	sCDXMLAttributeNames[kCDXProp_Transparent] = "Transparent";
	sCDXMLAttributeNames[kCDXProp_SupersededBy] = "SupersededBy";
	sCDXMLAttributeNames[kCDXProp_Sequence_Identifier] = "SequenceIdentifier";
	sCDXMLAttributeNames[kCDXProp_CrossReference_Container] = "CrossReferenceContainer";
	sCDXMLAttributeNames[kCDXProp_CrossReference_Document] = "CrossReferenceDocument";
	sCDXMLAttributeNames[kCDXProp_CrossReference_Identifier] = "CrossReferenceIdentifier";
	sCDXMLAttributeNames[kCDXProp_CrossReference_Sequence] = "CrossReferenceSequence";
	sCDXMLAttributeNames[kCDXProp_BracketedObjects] = "BracketedObjectIDs";
	sCDXMLAttributeNames[kCDXProp_Bracket_RepeatCount] = "RepeatCount";
	sCDXMLAttributeNames[kCDXProp_Bracket_ComponentOrder] = "ComponentOrder";
	sCDXMLAttributeNames[kCDXProp_Bracket_SRULabel] = "SRULabel";
	sCDXMLAttributeNames[kCDXProp_Bracket_GraphicID] = "GraphicID";
	sCDXMLAttributeNames[kCDXProp_Bracket_BondID] = "BondID";
	sCDXMLAttributeNames[kCDXProp_Bracket_InnerAtomID] = "InnerAtomID";
	sCDXMLAttributeNames[kCDXProp_GeometricFeature] = "GeometricFeature";
	sCDXMLAttributeNames[kCDXProp_RelationValue] = "RelationValue";
	sCDXMLAttributeNames[kCDXProp_BasisObjects] = "BasisObjects";
	sCDXMLAttributeNames[kCDXProp_ConstraintType] = "ConstraintType";
	sCDXMLAttributeNames[kCDXProp_ConstraintMin] = "ConstraintMin";
	sCDXMLAttributeNames[kCDXProp_ConstraintMax] = "ConstraintMax";
	sCDXMLAttributeNames[kCDXProp_IgnoreUnconnectedAtoms] = "IgnoreUnconnectedAtoms";
	sCDXMLAttributeNames[kCDXProp_DihedralIsChiral] = "DihedralIsChiral";
	sCDXMLAttributeNames[kCDXProp_PointIsDirected] = "PointIsDirected";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyType] = "ChemicalPropertyType";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyIsActive] = "ChemicalPropertyIsActive";
    sCDXMLAttributeNames[kCDXProp_ChemicalPropertyIsChemicallySignificant] = "ChemicallySignificant";
    sCDXMLAttributeNames[kCDXProp_ChemicalPropertyExternalBonds] = "ExternalBonds";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyDisplayID] = "ChemicalPropertyDisplayID";
	sCDXMLAttributeNames[kCDXProp_SG_DataValue] = "SGDataValue";
	sCDXMLAttributeNames[kCDXProp_SG_DataType] = "SGDataType";
	sCDXMLAttributeNames[kCDXProp_SG_PropertyType] = "SGPropertyType";
	sCDXMLAttributeNames[kCDXProp_SG_ComponentIsReactant] = "ComponentIsReactant";
	sCDXMLAttributeNames[kCDXProp_SG_ComponentIsHeader] = "ComponentIsHeader";
	sCDXMLAttributeNames[kCDXProp_IsHidden] = "IsHidden";
	sCDXMLAttributeNames[kCDXProp_IsReadOnly] = "IsReadOnly";
	sCDXMLAttributeNames[kCDXProp_IsEdited] = "IsEdited";
	sCDXMLAttributeNames[kCDXProp_SG_ComponentReferenceID] = "ComponentReferenceID";
 	sCDXMLAttributeNames[kCDXProp_BioShape_Type] = "BioShapeType";
	sCDXMLAttributeNames[kCDXProp_1SubstrateEnzyme_ReceptorSize]  = "EnzymeReceptorSize";
	sCDXMLAttributeNames[kCDXProp_Receptor_NeckWidth]    = "NeckWidth";
	sCDXMLAttributeNames[kCDXProp_HelixProtein_CylinderWidth]    = "CylinderWidth";
	sCDXMLAttributeNames[kCDXProp_HelixProtein_CylinderHeight]   = "CylinderHeight";
	sCDXMLAttributeNames[kCDXProp_HelixProtein_CylinderDistance] = "CylinderDistance";
	sCDXMLAttributeNames[kCDXProp_HelixProtein_PipeWidth]     = "PipeWidth";
	sCDXMLAttributeNames[kCDXProp_HelixProtein_Extra]     = "HelixProteinExtra";
	sCDXMLAttributeNames[kCDXProp_Membrane_ElementSize]   = "MembraneElementSize";
	sCDXMLAttributeNames[kCDXProp_Membrane_StartAngle]    = "MembraneStartAngle";
	sCDXMLAttributeNames[kCDXProp_Membrane_EndAngle]      = "MembraneEndAngle";
	sCDXMLAttributeNames[kCDXProp_DNA_WaveLength] = "DNAWaveLength";
	sCDXMLAttributeNames[kCDXProp_DNA_WaveWidth]  = "DNAWaveWidth";
	sCDXMLAttributeNames[kCDXProp_DNA_Offset]     = "DNAWaveOffset";
	sCDXMLAttributeNames[kCDXProp_DNA_WaveHeight] = "DNAWaveHeight";
	sCDXMLAttributeNames[kCDXProp_Gprotein_UpperHeight] = "GproteinUpperHeight";
	sCDXMLAttributeNames[kCDXProp_FadePercent] = "FadePercent";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_NumberBasePairs] = "NumberBasePairs";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_MarkerStart] = "MarkerStart";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_MarkerOffset] = "MarkerOffset";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_MarkerAngle] = "MarkerAngle";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_RegionStart] = "RegionStart";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_RegionEnd] = "RegionEnd";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_RegionOffset] = "RegionOffset";
	sCDXMLAttributeNames[kCDXProp_PlasmidMap_RingRadius] = "RingRadius";
    sCDXMLAttributeNames[kCDXProp_ChemicalPropertyName] = "ChemPropName";
    sCDXMLAttributeNames[kCDXProp_ChemicalPropertyID] = "ChemPropID";
    sCDXMLAttributeNames[kCDXProp_ChemicalPropertyFragmentLabel] = "ChemPropFragmentLabel";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyFormula] = "ChemPropFormula";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyExactMass] = "ChemPropExactMass";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyMolWeight] = "ChemPropMolWt";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyMOverZ] = "ChemPropMOverZ";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyAnalysis] = "ChemPropAnalysis";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyBoilingPoint] = "ChemPropBoilingPt";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyMeltingPoint] = "ChemPropMeltingPt";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyCriticalTemp] = "ChemPropCritTemp";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyCriticalPressure] = "ChemPropCritPres";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyCriticalVolume] = "ChemPropCritVol";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyGibbsEnergy] = "ChemPropGibbs";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyLogP] = "ChemPropLogP";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyMR] = "ChemPropMR";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyHenrysLaw] = "ChemPropHenry";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyHeatOfForm] = "ChemPropEForm";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertytPSA] = "ChemProptPSA";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyCLogP] = "ChemPropCLogP";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyCMR] = "ChemPropCMR";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyLogS] = "ChemPropLogS";
	sCDXMLAttributeNames[kCDXProp_ChemicalPropertyPKa] = "ChemPropPKa";
	sCDXMLAttributeNames[kCDXProp_RLogic_Group] = "RLogicGroup";
	sCDXMLAttributeNames[kCDXProp_RLogic_Occurrence] = "RLogicOccurrence";
	sCDXMLAttributeNames[kCDXProp_RLogic_RestH] = "RLogicRestH";
	sCDXMLAttributeNames[kCDXProp_RLogic_IfThenGroup] = "RLogicIfThenGroup";
	sCDXMLAttributeNames[kCDXProp_Annotation_Keyword] = "Keyword";
	sCDXMLAttributeNames[kCDXProp_Annotation_Content] = "Content";
    sCDXMLAttributeNames[kCDXProp_Property_Rule] = "Rule";
    sCDXMLAttributeNames[kCDXProp_Property_DataType] = "Type";
    sCDXMLAttributeNames[kCDXProp_StructurePerspective] = "ShowPerspective";



	// Now make the reverse map
	for (map<CDXTag, string>::const_iterator j = sCDXMLAttributeNames.begin();  j != sCDXMLAttributeNames.end();  ++j)
	{
		// Make sure no two IDs share the same name
		ASSERT (sCDXMLAttributeIDs.find(j->second.c_str()) == sCDXMLAttributeIDs.end());
		sCDXMLAttributeIDs[j->second.c_str()] = j->first;
	}
}

CDXTag CDXMLAttributeID(const char *s)
{
	InitAttributeMapsIfNeeded();

	map<const char *, CDXTag, XMLCharPtrLess<const char *> >::const_iterator i = sCDXMLAttributeIDs.find(s);
	return (i == sCDXMLAttributeIDs.end()) ? 0 : i->second;
}

const string &CDXMLAttributeName(CDXTag t)
{	
	InitAttributeMapsIfNeeded();	
	map<CDXTag, string>::const_iterator i = sCDXMLAttributeNames.find(t);
	//if (i == sXMLAttributeNames.end())
	//	OutputDebugString("Could not find name for tag in XMLAttributeName\n");
	//ASSERT (i != sXMLAttributeNames.end());
	if (i == sCDXMLAttributeNames.end())
	{
		// This stinks, but it allows us to return a const reference in the by-far-most-common
		// case of a recognized attribute, while still doing something acceptable for the very rare
		// unknown attributes.
		static map<CDXTag, string> unknownTagsMap;
		map<CDXTag, string>::const_iterator it = unknownTagsMap.find(t);
		if (it == unknownTagsMap.end())
		{
			ostringstream os;
			os << "attrib" << hex << setw(4) << setfill('0') << t;
			unknownTagsMap[t] = os.str();
			it = unknownTagsMap.find(t);
		}
		return it->second;
	}
	else
		return i->second;
}

void CDXMLPutAttribute(XMLDataSink &sink, CDXTag tag, bool v)
{
    if (v)
        sink.os << " " << CDXMLAttributeName(tag) << "=\"yes\"" << GetTextEOL();
    else
        sink.os << " " << CDXMLAttributeName(tag) << "=\"no\"" << GetTextEOL();
}

void CDXMLPutAttribute(XMLDataSink &sink, vector<CDXChemProp> props)
{
	for (vector<CDXChemProp>::iterator it = props.begin(); it != props.end(); ++it)
		CDXMLPutAttribute(sink, (*it).datumID, (*it).propString.str());
}

/**
 *  Writes a space separated basis object ID list to CDXML sink
 *
 *  @param sink XML data sink to write to
 *  @param objectIDList A object ID list
 */
void CDXMLPutAttributeForObjectIDList(XMLDataSink &sink, CDXTag tag, const vector<CDXObjectID>& objectIDList)
{
    if (!objectIDList.empty())
    {
        sink.os << string(" ") << CDXMLAttributeName(kCDXProp_BasisObjects) << "=\"";
        for (auto objectIterator = objectIDList.cbegin(); objectIterator != objectIDList.cend(); objectIterator++)
        {
            if (objectIterator != objectIDList.cbegin())
            {
                sink.os << " ";
            }

            sink.os << *objectIterator;
        }

        sink.os << "\"" << GetTextEOL();
    }
}

void XMLPut(XMLDataSink &sink, long  v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, double  v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, int   v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, short v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, char  v)
{
	sink.os << int(v);
}

void XMLPut(XMLDataSink &sink, unsigned long  v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, unsigned int   v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, unsigned short v)
{
	sink.os << v;
}

void XMLPut(XMLDataSink &sink, unsigned char  v)
{
	sink.os << (unsigned int)v;
}

void XMLPut(XMLDataSink &sink, bool  v)
{
	sink.os << (v ? "yes" : "no");
}

void XMLPut(XMLDataSink &sink, const string &v)
{
    // Strip characters that are not valid XML
    string s = RemoveIllegalXMLChars(v);    
    s = MakeStringSafe(s);
	string::size_type	curPos = 0,
						nextPos;
	while ((nextPos = s.find_first_of ("&<>\'\"\n\r", curPos)) != string::npos)
	{
		if (nextPos > curPos)
			sink.os.write(s.data() + curPos, (int)(nextPos - curPos));
		switch (s[nextPos])
		{
		case '&':
			sink.os << "&amp;";
			break;
		case '<':
			sink.os << "&lt;";
			break;
		case '>':
			sink.os << "&gt;";
			break;
		case '\'':
			sink.os << "&apos;";
			break;
		case '\"':
			sink.os << "&quot;";
			break;
		case '\n':
		case '\r':
			// TY: this should be "&#013;" and "&#010;"
			// But, we can't change this because ChemDraw has written lots of data already
			sink.os << GetTextEOL();
			break;
		default:
			ASSERT(false); // should never get here as the list of chars should all be in the switch
			break;
		}
		curPos = nextPos + 1;
	}
	if (s.size() > curPos)
		sink.os.write(s.data() + curPos, (int)(s.size() - curPos));
}

void XMLTranslateEscapedString(const string& in, string& out)
{
	// expatapp already takes care of them except "&#013;" and "&#010;"
	const string ampString		= "&amp;";
	const string ltString		= "&lt;";
	const string gtString		= "&gt;";
	const string aposString		= "&apos;";
	const string quotString		= "&quot;";
	const string crString		= "&#013;";
	const string nlString		= "&#010;";

	ostringstream os;

	string::size_type	curPos = 0,
						nextPos;
	while ((nextPos = in.find_first_of("&", curPos)) != string::npos)
	{
		if (nextPos > curPos)
		{
			os.write(in.data() + curPos, (int)(nextPos - curPos));
		}

		if (strcmp(&in[nextPos], ampString.c_str()) == 0)
		{
			os << '&';
			curPos = nextPos + ampString.length();
		}
		else if (strcmp(&in[nextPos], ltString.c_str()) == 0)
		{
			os << '<';
			curPos = nextPos + ltString.length();
		}
		else if (strcmp(&in[nextPos], gtString.c_str()) == 0)
		{
			os << '>';
			curPos = nextPos + gtString.length();
		}
		else if (strcmp(&in[nextPos], aposString.c_str()) == 0)
		{
			os << '\'';
			curPos = nextPos + aposString.length();
		}
		else if (strcmp(&in[nextPos], quotString.c_str()) == 0)
		{
			os << '\"';
			curPos = nextPos + quotString.length();
		}
		else if (strcmp(&in[nextPos], crString.c_str()) == 0)
		{
			os << '\r';
			curPos = nextPos + crString.length();
		}
		else if (strcmp(&in[nextPos], nlString.c_str()) == 0)
		{
			os << '\n';
			curPos = nextPos + nlString.length();
		}
		else
		{
			ASSERT(false);
			break;
		}
	}

	if (in.size() > curPos)
		os.write(in.data() + curPos, (int)(in.size() - curPos));

//	os.freeze(); 
	if (os.str().empty())
		out = string();
	else
		//out = string(os.str(), os.str().length()); AMT: CSBR 64449
		out = string(os.str());
}

string CDXValueToString(double v, int precision)
{
	// Convert a value to fixed point, leaving off unnecessary trailing zeroes.
	ostringstream os;
	os << fixed << setprecision(precision) << v;
	string result = os.str();
	if (precision > 0) // otherwise there won't be a decimal point
	{
		// Erase the trailing zeroes (and the decimal point, too, if all trailing digits are zeroes)
		string::size_type p = result.find_last_not_of('0');
		if (p != string::npos && result[p] == '.')
		{
			--p;
			if (p != string::npos)
				result.erase(p+1, string::npos);
		}
	}
	return result;
}

void XMLPut(XMLDataSink &sink, const CDXPoint2D &v)
{
	sink.os << CDXCoordinateToString(v.x) << " " << CDXCoordinateToString(v.y);
}

void XMLPut(XMLDataSink &sink, const Point32 &v)
{
	sink.os << v.x << " " << v.y;
}

void XMLPut(XMLDataSink &sink, const CDXPoint3D &v)
{
	sink.os << CDXCoordinateToString(v.x) << " " << CDXCoordinateToString(v.y) << " " << CDXCoordinateToString(v.z);
}

void XMLPut(XMLDataSink &sink, const CDXFloat64Point3D &v)
{
	sink.os << v.x << " " << v.y << " " << v.z;
}

void XMLPut(XMLDataSink &sink, const CDXFloat64Matrix3D &v)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			sink.os << " " << v.m[i][j];
}

void XMLPut(XMLDataSink &sink, const CDXRectangle &v)
{
	sink.os << CDXCoordinateToString(v.left) << " " << CDXCoordinateToString(v.top) << " " << CDXCoordinateToString(v.right) << " " << CDXCoordinateToString(v.bottom);
}

void XMLPut(XMLDataSink &sink, const CDXFontStyle &v)
{
	sink.os << v;
}

void CDXPut(CDXDataSink& sink, const string& s)
{
    const string stringToStore = MakeStringSafe(s);
    sink.Put(UINT16(stringToStore.size() + 2));
    sink.Put(UINT16(0));	// the styles
    sink.Put((INT8*)stringToStore.data(), stringToStore.size());
}
