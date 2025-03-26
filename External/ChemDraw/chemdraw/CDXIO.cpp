// CommonCS/LibCommon/Src/CDXIO.cpp
// Contains: Program-independent class library for managing CDX objects
//				This file contains the implementation of the low-level file i/o methods
//				for the CDXDataSource and CDXDataSink objects.
//				It also contains methods for handling low-level types such as CDXString.
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

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

#include "cs_univDefs.h"
#include "CDXIO.h"
#include "cs_swapBytes.h"
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "cs_auto_buffer.h"

using namespace std;

#define kCharCR		'\015'
#define kCharLF		'\012'

#define kEOLCR		"\015"
#define kEOLLF		"\012"
#define kEOLCRLF	kEOLCR kEOLLF

#if TARGET_OS_MAC || TARGET_WEB
static string textEOL(kEOLCR);
#elif TARGET_OS_WIN32
static string textEOL(kEOLCRLF);
#else
static string textEOL(kEOLLF);
#endif // TARGET_OS_MAC

 void SetTextEOL(EOLIndex index)
{
	switch (index) {
		case kEOLIndexCR:
			SetTextEOL(kEOLCR);
			break;
		case kEOLIndexLF:
			SetTextEOL(kEOLLF);
			break;
		case kEOLIndexCRLF:
			SetTextEOL(kEOLCRLF);
			break;
	}
}

void SetTextEOL(const string &s)
{
	textEOL = s;
}

string GetTextEOL()
{
	return textEOL;
}

void NormalizeEOL(std::string *s)
{
	string eol = GetTextEOL();
	bool wantCR = textEOL == kEOLCR;
	bool wantLF = textEOL == kEOLLF;
	bool wantCRLF = textEOL == kEOLCRLF;

	string::size_type i = 0;
	for (;;)
	{
		i = s->find_first_of(kEOLCRLF, i);
		if (i == string::npos)
			break;
		if ((*s)[i] == kCharLF) // A linefeed is always bare
		{
			if (wantCR)
				(*s)[i] = kCharCR;
			else if (wantCRLF)
			{
				s->insert(i, 1, kCharCR);
				++i;
			}
		}
		else if (i == s->size() - 1 || (*s)[i+1] != kCharLF) // then it's a bare CR
		{
			if (wantLF)
				(*s)[i] = kCharLF;
			else if (wantCRLF)
			{
				s->insert(i+1, 1, kCharLF);
				++i;
			}
		}
		else // Must be a CRLF
		{
			if (wantLF)
				s->erase(i, 1);
			else if (wantCR)
				s->erase(i+1, 1);
			else if (wantCRLF)
				++i;
		}
		++i;
	}
}

CDXString::CDXString(const string_t &s)
{
    SetText(s);
}

CDXString::CDXString(const char *s) : m_string(s)
{
    SetText(s);
}

void CDXString::SetText(const string_t &str)
{
    m_string = str;

    // Truncate at the first null character if there is one
    m_string.resize(strlen(m_string.c_str()));

    // Remove any styles that overflow the new text
    ClearOverflowStyles();
}

CDXDataSource & operator>>(CDXDataSource &src_arg, CDXStyle &style_arg)
{
	style_arg.startChar = src_arg.GetUINT16();
	style_arg.family    = src_arg.GetUINT16();
	style_arg.face      = src_arg.GetUINT16();
	style_arg.size      = src_arg.GetUINT16();
	style_arg.color     = src_arg.GetUINT16();
	return src_arg;
}

// ** function CDXStyle::Read **
// Read a CDX style from a data source
CDXStyle::CDXStyle(CDXDataSource &src_arg) :
	startChar(0),
	family(0xFFFF),
	face(0),
	size(0),
	color(3),
	alpha(0xFFFF),
	has_alpha(false)
{
	src_arg >> *this;
}

// ** function CDXStyle::Write **
// Write a CDX style to a data sink
CDXDataSink& operator<< (CDXDataSink& sink_arg, const CDXStyle &style_arg)
{
	return sink_arg << UINT16(style_arg.startChar)
					<< UINT16(style_arg.family)
					<< UINT16(style_arg.face)
					<< UINT16(style_arg.size)
					<< UINT16(style_arg.color);
}

// ** function CDXString::Read **
// Read a CDX string (including styles and text) from a data source
void CDXString::Read(CDXDataSource &src_arg, size_t size_arg)
{
	// A CDXString consists of
	//    UINT16 number of styles
	//    That many styles (each one is CDXStyle::CDXDataSize bytes long)
	//    zero or more bytes of string data
    string readString;

	// Get the number of styles, making sure there's enough data
	if (size_arg < kStyleCountSize)
	{
		// Rather than throwing in the absence of style, assume that the src_arg is a raw c_string
		readString.assign(src_arg.GetString(size_arg));
	}
	else
	{
		UINT16 numStyles = src_arg.GetUINT16();

		// Make sure there's enough data, then read the styles into the m_styles vector
		const size_t stylesSize = size_t(numStyles * CDXStyle::CDXDataSize);

		if (size_arg < kStyleCountSize + stylesSize)
		{
			// Rather than throwing in the absence of style, assume that the src_arg is a raw c_string
            // But first we copy over the 2 bytes we just read
            readString.assign((char*)&numStyles, kStyleCountSize);

            // Then read any remaining data
            const size_t remaining = size_arg - kStyleCountSize;
            if (remaining > 0)
            {
                readString += src_arg.GetString(remaining);
            }
		}
		else
		{
			m_styles.reserve(numStyles);
			generate_n(back_inserter(m_styles), numStyles, CDXDataSourceArray<CDXStyle>(src_arg));

			// Now read the actual string bytes
			const size_t stringSize = size_arg - kStyleCountSize - stylesSize;
			readString.assign(src_arg.GetString(stringSize));
		}
	}

    SetText(readString);
}

void CDXString::operator += (char c)
{
    m_string += c;
}

void CDXString::operator += (const CDXString &str)
{
	const unsigned int	oldLen = (unsigned int)m_string.size();
	m_string += str.m_string;
	// Append the new string's styles, offsetting the character position.
	for (unsigned int i = 0;  i < str.m_styles.size();  i++)
	{
		m_styles.push_back (str.m_styles [i]);
		m_styles.back().startChar += oldLen;
	}
}

void CDXString::operator += (const string &str)
{
	// add some text to the last style run
	m_string += str;
}

CDXDataSink& operator<< (CDXDataSink& sink_arg, const vector<CDXStyle> &styles_arg)
{
	vector<CDXStyle>::const_iterator	styleIter = styles_arg.begin();
	vector<CDXStyle>::const_iterator	styleTerminus = styles_arg.end();
#pragma warning (disable:4996)	// "'std::copy': Function call with parameters that may be unsafe."
	std::copy(styleIter, styleTerminus, CDXDataSinkIterator<CDXStyle>(sink_arg));
	return sink_arg;
}

void CDXString::PutTo(CDXDataSink& sink_arg) const
{
	sink_arg << UINT16(nstyles()) << styles();
	sink_arg.Put((INT8 *)c_str(), length());
}

/**
 *  @return size of data that would be written to a CDX stream
 */
size_t CDXString::CDXDataSize() const
{
    return kStyleCountSize + CDXStyle::CDXDataSize * nstyles() + length();
}

/**
 *  Add a style to our array
 */
void CDXString::AddStyle(const CDXStyle &style)
{
    m_styles.push_back(style);
}

/**
 *  Clear all our styles
 */
void CDXString::ClearStyles()
{
    m_styles.clear();
}

/**
 *  Clear any styles that overflow the end of our string
 */
void CDXString::ClearOverflowStyles()
{
    while (!m_styles.empty() && (m_styles.back().startChar != UINT16(-1)) && (m_styles.back().startChar >= m_string.size()))
    {
        m_styles.pop_back();
    }
}

// ** function CDXDataSource::GetUINT16 **
// Read an unsigned 16-bit little-endian integer from the data source
UINT16 CDXDataSource::GetUINT16()
{
	unsigned char p[2];
	GetBytes(p, 2);
	if (m_bigEndian)
		return (p[0] << 8) | p[1];
	else
		return (p[1] << 8) | p[0];
}

// ** function CDXDataSource::GetUINT32 **
// Read an unsigned 32-bit little-endian integer from the data source
UINT32 CDXDataSource::GetUINT32()
{
	unsigned char p[4];
	GetBytes(p, 4);
	if (m_bigEndian)
		return (p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
	else
		return (p[3] << 24) | (p[2] << 16) | (p[1] << 8) | p[0];
}

FLOAT64 CDXDataSource::GetFLOAT64()
{
	FLOAT64 f64;
	GetBytes((char *)&f64, sizeof(f64));
	return SwapBytes(f64);
}

UINT32 CDXDataSource::GetUINT(size_t size_arg)
{
	UINT32 returnValue = 0;
	switch(size_arg)
	{
	case 1:
		returnValue = GetUINT8();
		break;
	case 2:
		returnValue = GetUINT16();
		break;
	case 4:
		returnValue = GetUINT32();
		break;
	default:
	{
//		throw invalid_argument("Invalid integer size");
		while (size_arg > 4)
		{
			returnValue = GetUINT(4);
			size_arg -= 4;
		}
		if (size_arg > 0)
			return GetUINT(size_arg);
	}
	} // switch(size_arg)
	return returnValue;
}

INT32 CDXDataSource::GetINT(size_t size_arg)
{
	INT32 returnValue;
	switch(size_arg)
	{
	case 1:
		returnValue = GetINT8();
		break;
	case 2:
		returnValue = GetINT16();
		break;
	case 4:
		returnValue = GetINT32();
		break;
	default:
		throw std::invalid_argument("Invalid integer size");
	} // switch(size_arg)
	return returnValue;
}

std::string CDXDataSource::GetString(size_t n) 
{
	if (n == 0)
		return std::string();	//quite boundschecker
	auto_buffer<char> buf((unsigned int)n); 
	GetBytes(buf.ptr(), n); 
	return std::string(buf.ptr(), n); 
}

void CDXistream::GetBytes( char *p, size_t n )
{
	is.read(p, (int)n);
	if (!is)
		throw std::runtime_error("Unexpected end of input or other read error");
}

void CDXipointer::GetBytes( char *p, size_t n )
{
	if (m_pos + n > m_end)
		throw std::runtime_error("Unexpected end of input or other read error");
	memcpy(p, m_pos, n);
	m_pos += n;
}

void CDXistream::SkipBytes( size_t n )
{
	is.ignore((int)n);
	if (!is)
		throw std::runtime_error("Unexpected end of input or other read error");
}

void CDXipointer::SkipBytes( size_t n )
{
	if (m_pos + n > m_end)
		throw std::runtime_error("Unexpected end of input or other read error");
	m_pos += n;
}

// ** function CDXDataSink::PutUINT16 **
// Write an unsigned 16-bit little-endian integer to the data sink
void CDXDataSink::Put( UINT16 datum )
{
	unsigned char p[2];
	p[1] = (datum >>  8) & 0xFF;
	p[0] = datum & 0xFF;
	Put(p, 2);
}

// ** function CDXDataSink::Put( UINT32 datum ) **
// Write an unsigned 32-bit little-endian integer to the data sink
void CDXDataSink::Put( UINT32 datum )
{
	unsigned char p[4];
	p[3] = (unsigned char) (datum >> 24);
	p[2] = (unsigned char) (datum >> 16);
	p[1] = (unsigned char) (datum >>  8);
	p[0] = (unsigned char) (datum      );
	Put(p, 4);
}

// ** function CDXDataSink::Put( FLOAT64 datum ) **
// Write a 64-bit little-endian float to the data sink
void CDXDataSink::Put( FLOAT64 datum )
{
	char p[8];
	SwapBytes(&datum, p);
	Put((INT8 *)p, 8);
}

void CDXDataSink::PutAttributeForObjectIDList(CDXTag tag, const vector<CDXObjectID>& objectIDList)
{
    if (objectIDList.empty())
    {
        return;
    }

    const size_t numBytes = objectIDList.size() * sizeof(CDXObjectID);
    PutTag(kCDXProp_BasisObjects);
    Put(UINT16(numBytes));
    for (auto basisObject : objectIDList)
    {
        Put((UINT32)basisObject);
    }
}

// ** function CDXostream::Put( INT8 *data, size_t n ) **
// Write a sequence of bytes to the data sink
void CDXostream::Put( const INT8 *data, size_t n )
{
	os.write((const char *)data, (int)n);
}

ostream& operator<<(ostream& os, const CDXFontStyle &s)
{
	return os << "{" << s.family << "," << s.face << "," << s.size << "," << s.color << "}";
}

ostream& operator<<(ostream& os, const CDXStyle &s)
{
	os << "{";
	os << s.startChar;
	os << "," << s.family;
	os << "," << s.face;
	os << "," << s.size;
	os << "," << s.color;
	if (s.has_alpha)
		os << "," << s.alpha;
	os << "}";
	return os;
}

ostream& operator<<(ostream& os, const CDXString &s)
{
	return os << s.str() << s.styles();
}


ostream& operator<<(ostream& os, CDXAtomGeometry s)
{
	switch(s)
	{
	case kCDXAtomGeometry_Unknown:				os << "Unknown"; break;
	case kCDXAtomGeometry_1Ligand:				os << "1Ligand"; break;
	case kCDXAtomGeometry_Linear:				os << "Linear"; break;
	case kCDXAtomGeometry_Bent:					os << "Bent"; break;
	case kCDXAtomGeometry_TrigonalPlanar:		os << "TrigonalPlanar"; break;
	case kCDXAtomGeometry_TrigonalPyramidal:	os << "TrigonalPyramidal"; break;
	case kCDXAtomGeometry_SquarePlanar:			os << "SquarePlanar"; break;
	case kCDXAtomGeometry_Tetrahedral:			os << "Tetrahedral"; break;
	case kCDXAtomGeometry_TrigonalBipyramidal:	os << "TrigonalBipyramidal"; break;
	case kCDXAtomGeometry_SquarePyramidal:		os << "SquarePyramidal"; break;
	case kCDXAtomGeometry_5Ligand:				os << "5Ligand"; break;
	case kCDXAtomGeometry_Octahedral:			os << "Octahedral"; break;
	case kCDXAtomGeometry_6Ligand:				os << "6Ligand"; break;
	case kCDXAtomGeometry_7Ligand:				os << "7Ligand"; break;
	case kCDXAtomGeometry_8Ligand:				os << "8Ligand"; break;
	case kCDXAtomGeometry_9Ligand:				os << "9Ligand"; break;
	case kCDXAtomGeometry_10Ligand:				os << "10Ligand"; break;
	default:									os << "Geometry#" << int(s); break;
	}
	return os;
}

// *********************************************************
// ** ostream& operator<<(ostream& os, const CDXPoint2D &r) **
// *********************************************************
//
// Sends a CDXPoint2D to an output stream
ostream& operator<<(ostream& os, const CDXPoint2D &p)
{
	os << "(" << CDXCoordinateToString(p.x) << ", " << CDXCoordinateToString(p.y) << ")";
	return os;
}

// *********************************************************
// ** ostream& operator<<(ostream& os, const CDXPoint3D &r) **
// *********************************************************
//
// Sends a CDXPoint2D to an output stream
ostream& operator<<(ostream& os, const CDXPoint3D &p)
{
	os << "(" << CDXCoordinateToString(p.x) << ", " << CDXCoordinateToString(p.y) << ", " << CDXCoordinateToString(p.z) << ")";
	return os;
}

// *********************************************************
// ** ostream& operator<<(ostream& os, const CDXRectangle &r) **
// *********************************************************
//
// Sends a CDXRectangle to an output stream
ostream& operator<<(ostream& os, const CDXRectangle &r)
{
	os << "(ltrb = "
		<< CDXCoordinateToString(r.left)
		<< ", "
		<< CDXCoordinateToString(r.top)
		<< ", "
		<< CDXCoordinateToString(r.right)
		<< ", "
		<< CDXCoordinateToString(r.bottom)
		<< ")";
	return os;
}

// *********************************************************
// ** ostream& operator<<(ostream& os, CDXDatumID r) **
// *********************************************************
//
// Sends a CDXDatumID to an output stream
 ostream& operator<<(ostream& os, CDXDatumID tag)
{
	switch(tag)
	{
	case kCDXObj_Document:				os << "Document";				break;
	case kCDXObj_Page:					os << "Page";					break;
	case kCDXObj_Group:					os << "Group";					break;
	case kCDXObj_Fragment:				os << "Fragment";				break;
	case kCDXObj_Node:					os << "Node";					break;
	case kCDXObj_Bond:					os << "Bond";					break;
	case kCDXObj_Text:					os << "Text";					break;
	case kCDXObj_Graphic:				os << "Graphic";				break;
	case kCDXObj_Curve:					os << "Curve";					break;
	case kCDXObj_EmbeddedObject:		os << "EmbeddedObject";			break;
	case kCDXObj_NamedAlternativeGroup: os << "NamedAlternativeGroup";	break;
	case kCDXObj_TemplateGrid:			os << "TemplateGrid";			break;
	case kCDXObj_RegistryNumber:		os << "RegistryNumber";			break;
	case kCDXObj_ReactionScheme:		os << "ReactionScheme";			break;
	case kCDXObj_ReactionStep:			os << "ReactionStep";			break;
	case kCDXObj_ObjectDefinition:		os << "ObjectDefinition";		break;
	case kCDXObj_Spectrum:				os << "Spectrum";				break;
	case kCDXObj_ObjectTag:				os << "ObjectTag";				break;
	case kCDXObj_Sequence:				os << "Sequence";				break;
	case kCDXObj_CrossReference:		os << "CrossReference";			break;
	case kCDXObj_Splitter:				os << "Splitter";				break;
	case kCDXObj_Table:					os << "Table";					break;
	case kCDXObj_BracketedGroup:		os << "BracketedGroup";			break;
	case kCDXObj_BracketAttachment:		os << "BracketAttachment";		break;
	case kCDXObj_CrossingBond:			os << "CrossingBond";			break;
	case kCDXObj_Border:				os << "Border";					break;
	case kCDXObj_Geometry:				os << "Geometry";				break;
	case kCDXObj_Constraint:			os << "Constraint";				break;
	case kCDXObj_TLCPlate:				os << "TLCPlate";				break;
	case kCDXObj_TLCLane:				os << "TLCLane";				break;
	case kCDXObj_TLCSpot:				os << "TLCSpot";				break;
	case kCDXObj_ChemicalProperty:		os << "ChemicalProperty";		break;
	case kCDXObj_Arrow:					os << "Arrow";					break;
	case kCDXObj_StoichiometryGrid:		os << "StoichiometryGrid";		break;
	case kCDXObj_SGComponent:			os << "Compoenent";				break;
	case kCDXObj_SGDatum:				os << "StoichiometryDatum";		break;
	case kCDXObj_BioShape:				os << "BioShape";				break;
	case kCDXObj_PlasmidMap:			os << "PlasmidMap";				break;
	case kCDXObj_PlasmidMarker:			os << "PlasmidMarker";			break;
	case kCDXObj_PlasmidRegion:			os << "PlasmidRegion";			break;
	case kCDXObj_RLogic:				os << "RLogic";					break;
	case kCDXObj_RLogicItem:			os << "RLogicItem";				break;
	case kCDXObj_Annotation:			os << "Annotation";				break;
	case kCDXObj_GEPPlate:				os << "GEPPlate";				break;
	case kCDXObj_GEPLane:				os << "GEPLane";				break;
	case kCDXObj_GEPBand:				os << "GEPBand";				break;
	case kCDXObj_Marker:				os << "Marker";					break;
    case kCDXObj_DocumentProperties:	os << "DocumentProperties";		break;
    case kCDXObj_Property:			    os << "Property";				break;
    case kCDXObj_ColoredMolecularArea:  os << "ColoredMolecularArea";           break;
	default:
		{
			std::ios_base::fmtflags flags = os.flags();
			//std::streamsize osp = os.precision();
			std::streamsize osw = os.width();
			char ofc = os.fill();

			os << ((tag & kCDXTag_Object) ? "Obj#" : "Att#")
				<< right
				<< hex << setw(4) << setfill('0') << UINT16(tag) << dec;

			os.fill(ofc);
			os.width(osw);
			//os.precision(osp);
			os.flags(flags);
		}
		break;
	}
	return os;
}

CDXPoint2D StringToCDXPoint2D(const std::string &s)
{
	const char *c = s.c_str();
	char *rest = NULL;
	double x = strtod(c, &rest);
	double y = strtod(rest, &rest);
	return CDXPoint2D(CDXCoordinatefromPoints(x), CDXCoordinatefromPoints(y));
}

CDXPoint3D StringToCDXPoint3D(const std::string &s)
{
	const char *c = s.c_str();
	char *rest = NULL;
	double x = strtod(c, &rest);
	double y = strtod(rest, &rest);
	double z = strtod(rest, &rest);
	return CDXPoint3D(CDXCoordinatefromPoints(x), CDXCoordinatefromPoints(y), CDXCoordinatefromPoints(z));
}

CDXRectangle StringToCDXRectangle(const std::string &s)
{
	const char *c = s.c_str();
	char *rest = NULL;
	double left = strtod(c, &rest);
	double top = strtod(rest, &rest);
	double right = strtod(rest, &rest);
	double bottom = strtod(rest, &rest);
	return CDXRectangle(CDXCoordinatefromPoints(top), CDXCoordinatefromPoints(left),
						CDXCoordinatefromPoints(bottom), CDXCoordinatefromPoints(right));
}

// *********************************************************
// ***** CDXWriteHeader   Write a standard CDX header ******
// *********************************************************
//
void CDXWriteHeader (CDXDataSink &sink)
{
	sink.Put ((INT8 *)kCDX_HeaderString, 8);
	sink.Put ((UINT32)0x01020304L);	// Reserved magic #

	// TOC pointer, then six bytes reserved.
	// Used to be 12 bytes, but the object-oriented code writes 6 bytes worth of tags for the top-level document object
	// and we want older versions to ignored those.  (CSBR-34111:Files not compatible with ChemDraw 5.0 or earlier)
	sink.Put ((UINT32)0);	// TOC pointer
	sink.Put ((UINT32)0);	// reserved
	sink.Put ((UINT16)0);	// reserved
}


#if CDX_DUMP
//---------------------------------------------------------------------
void CDXFileDumper::On()
{
	m_dump = 1;
	m_file = fopen("C:\\DBG.TXT", "w");
}
//---------------------------------------------------------------------
void CDXFileDumper::Off()
{
	if (m_dump) {
		m_dump = 0;
		fclose(m_file);
		m_file = 0;
	}
}
//---------------------------------------------------------------------
void CDXFileDumper::Suspend(bool on)
{
	m_dump = !on;
}
//---------------------------------------------------------------------
void CDXFileDumper::DumpEndObject()
{
	strcpy(m_str, "---------------------------------");
	Out();
}
//---------------------------------------------------------------------
static char *objname[] = {
	"Document", "Page", "Group", "Fragment", "Node", "Bond", "Text",						// 0x8006
	"Graphic", "Curve", "EmbeddedObject", "NamedAlternativeGroup", "TemplateGrid",				// 0x800b
	"RegistryNumber", "ObjectDefinition", "Spectrum", "ObjectTagDictionaryEntry", "ObjectTag"
};

//---------------------------------------------------------------------
void CDXFileDumper::DumpObj(int tag)
{
	sprintf(m_str, "Obj 0x0%x [%s]", tag, objname[tag - 0x8000]);
	Out();
}
//---------------------------------------------------------------------
void CDXFileDumper::DumpAttrib(int tag, int len)
{
	sprintf(m_str, "Att 0x0%x len=%d ", tag, len);
	char *str = NULL;
	switch (tag) {
	case kCDXProp_2DPosition:			str = "2DPosition";	 			break;	
	case kCDXProp_BoundingBox:			str = "BoundingBox"; 			break;	
	case kCDXProp_Node_LabelDisplay:	str = "LabelDisplay	";			break;	
	case kCDXProp_Node_Element:			str = "Element";  				break;	
	case kCDXProp_Atom_ElementList:		str = "ElementList"; 			break;	
	case kCDXProp_Formula:				str = "Formula";  				break;	
	case kCDXProp_Atom_Isotope:			str = "Isotope";				break;	
	case kCDXProp_Atom_Charge:			str = "Charge";					break;	
	case kCDXProp_Atom_Radical:			str = "Radical";  				break;	
	case kCDXProp_Bond_Order:			str = "Order";					break;			
    case kCDXProp_Bond_Connectivity:	str = "Connectivity";			break;
    case kCDXProp_Bond_BeginExternalNum:    str = "BeginExternalNum";	break;
    case kCDXProp_Bond_EndExternalNum:  str = "EndExternalNum";         break;
    case kCDXProp_Bond_Connectivity:	str = "Connectivity";			break;
	case kCDXProp_Bond_Display:			str = "Display";				break;
	case kCDXProp_Bond_Display2:		str = "Display2";				break;			
	case kCDXProp_Bond_DoublePosition:	str = "DoublePosition";			break;	
	case kCDXProp_Bond_Begin:			str = "Begin";					break;			
	case kCDXProp_Bond_End:				str = "End";					break;				
	case kCDXProp_Text:					str = "Text";					break;					
	case kCDXProp_Justification:		str = "Justification";			break;		
	case kCDXProp_LabelJustification:	str = "LabelJustification";		break;		
	case kCDXProp_CaptionJustification:	str = "CaptionJustification";	break;		
	case kCDXProp_LineHeight:			str = "LineHeight";				break;		
	case kCDXProp_LabelLineHeight:		str = "LabelLineHeight";		break;		
	case kCDXProp_CaptionLineHeight:	str = "CaptionLineHeight";		break;		
	case kCDXProp_WordWrapWidth:		str = "WordWrapWidth";			break;		
	case kCDXProp_LineStarts:			str = "LineStarts";				break;		
	case kCDXProp_LabelAlignment:		str = "LabelAlignment";			break;	
	case kCDXProp_InterpretChemically:	str = "InterpretChemically";	break;		
	case kCDXProp_Graphic_Type:			str = "Graphic_Type";			break;				
	case kCDXProp_Line_Type:			str = "Line_Type";				break;				
	case kCDXProp_Arrow_Type:			str = "Arrow_Type";				break;			
	case kCDXProp_Rectangle_Type:		str = "Rectangle_Type";			break;		
	case kCDXProp_Oval_Type:			str = "Oval_Type";				break;			
	case kCDXProp_Orbital_Type:			str = "Orbital_Type";			break;		
	case kCDXProp_Bracket_Type:			str = "Bracket_Type";			break;			
	case kCDXProp_Symbol_Type:			str = "Symbol_Type";			break;			
	case kCDXProp_Curve_Type:			str = "Curve_Type";				break;			
	case kCDXProp_Arrow_HeadSize:		str = "Arrow_HeadSize";			break;			
	case kCDXProp_Arc_AngularSize:		str = "Arc_AngularSize";		break;	
	case kCDXProp_Bracket_LipSize:		str = "Bracket_LipSize";		break;			
	case kCDXProp_Curve_Points:			str = "Curve_Points";			break;			
	}
	if (str) {
		strcat(m_str, "[");	strcat(m_str, str); strcat(m_str, "]");
	}
	Out();
}
//---------------------------------------------------------------------
void CDXFileDumper::Dump(void *buf, int n)
{
	if (n == 1)
		sprintf(m_str, "byte %d", *(INT8*)buf);
	else if (n == 2)
		sprintf(m_str, "word %d [0x0%x]", *(UINT16*)buf, *(UINT16*)buf);
	else if (n == 4)
		sprintf(m_str, "long %ld [0x0%lx]", *(UINT32*)buf, *(UINT32*)buf);
	else {
		sprintf(m_str, "arr %d [", n);
		char val[10];
		int i, m; INT8 *u;
		for (i = 0, m = min(n, 10), u = (INT8*)buf; i < m; ++i, ++u) {
			sprintf(val, "%x ", *u);
			strcat(m_str, val);
		}
		if (m == 10) strcat(m_str, "...");
		strcat(m_str, "]");
	}
	Out();
}
//---------------------------------------------------------------------
void CDXFileDumper::Out()
{
	if (!m_dump) return;
	strcat(m_str, GetTextEOL().c_str());
	fputs(m_str, m_file);
}
//---------------------------------------------------------------------
#endif // CDX_DUMP
