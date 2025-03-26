// CommonCS/LibCommon/Hdr/CDXIO.h
// Contains: Program independent class library for managing CDX objects
// Copyright (c) 1999-2008, CambridgeSoft Corp., All Rights Reserved

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
// Don't warn about casts that truncate constant values
#pragma warning (disable:4310)
// turn off warnings about unreferenced formal parameter
#pragma warning (disable:4100)
#endif

#include "CoreChemistryAPI.h"
//#include "cs_univDefs.h"
#include "CDXConstants.h"

#ifdef __hpux
#include <iostream>
#endif

#include <stdlib.h> // for atoi
#include <string.h> // for strcmp
#include <string>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <set>
#include <functional>
#include <map>

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
using std::streambuf;
using std::endl;
using std::hex;
using std::vector;
using std::ostream;
using std::min;
using std::max;
using std::string;

#if !TARGET_OS_WIN32 // 'cause it's not in std:: on windows yet
using std::size_t;
#endif

#ifdef __linux
#include <stdexcept>
#endif
#ifndef TYPENAME
#if (TARGET_OS_WIN32  &&  defined(_MSC_VER)  &&  _MSC_VER > 1300) || TARGET_OS_MAC  ||  TARGET_OS_LINUX || __hpux
#define TYPENAME	typename
#else
#define TYPENAME
#endif
#endif

typedef enum
{
	kEOLIndexCR = 0,
	kEOLIndexLF,
	kEOLIndexCRLF
} EOLIndex;

CORE_CHEMISTRY_API extern void SetTextEOL(EOLIndex index);
CORE_CHEMISTRY_API extern void SetTextEOL(const std::string &s);
CORE_CHEMISTRY_API extern std::string GetTextEOL();
CORE_CHEMISTRY_API extern void NormalizeEOL(std::string *);

// Classes defined in this file:

class CDXDataSource;			// Abstract class describing interface for input
class CDXistream;				// Implementation of CDXDataSource for input from C++ streams
class CDXDataSink;				// Abstract class describing interface for output
class CDXostream;				// Implementation of CDXDataSink for output to C++ streams

typedef INT32 CDXCoordinate;	// simple classes for geometry interfaces
inline double CDXCoordinateToPoints(CDXCoordinate x) { return x / double(0x10000); }
inline CDXCoordinate CDXCoordinatefromPoints(double x) { return (CDXCoordinate)(x * double(0x10000)); }
std::string CDXValueToString(double v, int precision);
inline std::string CDXCoordinateToString(CDXCoordinate v) { return CDXValueToString(CDXCoordinateToPoints(v), 2); }

class CDXObject;
class CDXPoint2D;
class CDXPoint3D;
class CDXRectangle;
class CDXString;
class CDXStyle;

// C3XColorRGBA is defined in ChemDraw-Chem3D/Src/C3X/Hdr/C3XObjects.h and used in 
// ChemDraw-Chem3D/Src/C3X/Src/C3XDataUtils.cpp

#if !defined(CORE_CHEMISTRY_API_DLL_BUILD)
    class C3XColorRGBA;  // TODO: C3XColorRGBA leakage from Chem3D into Core Chemistry. This should be fixed in Chem3D. 
#endif

inline double CdxAngleToDegrees (INT32 cdxAngle)	{ return cdxAngle / 65536.; }
inline INT32  DegreesToCdxAngle (double degrees)	{ while (degrees < 0) degrees += 360; return (INT32)(degrees  * 65536. + 0.5); }

// ********************
// * class CDXFontStyle
// ********************
//
// This contains the style for a single run of identically-styled text in a CDXString
class CDXFontStyle {
public:
	UINT16			family;		// index into the font table
	UINT16			face;		// bit-encoded
	UINT16			size;
	UINT16			color;		// index into the palette

	CDXFontStyle (UINT16 fam = -1, UINT16 fac = 0, UINT16 siz = 0, UINT16 col = -1)
		:	family (fam), face (fac), size (siz), color (col)	{}
};

// ********************
// ** class CDXStyle **
// ********************
//
// This contains the style for a single run of identically-styled text in a CDXString
class CDXStyle {
public:
	enum { CDXDataSize = 10};	// #bytes stored in a CDX file for a CDXStyle: 5 UINT16's - note alpha value deliberately omitted for
								// compatibility with releases prior to version 12. Modifying this number will cause old versions to fail.

	UINT16			startChar;	// index into the string
	UINT16			family;		// index into the font table
	UINT16			face;		// bit-encoded
	UINT16			size;		// in units of 1/20 point.
	UINT16			color;		// index into the palette
	UINT16			alpha;		// alpha value
	bool			has_alpha;	// is alpha value defined


	CDXStyle(UINT16 sch = 0, UINT16 fam = 0xFFFF, UINT16 fac = 0, UINT16 siz = 0, UINT16 col = 3, UINT16 alf = 0xFFFF, bool hasalf = false)
		:    startChar(sch), family (fam), face (fac), size (siz), color (col), alpha (alf), has_alpha (hasalf)
		{ }
	CDXStyle(CDXDataSource &);
	bool	operator == (const CDXStyle& op2) const	{ return 0 == memcmp (this, &op2, sizeof (CDXStyle)); }
};
const int	kNumUnitsPerPoint = 20;	// There are 20 CDX units per (font) point unit.

CDXDataSink&   operator<< (CDXDataSink&  sink_arg, const CDXStyle &style_arg);
CDXDataSource& operator>> (CDXDataSource& src_arg,  CDXStyle &style_arg);
CDXDataSink&   operator<< (CDXDataSink&  sink_arg, const std::vector<CDXStyle> &styles_arg);



// *********************
// ** class CDXString **
// *********************
//
// This contains a string together with style information
class CORE_CHEMISTRY_API CDXString
{
public:

	typedef std::string				string_t;
	typedef std::vector<CDXStyle>	styles_t;

private:

	string_t m_string;
	styles_t m_styles;
    static const size_t kStyleCountSize = sizeof(UINT16);

public:

	CDXString() {}
    explicit CDXString(const string_t &s);
    explicit CDXString(const char *s);
	bool	operator == (const CDXString& op2) const	{ return m_string == op2.m_string && m_styles == op2.m_styles; }    
	string_t::size_type			length()			const { return m_string.length(); }
	const char *				c_str()				const { return m_string.c_str(); }
	bool						empty()				const { return m_string.empty(); }
	styles_t::size_type			nstyles()			const { return m_styles.size(); }
	const styles_t &			styles()			const { return m_styles; }
	styles_t &					styles()				  { return m_styles; }
	const CDXStyle &			style(size_t i)		const { return m_styles[i]; }
	CDXStyle &					style(size_t i)			  { return m_styles[i]; }
	const string_t &			str()				const { return m_string; }
	size_t						stylestart(size_t i)	const { return m_styles[i].startChar; }
	size_t						styleend(size_t i)		const { return i+1 < nstyles() ? style((int)i+1).startChar : length(); }
	size_t						stylelength(size_t i)	const { return styleend(i) - stylestart(i); }
		// styleNumAtChar returns the index of the style in effect at character number "charNum",
		// or -1 if no style is defined for the character.
	int							styleNumAtChar (int charNum)	const { for (int j = (int)(nstyles() - 1);  j >= 0;  --j)
																	if (style (j).startChar <= charNum)
																		return j;
																return -1;
															  }
	const string_t				str(size_t i)			const { return m_string.substr(stylestart(i), stylelength(i)); }
    size_t						CDXDataSize() const;
	void						Clear()						{ m_string.clear(); m_styles.clear(); }
	void						Read(CDXDataSource &src_arg, size_t size_arg);
    void						SetText(const string_t &str);
    void						AddStyle(const CDXStyle &style);
    void						ClearStyles();
    void						ClearOverflowStyles();
	void						operator += (const CDXString &str);
	void						operator += (const string_t &str);
    void						operator += (char c);
	void						PutTo(CDXDataSink& sink_arg) const;
};



class Point16
{
public:
	INT16 x;
	INT16 y;
	Point16(INT16 x_arg = 0,
			INT16 y_arg = 0)
			  : x(x_arg),
			    y(y_arg)
			  { }
};

class Point32
{
public:
	INT32 x;
	INT32 y;
	Point32(INT32 x_arg = 0,
			INT32 y_arg = 0)
			  : x(x_arg),
			    y(y_arg)
			  { }
};

// ***********************
// ** class CDXPoint2D  **
// ***********************
//
// CDXPoint2D is a simple 2 dimensional point, used in the interface to the geometry
// of CDXGraphicObjects

class CDXPoint2D
{
public:
	CDXCoordinate x;
	CDXCoordinate y;
	CDXPoint2D(CDXCoordinate x_arg = 0,
			   CDXCoordinate y_arg = 0)
			  : x(x_arg),
			    y(y_arg)
			  { }

	void Offset(CDXCoordinate x_arg, CDXCoordinate y_arg) { x += x_arg; y += y_arg; }
};

CDXPoint2D StringToCDXPoint2D(const std::string &s);

// ***********************
// ** class CDXPoint3D  **
// ***********************
//
// CDXPoint3D is a simple 3 dimensional point, used in the interface to the geometry
// of CDXGraphicObjects

class CDXPoint3D
{
public:
	CDXCoordinate x;
	CDXCoordinate y;
	CDXCoordinate z;
	CDXPoint3D(CDXCoordinate x_arg = 0,
			   CDXCoordinate y_arg = 0,
			   CDXCoordinate z_arg = 0)
			  : x(x_arg),
			    y(y_arg),
				z(z_arg)
			  { }

	void Offset(CDXCoordinate x_arg, CDXCoordinate y_arg, CDXCoordinate z_arg) { x += x_arg; y += y_arg; z += z_arg; }
};

CDXPoint3D StringToCDXPoint3D(const std::string &s);

// ******************************
// ** class CDXFloat64Point3D  **
// ******************************
//
// CDXFloat64Point3D is a simple 3 dimensional point defined using FLOAT64

class CDXFloat64Point3D
{
public:
	FLOAT64 x;
	FLOAT64 y;
	FLOAT64 z;

	CDXFloat64Point3D(
		FLOAT64 x_arg = 0,
		FLOAT64 y_arg = 0,
		FLOAT64 z_arg = 0) :
	x(x_arg),
	y(y_arg),
	z(z_arg)
	{
	}

	CDXFloat64Point3D(const CDXFloat64Point3D& src)
	{
		x = src.x;  y = src.y;  z = src.z;
	}

	void FromString(const std::string &s)
	{
		std::istringstream is(s);
		is >> x >> y >> z;
	}
};

CDXFloat64Point3D StringToCDXFloat64Point3D(const std::string &s);

// ******************************
// ** class CDXFloat64Matrix3D  **
// ******************************
//
// CDXFloat64Matrix3D is a simple 4*4 matrix defined using FLOAT64

class CDXFloat64Matrix3D
{
public:
	FLOAT64 m[4][4];

	CDXFloat64Matrix3D()
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] = 0.0;

		for (int i = 0; i < 4; ++i)
			m[i][i] = 1.0;
	}

	CDXFloat64Matrix3D(const CDXFloat64Matrix3D& src)
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] = src.m[i][j];
	}

	void FromString(const std::string &s)
	{
		std::istringstream is(s);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				is >> m[i][j];
	}
};

CDXFloat64Matrix3D StringToCDXFloat64Matrix3D(const std::string &s);


// *************************
// ** class CDXRectangle  **
// *************************
//
// CDXRectangle is a simple 2 dimensional rectangle, used in the interface to the geometry
// of CDXGraphicObjects

class CDXRectangle
{
public:
	CDXCoordinate top;
	CDXCoordinate left;
	CDXCoordinate bottom;
	CDXCoordinate right;
	CDXRectangle(CDXCoordinate top_arg = 0,
				 CDXCoordinate left_arg = 0,
				 CDXCoordinate bottom_arg = 0,
				 CDXCoordinate right_arg = 0)
				: top(top_arg),
				  left(left_arg),
				  bottom(bottom_arg),
				  right(right_arg)
				{ }
	const CDXCoordinate	Width() const	{ return right - left; }
	const CDXCoordinate Height() const	{ return bottom - top; }
	bool				IsEmpty() const	{ return Width() == 0  ||  Height() == 0; }
	CDXRectangle	operator | (const CDXRectangle& op2) const
					{
						CDXRectangle	tmp (*this);
						tmp |= op2;
						return tmp;
					}
	const CDXRectangle&	operator |= (const CDXRectangle& op2)
					{
						left	= min (left,	op2.left);
						top		= min (top,		op2.top);
						right	= max (right,	op2.right);
						bottom	= max (bottom,	op2.bottom);
						return *this;
					}

	bool	Contains (const CDXRectangle &r2) const
					{
						return left <= r2.left  &&  right >= r2.right  &&  top <= r2.top  &&  bottom >= r2.bottom;
					}
	void	Inflate (const CDXCoordinate &dx, const CDXCoordinate &dy)
					{
						left -= dx;
						right += dx;
						top -= dy;
						bottom += dy;
					}
	void	Offset (const CDXCoordinate &dx, const CDXCoordinate &dy) 
					{ 
						left += dx;
						right += dx;
						top += dy;
						bottom += dy;
					}
};

CDXRectangle StringToCDXRectangle(const std::string &s);

// ***********************
// ** struct CDXPropRep **
// ***********************
//
// A CDXPropRep is stored for objects which represent properties of other objects.
// For example, a charge object might represent a charge property on a node.
// This is necessary so that the node doesn't try to display the charge itself.
#if TARGET_OS_WIN32
#pragma mark
#endif

class CDXPropRep
{
public:
	CDXObjectID objectID;
	CDXTag		propTag;
	CDXPropRep(CDXObjectID obj, CDXTag t) : objectID(obj), propTag(t) {}
	CDXPropRep() : objectID(kCDXUndefinedId), propTag(kCDXTag_UserDefined) {}
};

class CDXChemProp
{
public:
	CDXTag		datumID;
	CDXString	propString;
	CDXChemProp(CDXTag obj, CDXString t) : datumID(obj), propString(t) {}
	CDXChemProp() : datumID(kCDXTag_UserDefined) {}
};

// *************************
// ** class CDXDataSource **
// *************************
//
// CDXDataSource is an abstract interface for a data source.
// You'll probably want to multiply inherit from them and your IO object so you can
// pass your IO object to the CDX IO methods.  The GetBytes method must be overridden
// to get the bytes from wherever your data is coming from.

class CORE_CHEMISTRY_API CDXDataSource
{
public:
	CDXDataSource() : m_bigEndian(false) {}
	virtual ~CDXDataSource() {}

	void		SetBigEndian()	{ m_bigEndian = true; }

	// These are the only two functions which must be overridden to implement the data source.
	virtual void   GetBytes (char *p, size_t n) = 0;
	virtual void   SkipBytes (size_t n) = 0;

	// The rest of the member functions use GetBytes(char *, size_t) to get the bytes.
	// They should not be overridden.

	void   GetBytes  (unsigned char *p, size_t n) { GetBytes((char *)p, n); }

	// These use GetBytes to get various integral types, and perform byte swapping if necessary.
	// They should not be overridden.

	UINT8		GetUINT8  ()			{ UINT8 p; GetBytes(&p, 1); return p; }
	UINT16		GetUINT16 ();
	UINT32		GetUINT32 ();

	INT8		GetINT8  ()		{ return INT8 (GetUINT8 ()); };
	INT16		GetINT16 ()		{ return INT16(GetUINT16()); };
	INT32		GetINT32 ()		{ return INT32(GetUINT32()); };

	FLOAT64		GetFLOAT64();

	// These read the specified number of bytes
	UINT32		GetUINT(size_t size_arg);
	INT32		GetINT(size_t size_arg);

	// This reads n bytes as a standard C++ string
	std::string	GetString(size_t n);

	CDXTag		GetTag ()		{ return GetUINT16(); };
	CDXObjectID GetObjectID ()	{ return GetUINT32(); };

	virtual bool IgnoreTag(CDXTag tag, const CDXObject *parent) const { return false; }
	virtual bool IgnoreObject(CDXTag objTag, const CDXObject *parent) const { return false; }
	
private:
	bool		m_bigEndian;
};

inline CDXDataSource& operator>> (CDXDataSource& src_arg, INT8   &value_arg) { value_arg = src_arg.GetINT8();   return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, INT16  &value_arg) { value_arg = src_arg.GetINT16();  return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, INT32  &value_arg) { value_arg = src_arg.GetINT32();  return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, UINT8  &value_arg) { value_arg = src_arg.GetUINT8();  return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, UINT16 &value_arg) { value_arg = src_arg.GetUINT16(); return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, UINT32 &value_arg) { value_arg = src_arg.GetUINT32(); return src_arg; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, CDXPoint2D &value_arg) { return src_arg >> value_arg.y >> value_arg.x; }
inline CDXDataSource& operator>> (CDXDataSource& src_arg, CDXPoint3D &value_arg) { return src_arg >> value_arg.x >> value_arg.y >> value_arg.z; }

inline CDXDataSource& operator>> (CDXDataSource& src_arg, CDXFloat64Point3D &value_arg)
{
	value_arg.x = src_arg.GetFLOAT64();
	value_arg.y = src_arg.GetFLOAT64();
	value_arg.z = src_arg.GetFLOAT64();
	return src_arg;
}

inline CDXDataSource& operator>> (CDXDataSource& src_arg, CDXFloat64Matrix3D &value_arg)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			value_arg.m[i][j] = src_arg.GetFLOAT64();
	return src_arg;
}

// ***********************
// ** class CDXistream **
// ***********************
//
// CDXistream is an implementation of CDXDataSource which you can use if your
// input is coming from a C++ istream.  With this, one can also easily do:
//		CDXistream inp(istrstream(ptr, siz));
//	or	CDXistream inp(istringstream(str));
//	or	CDXistream inp(ifstream(infile, ios::binary));
// to get your input from a blob of memory or a C++ string or a stdio FILE *.

class CORE_CHEMISTRY_API CDXistream : public CDXDataSource
{
	std::istream &is;
public:
	CDXistream(std::istream &input) : is(input) {}
	virtual void GetBytes  (char *p, size_t n);
	virtual void SkipBytes  (size_t n);
};

// ***********************
// ** class CDXipointer **
// ***********************
//
// CDXipointer is an implementation of CDXDataSource which you can use if your
// input is coming from an in-memory pointer.

class CDXipointer : public CDXDataSource
{
	const char *m_pos;
	const char *m_end;
public:
	CDXipointer(const char *pos, size_t n) : m_pos(pos), m_end(pos + n) {}
	virtual void GetBytes  (char *p, size_t n);
	virtual void SkipBytes  (size_t n);
};

// ***********************
// ** class CDXDataSink **
// ***********************
//
// CDXDataSink is the abstract interface for the data sink.
// You'll probably want to multiply inherit from this and your output
// stream object(s).  The PutBytes method must be overridden to put
// the output where you want it.

class CORE_CHEMISTRY_API CDXDataSink
{
public:
	virtual ~CDXDataSink(){}
	// This is the one and only function which must be overridden to implement the data sink.
	virtual void Put ( const INT8 *data, size_t n ) = 0;

	// These use Put(INT8 *, size_t) to write various other types, and perform byte swapping if necessary.
	// They should not be overridden.
	void Put ( const UINT8 *data, size_t n ) { Put((const INT8 *)data, n); }

	void Put ( UINT8  datum ) { Put(&datum, 1); }
	void Put ( UINT16 datum );
	void Put ( UINT32 datum );
	void Put ( FLOAT64 datum );

	void Put ( INT8  datum ) { Put( UINT8 (datum) ); }
	void Put ( INT16 datum ) { Put( UINT16(datum) ); }
	void Put ( INT32 datum ) { Put( UINT32(datum) ); }
	void Put ( const CDXPoint2D &datum ) { Put(datum.y); Put(datum.x); }
	void Put ( const CDXPoint3D &datum ) { Put(datum.x); Put(datum.y); Put(datum.z); }

	void Put ( const CDXFloat64Point3D &datum )
	{
		Put(datum.x);
		Put(datum.y);
		Put(datum.z);
	}

	void Put ( const CDXFloat64Matrix3D &datum )
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				Put(datum.m[i][j]);
	}

#if !defined(CORE_CHEMISTRY_API_DLL_BUILD)
    void Put(const C3XColorRGBA &datum);  // TODO: C3XColorRGBA leakage from Chem3D into Core Chemistry.
#endif

	void PutTag      ( CDXTag      datum ) { Put( UINT16(datum) ); }
	void PutObjectID ( CDXObjectID datum ) { Put( UINT32(datum) ); }

	// Put an object header
	void PutObject   ( CDXTag tag, CDXObjectID id ) { PutTag(tag); PutObjectID(id); }

	// Put a complete attribute, consisting of the tag, length, and data
	// Someday some of this will be a template member function
	void PutAttribute ( CDXTag tag ) { PutTag(tag); Put(UINT16(0)); }
	void PutAttribute ( CDXTag tag, UINT8  datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, UINT16 datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, UINT32 datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, INT8   datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, INT16  datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, INT32  datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, FLOAT32 datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, FLOAT64 datum ) { PutTag(tag); Put(UINT16(sizeof datum)); Put(datum); }
	void PutAttribute ( CDXTag tag, const CDXPoint2D &datum )
		{ PutTag(tag); Put(UINT16(sizeof datum.x + sizeof datum.y)); Put(datum); }
	void PutAttribute ( CDXTag tag, const CDXPoint3D &datum )
		{ PutTag(tag); Put(UINT16(sizeof datum.x + sizeof datum.y + sizeof datum.z)); Put(datum); }
	void PutAttribute ( CDXTag tag, const CDXFloat64Point3D &datum )
		{ PutTag(tag); Put(UINT16(3 * sizeof(double))); Put(datum); }
	void PutAttribute ( CDXTag tag, const CDXFloat64Matrix3D &datum )
		{ PutTag(tag); Put(UINT16(16 * sizeof(double))); Put(datum); }

#if !defined(CORE_CHEMISTRY_API_DLL_BUILD)
    void PutAttribute ( CDXTag tag, const C3XColorRGBA &datum ) // TODO: C3XColorRGBA leakage from Chem3D into Core Chemistry.
        { PutTag(tag); Put(UINT16(4 * sizeof(UINT8))); Put(datum); }
#endif

    void PutAttribute ( CDXTag tag, const Point32 &datum )
		{ PutTag(tag); Put(UINT16(sizeof datum.x + sizeof datum.y)); Put(datum.y); Put(datum.x); }
	void PutAttribute ( CDXTag tag, const CDXRectangle &datum )
		{ PutTag(tag); Put(UINT16(4 * sizeof datum.top));
		  Put(datum.top); Put(datum.left); Put(datum.bottom); Put(datum.right); }
	void PutAttribute ( CDXTag tag, const std::vector<CDXPropRep> &datum )
		{ PutTag(tag);
		  Put(UINT16(datum.size() * (sizeof datum[0].objectID + sizeof datum[0].propTag)));
		  for (std::vector<CDXPropRep>::const_iterator i = datum.begin(); i != datum.end();  ++i)
		  {
			  Put(i->objectID);
			  Put(i->propTag);
		  }
		}
	void PutAttribute ( const std::vector<CDXChemProp> &datum )
	{ 
		for (std::vector<CDXChemProp>::const_iterator i = datum.begin(); i != datum.end();  ++i)
		{
			PutAttributeCDXString((*i).datumID, (*i).propString);

		}
	}
	void PutAttribute ( CDXTag tag, const CDXFontStyle &datum )
		{ PutTag(tag); Put(UINT16(sizeof datum.family + sizeof datum.face + sizeof datum.size + sizeof datum.color));
		  Put(datum.family); Put(datum.face); Put(datum.size); Put(datum.color); }
	virtual void PutAttributeCDXString ( CDXTag tag, const CDXString &datum )
		{ PutTag( tag );
		  Put( UINT16( datum.CDXDataSize() ) );
		  datum.PutTo(*this); 
	}
	void PutAttributeString ( CDXTag tag, const std::string &datum )
		{ PutTag( tag );
		  Put( UINT16( datum.size() ) );
		  Put((INT8 *)datum.data(), datum.size()); }
	void PutAttribute ( CDXTag tag, const INT8 *data, size_t n )
		{ PutTag(tag); 
		  if (sizeof(*data) * n < 0xFFFF) 
			  Put(UINT16(sizeof(*data) * n)); 
		  else 
		  {
			  Put(INT16(0xFFFF));
			  Put(UINT32(sizeof(*data) * n)); 
		  }
		  Put(data, n); }
	void PutAttribute ( CDXTag tag, const UINT32 *data, size_t n )
		{ PutTag(tag); 
		  if (sizeof(*data) * n < 0xFFFF) 
			  Put(UINT16(sizeof(*data) * n)); 
		  else 
		  {
			  Put(INT16(0xFFFF));
			  Put(UINT32(sizeof(*data) * n)); 
		  }
		  while (n-- > 0) Put(*(data++)); }
	void PutAttribute ( CDXTag tag, const FLOAT64 *data, size_t n )
		{ PutTag(tag); 
		  if (sizeof(*data) * n < 0xFFFF) 
			  Put(UINT16(sizeof(*data) * n)); 
		  else 
		  {
			  Put(INT16(0xFFFF));
			  Put(UINT32(sizeof(*data) * n)); 
		  }
		  while (n-- > 0) Put(*(data++)); }

    void PutAttributeForObjectIDList(CDXTag tag, const std::vector<CDXObjectID>& objectIDList);

	void PutFileHeader ()
		{ Put((INT8 *)kCDX_HeaderString, 8);
		  Put(INT32(0x01020304L));
		  Put(INT32(0));
		  Put(INT32(0));
		  Put(INT16(0)); }
};

// This ought to be a template member function, too.
template <class InputIterator>
inline void Put (CDXDataSink &sink, InputIterator first, InputIterator last)
{
	while (first != last)
		sink.Put(*(first++));
}

template <class T>
inline void Put (CDXDataSink &sink, const std::vector<T> &c)
{
	Put(sink, c.begin(), c.end());
}

inline CDXDataSink& operator<< (CDXDataSink& sink_arg, INT8   value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, INT16  value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, INT32  value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, UINT8  value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, UINT16 value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, UINT32 value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, const CDXPoint2D &value_arg) { sink_arg.Put(value_arg); return sink_arg; }
inline CDXDataSink& operator<< (CDXDataSink& sink_arg, const CDXPoint3D &value_arg) { sink_arg.Put(value_arg); return sink_arg; }


// **********************
// ** class CDXostream **
// **********************
//
// CDXostream is an implementation of CDXDataSink which you can use if your
// output is going to a C++ ostream.  With this, one can also easily do:
//		CDXostream s(ostrstream(ptr, siz));
//	or	CDXostream s(ostringstream(outstr));
//	or	CDXostream s(ofstream(outfile, ios::out | ios::trunc | ios::binary));
// to send your output to a blob of memory or a C++ string, or to a stdio FILE *.

class CDXostream : public CDXDataSink
{
	std::ostream &os;
public:
	CDXostream(std::ostream &output) : os(output) {}
	
	CORE_CHEMISTRY_API void Put ( const INT8 *data, size_t n );
};

// ***********************
// ** class XMLDataSink **
// ***********************
//

class CORE_CHEMISTRY_API XMLDataSink
{
public:
	std::ostream &os;
	XMLDataSink(std::ostream &s) : os(s) {}
};

// *************************
// ** Low-level utilities **
// *************************

// ** function CDXIsObjectTag **
// Determine whether a particular tag represents an object or attribute
inline bool CDXIsObjectTag(CDXTag tag)
{
	return (tag & kCDXTag_Object) != 0;
}

// *********************************************************
// ** ostream& operator<<(ostream& os, const vector<T> &v) **
// *********************************************************
//
// Sends a vector of T's to an output stream, separated by spaces
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
	for (typename vector<T>::const_iterator i = v.begin();  i != v.end();  ++i)
	{
		if (i != v.begin())
			os << " ";
		os << *i;
	}
	return os;
}

// These additional output stream operators are used to dump various objects
std::ostream& operator<<(std::ostream& os, const CDXFontStyle &s);
std::ostream& operator<<(std::ostream& os, const CDXStyle &s);
CORE_CHEMISTRY_API std::ostream& operator<<(std::ostream& os, const CDXString &s);
std::ostream& operator<<(std::ostream& os, CDXAtomGeometry s);
std::ostream& operator<<(std::ostream& os, const CDXPoint2D &p);
std::ostream& operator<<(std::ostream& os, const CDXPoint3D &p);
std::ostream& operator<<(std::ostream& os, const CDXRectangle &r);
CORE_CHEMISTRY_API std::ostream& operator<<(std::ostream& os, CDXDatumID tag);

// ** class CDXDataSourceArray<T> **
// This is a helper class, used to turn a CDXDataSource into a generator function
// for a given type.  This is useful when we want to read N of some type T
// into a container C using code like:
//    generate_n(back_inserter(C), N, CDXDataSourceArray<T>(src));
template<class T>
class CDXDataSourceArray {
	CDXDataSource &src;
public:
	CDXDataSourceArray(CDXDataSource &s) : src(s) {}
	T operator()() { T inp; src >> inp; return inp; }
};

template<class Container>
void CDXReadItems(Container &container_arg, size_t nItems_arg, CDXDataSource &src_arg)
{
	container_arg.reserve(nItems_arg);
	std::generate_n(std::back_inserter(container_arg), nItems_arg, CDXDataSourceArray<typename Container::value_type>(src_arg));
}

// ** class CDXDataSinkIterator<T> **
// This is a helper class, used to turn a CDXDataSink into an output iterator
// for a given type.  This is useful when we want to write a bunch of some type T
// from a container C using code like:
//    copy(C.begin(), C.end(), CDXDataSinkIterator<T>(sink));
// It relies on the existence of operator<<(CDXDataSink&, const T&).
template<class T>
class CDXDataSinkIterator
{
private:
    
	CDXDataSink *m_sink;

public:

    using iterator_category = std::output_iterator_tag;
    using value_type = T;
    using difference_type = T;
    using pointer = T*;
    using reference = T&;

    CDXDataSinkIterator(CDXDataSink &s) : m_sink(&s) 
    {
    }

    CDXDataSinkIterator<T>& operator=(const T& value) 
    { 
        *m_sink << value;
        return *this;
    }

    CDXDataSinkIterator<T>& operator*()     { return *this; }
    CDXDataSinkIterator<T>& operator++()    { return *this; } 
    CDXDataSinkIterator<T>& operator++(int) { return *this; } 
};

void CDXWriteHeader (CDXDataSink &sink);

