// CommonCS/LibCommon/Hdr/CDXFontTable.h
// Copyright © 1999-2004, CambridgeSoft Corp., All Rights Reserved

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

#include "CoreChemistryAPI.h"
#include "CDXIO.h"
#include "CDXUnicode.h"
#if defined _CFW_  &&  !defined NO_CDXML_PARSER
#	define NO_CDXML_PARSER 1
#endif
#ifndef NO_CDXML_PARSER
#include "CDXMLParser.h"
#endif
#include <string>
#include <vector>

// *******************************
// ***** class CDXFontTable ******
// *******************************
//
// Helper class for manipulating a CDX font table

struct CORE_CHEMISTRY_API CDXFontTableEntry
{
	CDXFontTableEntry() : m_internalFontNum(0), m_type(CDXCharSet::kCDXCharSetUnknown) {}

	UINT16		m_internalFontNum;
	CDXCharSet	m_type;
	std::string	m_family;
};

inline bool operator==(const CDXFontTableEntry &f, int i) { return f.m_internalFontNum == i; }
inline bool operator==(const CDXFontTableEntry &f, const std::string &s) { return f.m_family == s; }

class CORE_CHEMISTRY_API CDXFontTable
{
public:
	std::vector<CDXFontTableEntry>	m_fonts;

	CDXFontTable () {}
	
	void Read(CDXDataSource &src_arg, size_t size_arg);
		// Read from binary data stream

	void Write(CDXDataSink &sink_arg) const;
		// Write to binary data stream

	void XMLRead(const std::string &);
		// Read from XML attribute string

	void XMLWrite(XMLDataSink &) const;
		// Write as XML attribute string

	void push_back(const CDXFontTableEntry &f)
			{ m_fonts.push_back(f); }

	int NumFonts() const
			{ return (int)m_fonts.size(); }
		// Return the current size of the font table

	const CDXFontTableEntry &LookupFont(int family) const;

    CDXCharSet EncodingForFontID(int inFamily) const;
};

#ifndef NO_CDXML_PARSER

// *************************
// class CDXMLFontHandler
//
// element handler for font elements
// *************************

class CDXMLFontHandler : public XMLElementHandler
{
public:
	CDXMLFontHandler(CDXMLParser *parser, const XML_Char **atts);
	~CDXMLFontHandler() {}
};

// *************************
// class CDXMLFontTableHandler
//
// element handler for fonttable elements
// *************************

class CORE_CHEMISTRY_API CDXMLFontTableHandler : public XMLElementHandler
{
	CDXDocument *m_document;

public:
	~CDXMLFontTableHandler() {}
	CDXMLFontTableHandler(CDXMLParser *parser, const XML_Char **atts);
	void startElement(const XML_Char* name, const XML_Char** atts);
};

#endif // NO_CDXML_PARSER

