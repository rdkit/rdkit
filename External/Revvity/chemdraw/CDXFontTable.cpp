// CommonCS/LibCommon/Src/CDXFontTable.cpp
// Contains: Program-independent class library for managing CDX objects
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

#include "CDXStdObjects.h"
#include "CDXFontTable.h"
#include "CDXMLNames.h"
#include "CDXUnicode.h"
#include "XMLPutEnum.h"
#include "cs_mapUtils.h"
#include <stdlib.h>	// for atoi and friends
#include <string.h>	// for strncmp and friends
#include <sstream>
#include <stdexcept>

XMLPutEnum<CDXCharSet>::Value sXMLCharSetValues[] = {
	{kCDXCharSetUnknown,				"Unknown"},
	{kCDXCharSetEBCDICOEM,				"EBCDICOEM"},
	{kCDXCharSetMSDOSUS,				"MSDOSUS"},
	{kCDXCharSetEBCDIC500V1,			"EBCDIC500V1"},
	{kCDXCharSetArabicASMO708,			"ASMO-708"},
	{kCDXCharSetArabicASMO449P,			"ArabicASMO449P"},
	{kCDXCharSetArabicTransparent,		"ArabicTransparent"},
	{kCDXCharSetArabicTransparentASMO,	"DOS-720"},
	{kCDXCharSetGreek437G,				"Greek437G"},
	{kCDXCharSetBalticOEM,				"cp775"},
	{kCDXCharSetMSDOSLatin1,			"windows-850"},
	{kCDXCharSetMSDOSLatin2,			"ibm852"},
	{kCDXCharSetIBMCyrillic,			"cp855"},
	{kCDXCharSetIBMTurkish,				"cp857"},
	{kCDXCharSetMSDOSPortuguese,		"cp860"},
	{kCDXCharSetMSDOSIcelandic,			"cp861"},
	{kCDXCharSetHebrewOEM,				"DOS-862"},
	{kCDXCharSetMSDOSCanadianFrench,	"cp863"},
	{kCDXCharSetArabicOEM,				"cp864"},
	{kCDXCharSetMSDOSNordic,			"cp865"},
	{kCDXCharSetMSDOSRussian,			"cp866"},
	{kCDXCharSetIBMModernGreek,			"cp869"},
	{kCDXCharSetThai,					"windows-874"},
	{kCDXCharSetEBCDIC,					"EBCDIC"},
	{kCDXCharSetJapanese,				"shift_jis"},
	{kCDXCharSetChineseSimplified,		"gb2312"}, // PRC, Singapore
	{kCDXCharSetKorean,					"ks_c_5601-1987"},
	{kCDXCharSetChineseTraditional,		"big5"}, // Taiwan, Hong Kong
	{kCDXCharSetUnicodeISO10646,		"iso-10646"},
	{kCDXCharSetWin31EasternEuropean,	"windows-1250"},
	{kCDXCharSetWin31Cyrillic,			"windows-1251"},
	{kCDXCharSetWin31Latin1,			"iso-8859-1"},
	{kCDXCharSetWin31Greek,				"iso-8859-7"},
	{kCDXCharSetWin31Turkish,			"iso-8859-9"},
	{kCDXCharSetHebrew,					"windows-1255"},
	{kCDXCharSetArabic,					"windows-1256"},
	{kCDXCharSetBaltic,					"windows-1257"},
	{kCDXCharSetVietnamese,				"windows-1258"},
	{kCDXCharSetKoreanJohab,			"windows-1361"},
	{kCDXCharSetMacRoman,				"x-mac-roman"},
	{kCDXCharSetMacJapanese,			"x-mac-japanese"},
	{kCDXCharSetMacTradChinese,			"x-mac-tradchinese"},
	{kCDXCharSetMacKorean,				"x-mac-korean"},
	{kCDXCharSetMacArabic,				"x-mac-arabic"},
	{kCDXCharSetMacHebrew,				"x-mac-hebrew"},
	{kCDXCharSetMacGreek,				"x-mac-greek"},
	{kCDXCharSetMacCyrillic,			"x-mac-cyrillic"},
	{kCDXCharSetMacReserved,			"x-mac-reserved"},
	{kCDXCharSetMacDevanagari,			"x-mac-devanagari"},
	{kCDXCharSetMacGurmukhi,			"x-mac-gurmukhi"},
	{kCDXCharSetMacGujarati,			"x-mac-gujarati"},
	{kCDXCharSetMacOriya,				"x-mac-oriya"},
	{kCDXCharSetMacBengali,				"x-mac-nengali"},
	{kCDXCharSetMacTamil,				"x-mac-tamil"},
	{kCDXCharSetMacTelugu,				"x-mac-telugu"},
	{kCDXCharSetMacKannada,				"x-mac-kannada"},
	{kCDXCharSetMacMalayalam,			"x-mac-Malayalam"},
	{kCDXCharSetMacSinhalese,			"x-mac-sinhalese"},
	{kCDXCharSetMacBurmese,				"x-mac-burmese"},
	{kCDXCharSetMacKhmer,				"x-mac-khmer"},
	{kCDXCharSetMacThai,				"x-mac-thai"},
	{kCDXCharSetMacLao,					"x-mac-lao"},
	{kCDXCharSetMacGeorgian,			"x-mac-georgian"},
	{kCDXCharSetMacArmenian,			"x-mac-armenian"},
	{kCDXCharSetMacSimpChinese,			"x-mac-simpChinese"},
	{kCDXCharSetMacTibetan,				"x-mac-tibetan"},
	{kCDXCharSetMacMongolian,			"x-mac-mongolian"},
	{kCDXCharSetMacEthiopic,			"x-mac-ethiopic"},
	{kCDXCharSetMacCentralEuroRoman,	"x-mac-ce"},
	{kCDXCharSetMacVietnamese,			"x-mac-vietnamese"},
	{kCDXCharSetMacExtArabic,			"x-mac-extArabic"},
	{kCDXCharSetMacUninterpreted,		"x-mac-uninterpreted"},
	{kCDXCharSetMacSymbol,				"x-mac-symbol"},
	{kCDXCharSetMacDingbats,			"x-mac-dingbats"},
	{kCDXCharSetMacCroatian,			"x-mac-croatian"},
	{kCDXCharSetMacRomanian,			"x-mac-romanian"},
	{kCDXCharSetMacCeltic,				"x-mac-celtic"},
	{kCDXCharSetMacGaelic,				"x-mac-gaelic"},
	{kCDXCharSetMacKeyboardGlyphs,		"x-mac-keyboard-glyphs"},
	{kCDXCharSetMacIcelandic,			"x-mac-icelandic"},
    {kCDXCharSetMacTurkish,				"x-mac-turkish"},
    {kCDXCharSetUTF8,                   "utf-8"}
};

XMLPutEnum<CDXCharSet> sXMLCharSet(sXMLCharSetValues, sizeof sXMLCharSetValues, kCDXCharSetUnknown);
void XMLPut(XMLDataSink &sink, CDXCharSet  v);
void XMLPut(XMLDataSink &sink, CDXCharSet  v) {	sink.os << sXMLCharSet.lookup(v); }

static bool ConvertUnreadableFontName(string &font)
{
    // These Japanese font names were written incorrectly in older versions.
    // Here we convert the incorrect character sequences to the correct Mac UTF8 sequence for the affected font names.
    // The map is: bad sequence / good sequence
    // The sequences are not stored here in UTF8 as that would require saving the source file in that format.
    // Note:  Some sequence strings are split and rejoined from strings due to a Mac lexical issue.
    static cs::csInitMap<std::string, std::string> smMacFontNameMap = cs::csInitMap<std::string, std::string>()
        << std::make_pair("\x82l\x82r \x96\xbe\x92\xa9",
                          "\xef\xbc\xad\xef\xbc\xb3 \xe6\x98\x8e\xe6\x9c\x9d") // MS Mincho
        << std::make_pair("\x82l\x82r \x82o\x96\xbe\x92\xa9",
                          "\xef\xbc\xad\xef\xbc\xb3 \xef\xbc\xb0\xe6\x98\x8e\xe6\x9c\x9d") // MS Pincho
        << std::make_pair(std::string("\x82l\x82r \x83S\x83V\x83") + std::string("b\x83N"),
                          "\xef\xbc\xad\xef\xbc\xb3 \xe3\x82\xb4\xe3\x82\xb7\xe3\x83\x83\xe3\x82\xaf") // MS Mincho Bold
        << std::make_pair(std::string("\x82l\x82r \x82o\x83S\x83V\x83") + std::string("b\x83N"),
                          "\xef\xbc\xad\xef\xbc\xb3 \xef\xbc\xb0\xe3\x82\xb4\xe3\x82\xb7\xe3\x83\x83\xe3\x82\xaf") // MS Pincho Bold
        << std::make_pair("Osaka\x81|\x93\x99\x95\x9d",
                          "Osaka\xe2\x88\x92\xe7\xad\x89\xe5\xb9\x85")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S Pro W3",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 Pro W3")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8a\xdb\x83S Pro W4",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe4\xb8\xb8\xe3\x82\xb4 Pro W4")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S Pro W6",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 Pro W6")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S ProN W3",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 ProN W3")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8a\xdb\x83S ProN W4",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe4\xb8\xb8\xe3\x82\xb4 ProN W4")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S ProN W6",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 ProN W6")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S Std W8",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 Std W8")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x8ap\x83S StdN W8",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe8\xa7\x92\xe3\x82\xb4 StdN W8")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x96\xbe\x92\xa9 Pro W3",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe6\x98\x8e\xe6\x9c\x9d Pro W3")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x96\xbe\x92\xa9 Pro W6",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe6\x98\x8e\xe6\x9c\x9d Pro W6")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x96\xbe\x92\xa9 ProN W3",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe6\x98\x8e\xe6\x9c\x9d ProN W3")
        << std::make_pair("\x83q\x83\x89\x83M\x83m\x96\xbe\x92\xa9 ProN W6",
                          "\xe3\x83\x92\xe3\x83\xa9\xe3\x82\xae\xe3\x83\x8e\xe6\x98\x8e\xe6\x9c\x9d ProN W6")
        << std::make_pair(std::string("\x83\x81\x83") + std::string("C\x83\x8a\x83I"),
                          "\xe3\x83\xa1\xe3\x82\xa4\xe3\x83\xaa\xe3\x82\xaa")
        << std::make_pair(std::string("\x83\x81\x83") + std::string("C\x83\x8a\x83I \x83") + std::string("C\x83^\x83\x8a\x83") + std::string("b\x83N"),
                          "\xe3\x83\xa1\xe3\x82\xa4\xe3\x83\xaa\xe3\x82\xaa \xe3\x82\xa4\xe3\x82\xbf\xe3\x83\xaa\xe3\x83\x83\xe3\x82\xaf")
        << std::make_pair(std::string("\x83\x81\x83") + std::string("C\x83\x8a\x83I \x83{\x81[\x83\x8b\x83h \x83") + std::string("C\x83^\x83\x8a\x83") + std::string("b\x83N"),
                          "\xe3\x83\xa1\xe3\x82\xa4\xe3\x83\xaa\xe3\x82\xaa \xe3\x83\x9c\xe3\x83\xbc\xe3\x83\xab\xe3\x83\x89 \xe3\x82\xa4\xe3\x82\xbf\xe3\x83\xaa\xe3\x83\x83\xe3\x82\xaf")
        << std::make_pair(std::string("\x83\x81\x83") + std::string("C\x83\x8a\x83I \x83{\x81[\x83\x8b\x83h"),
                          "\xe3\x83\xa1\xe3\x82\xa4\xe3\x83\xaa\xe3\x82\xaa \xe3\x83\x9c\xe3\x83\xbc\xe3\x83\xab\xe3\x83\x89");

    if (smMacFontNameMap.find(font) != smMacFontNameMap.end())
    {
        font = smMacFontNameMap[font];
        return true;
    }
    
    return false;
}

// *******************************
// ***** class CDXFontTable *****
// *******************************
//

void CDXFontTable::Read(CDXDataSource &src_arg, size_t size_arg)
{
	if (size_arg < 4)
    {
		throw std::runtime_error ("Bad font table");
    }

	(void) src_arg.GetUINT16();	// used to be platform, but we don't use it anymore
	UINT16 numFonts = src_arg.GetUINT16();

	int	numBytesRead = 4;
	for (int i = 0;  i < numFonts;  i++)
	{
		CDXFontTableEntry fte;
		fte.m_internalFontNum = src_arg.GetUINT16();
		fte.m_type = (CDXCharSet) src_arg.GetUINT16();
		size_t nameSize = src_arg.GetUINT16();
        std::string family = src_arg.GetString(nameSize);

        #if UTF8_STD_STRING
            // Guess that our font name is in the same character set as the font, and transcode from
            // that character set to UTF-8
            family = TranscodeString(family, fte.m_type, kCDXCharSetUTF8);
        #endif

        family = UnicodeToString(family);
        ConvertUnreadableFontName(family);

        fte.m_family.assign(family);
       
		m_fonts.push_back(fte);
		numBytesRead += 6 + nameSize;
	}

	if (numBytesRead != size_arg)
    {
		throw std::runtime_error ("Bad font table");
    }
}

void CDXFontTable::Write(CDXDataSink &sink_arg) const
{
	sink_arg.Put(INT16(kCDXProp_FontTable));

    // Create a copy of the map that contains font names in the format that we will write out
    vector<CDXFontTableEntry> fontMapCopy;

    for (const auto& fontEntry : m_fonts)
    {
        CDXFontTableEntry newEntry = fontEntry;
        newEntry.m_family = StringToUnicode(newEntry.m_family, newEntry.m_family);

        #if UTF8_STD_STRING
            // Transcode the font name from UTF-8 to the encoding that we store with the font. This allows
            // backwards-compatibility with non-UTF-8 versions of ChemDraw.
            newEntry.m_family = TranscodeString(newEntry.m_family, kCDXCharSetUTF8, newEntry.m_type);
        #endif

        fontMapCopy.push_back(newEntry);
    }

	// Add up the size of the attribute. We start on 4 to account for platform and numFonts.
	size_t numBytes = 4;

    for (const auto& fontEntry : fontMapCopy)
    {
        // internalFontNum, type, nameLength, name
        numBytes += 6 + fontEntry.m_family.size();
    }

    // Now put the numBytes, platform, numFonts
	sink_arg.Put(INT16(numBytes));

#if TARGET_OS_MAC
	sink_arg.Put(INT16(1)); // 1 for mac, 0 for windows, but we don't use this anymore
#else // TARGET_OS_MAC
	sink_arg.Put(INT16(0)); // 1 for mac, 0 for windows, but we don't use this anymore
#endif // TARGET_OS_MAC

	sink_arg.Put(INT16(fontMapCopy.size()));
	
	// Finally, each font has a fontnum, type, the name length, and the name
    for (const auto& fontEntry : fontMapCopy)
    {
        sink_arg.Put(INT16(fontEntry.m_internalFontNum));
        sink_arg.Put(INT16(fontEntry.m_type));
        sink_arg.Put(INT16(fontEntry.m_family.size()));
        sink_arg.Put((INT8 *)fontEntry.m_family.data(), fontEntry.m_family.size());
    }
}

std::ostream & operator<<(std::ostream &os, const CDXFontTableEntry &f)
{
	os << GetTextEOL() << "<" << kCDXML_fontElement
		<< " " << kCDXML_id << "=\"" << f.m_internalFontNum << "\""
		<< " " << kCDXML_charset << "=\"" << sXMLCharSet.lookup(f.m_type) << "\""
        << " " << kCDXML_name << "=\"" << StringToUnicode("", f.m_family) << "\""
		<< "/>";

    return os;
}

void CDXFontTable::XMLWrite(XMLDataSink &sink) const
{
	sink.os << "<" << kCDXML_fonttable << ">";
	std::copy(m_fonts.begin(), m_fonts.end(), std::ostream_iterator<CDXFontTableEntry>(sink.os));
	sink.os << GetTextEOL() << "</" << kCDXML_fonttable << ">";
}

const CDXFontTableEntry &CDXFontTable::LookupFont(int i) const
{
	std::vector<CDXFontTableEntry>::const_iterator it = std::find(m_fonts.begin(), m_fonts.end(), i);
	if (it == m_fonts.end())
		throw std::runtime_error("font not found in font table");
	return *it;
}

/**
 *  Return the encoding for the given font ID. Returns kCDXCharSetUnknown if the font is
 *  unknown.
 *
 *  @param inFamily Font family ID
 *
 *  @return Encoding, or kCDXCharSetUnknown if unknown
 */
CDXCharSet CDXFontTable::EncodingForFontID(int inFamily) const
{
    CDXCharSet encoding = kCDXCharSetUnknown;
    
    std::vector<CDXFontTableEntry>::const_iterator it = std::find(m_fonts.begin(), m_fonts.end(), inFamily);
    if (it != m_fonts.end())
    {
        encoding = it->m_type;
    }
    
    return encoding;
}


#ifndef NO_CDXML_PARSER
// *************************
// class CDXMLFontTableHandler
//
// element handler for fonttable elements
// *************************
CDXMLFontTableHandler::CDXMLFontTableHandler(CDXMLParser *parser, const XML_Char **atts)
	: XMLElementHandler(parser)
	, m_document(dynamic_cast<CDXDocument *>(Parser()->CurrentObject()))
{
	if (m_document == 0)
		throw std::runtime_error("fonttable not contained in a document");
}

void CDXMLFontTableHandler::startElement(const XML_Char* name, const XML_Char** atts)
{
	XMLElementHandler* handler;
	CDXMLParser* parser;
	if (strcmp(name, kCDXML_fontElement) == 0 && ((parser = dynamic_cast<CDXMLParser*>(Parser())) != NULL) )
		handler = new CDXMLFontHandler(parser, atts);
	else
		handler = new XMLUnknownHandler(Parser(), atts);

	Parser()->PushHandler(handler);
}

// *************************
// class CDXMLFontHandler
//
// element handler for font elements
// *************************
CDXMLFontHandler::CDXMLFontHandler(CDXMLParser *parser, const XML_Char **atts)
	: XMLElementHandler(parser)
{
	CDXFontTableEntry f;
	f.m_internalFontNum = 0;
	f.m_type = kCDXCharSetUnknown;
	for (const XML_Char **p = atts;  *p != 0;  p += 2)
	{
		if (strcmp(p[0], kCDXML_id) == 0)
			f.m_internalFontNum = atoi(p[1]);
		else if (strcmp(p[0], kCDXML_charset) == 0)
		{
			try {
				f.m_type = sXMLCharSet.lookup(p[1]);
			}
			catch (...) { // It might be a number
				f.m_type = CDXCharSet(atoi(p[1]));
			}
		}
		else if (strcmp(p[0], kCDXML_name) == 0)
        {
            f.m_family = UnicodeToString(p[1]);
#if TARGET_OS_MAC
            if (!ConvertUnreadableFontName(f.m_family))
            {
                // Might be using a non-English primary language, in which case we have to try decoding the string explicitly
                CDXCharSet actualCharSet = kCDXCharSetWin31Latin1;
                std::string family = UnicodeToString("", &actualCharSet, p[1]);
                if (ConvertUnreadableFontName(family))
                {
                    f.m_family = family;
                }
            }
#endif // TARGET_OS_MAC
        }
	}
	parser->GetDocument()->m_fontTable.push_back(f);
}

#endif // NO_CDXML_PARSER
