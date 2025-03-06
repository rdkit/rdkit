// CommonCS/LibCommon/Src/CDXText.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2004, CambridgeSoft Corp., All Rights Reserved

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
#include "XMLPutEnum.h"
#include "CDXFontTable.h"
#include "CDXMLNames.h"
#include "UTF8Iterator.h"
#include "cs_charUtils.h"
#include "cs_specialchardefs.h"
#include "cs_assert.h"
#include <stdlib.h>	// for atoi and friends
#include <sstream>
#include <list>

extern void XMLPut(XMLDataSink &sink, CDXLabelDisplay  v);
extern XMLPutEnum<CDXLabelDisplay> sXMLLabelDisplay;
static void XMLWriteCDXString(XMLDataSink &sink_arg, const CDXString &theText);

namespace
{
    /**
     *  Convert text and styles a run at a time from one encoding to another
     *
     *  @param inString string to convert
     *  @param inConverter converter function
     *  @return converted string
     */
    CDXString ConvertTextEncodingWithConverter(const CDXString& inString, const std::function<string(const string&, const CDXStyle&)>& inConverter)
    {
        CDXString utf8Text;

        if (inString.nstyles() == 0)
        {
            // No inline styles in the text - we can do this in one go
            string convertedFragment = inConverter(inString.str(), CDXStyle());
            utf8Text.SetText(convertedFragment);
        }
        else
        {
            // Detect any unstyled 1st fragment and include it
            size_t firstStart = inString.stylestart(0);
            if (firstStart > 0)
            {
                string legacyFragment(inString.str().begin(), inString.str().begin() + firstStart);
                string convertedFragment = inConverter(legacyFragment, CDXStyle());

                CDXStyle style(0);
                utf8Text += convertedFragment;
                utf8Text.AddStyle(style);
            }

            // We need to work through our styles, encoding a run at a time
            for (size_t styleIndex=0; styleIndex<inString.nstyles(); styleIndex++)
            {
                size_t runStart = inString.stylestart(styleIndex);
                size_t runEnd = inString.styleend(styleIndex);
                CDXStyle style = inString.style(styleIndex);

                string legacyFragment(inString.str().begin() + runStart, inString.str().begin() + runEnd);
                string convertedFragment = inConverter(legacyFragment, style);

                style.startChar = utf8Text.length();
                utf8Text += convertedFragment;
                utf8Text.AddStyle(style);
            }
        }

        return utf8Text;
    }
}

XMLPutEnum<CDXTextJustification>::Value sXMLTextJustificationValuesText[] = {
	{kCDXTextJustification_Left,		"Left"},
	{kCDXTextJustification_Center,		"Center"},
	{kCDXTextJustification_Right,		"Right"},
	{kCDXTextJustification_Full,		"Full"},
	{kCDXTextJustification_Above,		"Above"},
	{kCDXTextJustification_Below,		"Below"},
	{kCDXTextJustification_Auto,		"Auto"},
	{kCDXTextJustification_BestInitial,	"Best"}
};

XMLPutEnum<CDXTextJustification> sXMLTextJustificationText(sXMLTextJustificationValuesText, sizeof sXMLTextJustificationValuesText, kCDXTextJustification_Left);
void XMLPut(XMLDataSink &sink, CDXTextJustification  v) {	sink.os << sXMLTextJustificationText.lookup(v); }

// *******************
// ** class CDXText **
// *******************
//
// Specialization of CDXObject for CDXText objects

CDXText::CDXText(CDXObjectID id_arg, const std::string &text /* = "" */)
	:	CDXGraphicObject(kCDXObj_Text, id_arg)
	,	m_lineStarts(NULL)
	,	m_rotationAngle(0)
	,	m_lineHeight(kCDXLineHeight_Variable)
	,	m_wrapWidth(0)
	,	m_justification(kCDXTextJustification_Left)
	,	m_labelAlignment(kCDXLabelDisplay_Auto)
{
    SetText(CDXString(text));
}

CDXText::CDXText (const CDXText &src)
	:	CDXGraphicObject(src)
	,	m_lineStarts((src.m_lineStarts == NULL) ? NULL : new vector<UINT16>(*src.m_lineStarts))
	,	CDX_INHERIT_SRC (m_rotationAngle)
	,	CDX_INHERIT_SRC (m_lineHeight)
	,	CDX_INHERIT_SRC (m_wrapWidth)
	,	CDX_INHERIT_SRC (m_justification)
	,	CDX_INHERIT_SRC (m_labelAlignment)
    ,	CDX_INHERIT_SRC (m_legacyText)
    ,	CDX_INHERIT_SRC (m_utf8Text)
{
}

CDXText::~CDXText()
{
	DeleteAndNull(m_lineStarts);
}

CDXObject*	CDXText::Clone() const
{
	return new CDXText (*this);
}

std::string CDXText::XMLObjectName() const
{
	return kCDXML_text;
}

void CDXText::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
        case kCDXProp_Text:
        {
            m_legacyText.Read(src_arg, size_arg);
            break;
        }
        
#if UTF8_STD_STRING
        case kCDXProp_UTF8Text:
        {
            m_utf8Text.Read(src_arg, size_arg);
            m_legacyText.Clear();
            break;
        }
#endif
        
        case kCDXProp_RotationAngle:
        {
            m_rotationAngle = src_arg.GetUINT(size_arg);
            break;
        }
        
        case kCDXProp_Justification:
        {
            m_justification = (CDXTextJustification) src_arg.GetINT(size_arg);
            break;
        }
        
        case kCDXProp_LineHeight:
        {
            m_lineHeight = src_arg.GetUINT(size_arg);
            if (size_arg == 2 && m_lineHeight != kCDXLineHeight_Automatic && m_lineHeight != kCDXLineHeight_Variable)
            {
                m_lineHeight = (m_lineHeight * 0x10000) / 20; // pre-6.0 files have this in 1440'ths of an inch
            }
            break;
        }
        
        case kCDXProp_WordWrapWidth:
        {
            m_wrapWidth = src_arg.GetUINT(size_arg);
            if (size_arg == 2)
            {
                m_wrapWidth = m_wrapWidth * 0x10000; // pre-6.0 files have this in 72's of an inch
            }
            break;
        }
        
        case kCDXProp_LabelAlignment:
        {
            m_labelAlignment = (CDXLabelDisplay) src_arg.GetUINT(size_arg);
            if (int(m_labelAlignment) > kCDXLabelDisplay_BestInitial) // some beta versions wrote the wrong codes
            {
                switch(int(m_labelAlignment))
                {
                case  7: m_labelAlignment = kCDXLabelDisplay_Left;   break;
                case  8: m_labelAlignment = kCDXLabelDisplay_Center; break;
                case  9: m_labelAlignment = kCDXLabelDisplay_Right;  break;
                case 13: m_labelAlignment = kCDXLabelDisplay_Above;  break;
                case 14: m_labelAlignment = kCDXLabelDisplay_Below;  break;
                case 15: m_labelAlignment = kCDXLabelDisplay_Auto;   break;
                }
            }
            break;
        }
        
        case kCDXProp_LineStarts:
        {
            UINT16 nLines = src_arg.GetUINT16();
            if (size_arg != size_t(2 * (nLines + 1)))
                throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
            DeleteAndNull(m_lineStarts);
            m_lineStarts = new vector<UINT16>;
            CDXReadItems(*m_lineStarts, nLines, src_arg);
            break;
        }
        
        default:
        {
            CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
            break;
        }
	}
}

void CDXText::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_RotationAngle:
	{
        m_rotationAngle = (UINT32) std::strtoul(attribValue_arg.c_str(),NULL,10);
		break;
	}
	case kCDXProp_Justification:
		m_justification = sXMLTextJustificationText.lookup(attribValue_arg);
		break;
	case kCDXProp_LineHeight:
		if (attribValue_arg == kCDXML_auto)
			m_lineHeight = kCDXLineHeight_Automatic;
		else if (attribValue_arg == kCDXML_variable)
			m_lineHeight = kCDXLineHeight_Variable;
		else
			m_lineHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		break;
	case kCDXProp_WordWrapWidth:
		m_wrapWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		break;
	case kCDXProp_LabelAlignment:
		m_labelAlignment = sXMLLabelDisplay.lookup(attribValue_arg);
		break;
	case kCDXProp_LineStarts:
	{
		std::istringstream is(attribValue_arg);
		DeleteAndNull(m_lineStarts);
		m_lineStarts = new vector<UINT16>;
		std::copy(std::istream_iterator<UINT16>(is),
				  std::istream_iterator<UINT16>(),
				  std::back_inserter(*m_lineStarts));
		break;
	}
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXText::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_rotationAngle != 0)
    {
    	sink_arg.PutAttribute( kCDXProp_RotationAngle, m_rotationAngle );
    }
    
	if (m_justification != kCDXTextJustification_Left)
    {
    	sink_arg.PutAttribute( kCDXProp_Justification, (UINT8) m_justification );
    }
    
	// Put this out in 1440'ths of an inch, for compatibility with older ChemDraw versions
	if (m_lineHeight == kCDXLineHeight_Automatic)
    {
		sink_arg.PutAttribute( kCDXProp_LineHeight, INT16(kCDXLineHeight_Automatic) );
    }
    else
    {
    	sink_arg.PutAttribute( kCDXProp_LineHeight, INT16(m_lineHeight * 20.0 / 65536.0 + 0.5) );
    }
    
	// Put this out in 72's of an inch, for compatibility with older ChemDraw versions
	if (m_wrapWidth != 0)
    {
		sink_arg.PutAttribute( kCDXProp_WordWrapWidth, INT16(m_wrapWidth / 65536.0 + 0.5) );
    }
    
	if (m_labelAlignment != kCDXLabelDisplay_Auto)
    {
		sink_arg.PutAttribute( kCDXProp_LabelAlignment, (UINT8) m_labelAlignment );
    }
    
	if (m_lineStarts != NULL && !m_lineStarts->empty())
	{
		sink_arg.PutTag( kCDXProp_LineStarts );
		sink_arg.Put( UINT16( 2 * (1 + m_lineStarts->size()) ) );
		sink_arg.Put( UINT16( m_lineStarts->size() ) );
		Put( sink_arg, *m_lineStarts );
	}
    
    if (m_legacyText.length() != 0)
    {
        sink_arg.PutAttributeCDXString( kCDXProp_Text, m_legacyText );
    }
    
#if UTF8_STD_STRING
    if (m_utf8Text.length() != 0)
    {
        sink_arg.PutAttributeCDXString( kCDXProp_UTF8Text, m_utf8Text );
    }
#endif
}

void CDXText::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_rotationAngle != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_RotationAngle, m_rotationAngle );

	if (m_justification != kCDXTextJustification_Left)
		CDXMLPutAttribute(sink_arg, kCDXProp_Justification, m_justification );

	if (m_lineHeight == kCDXLineHeight_Automatic)
		CDXMLPutAttribute(sink_arg, kCDXProp_LineHeight, std::string(kCDXML_auto) );
	else if (m_lineHeight != kCDXLineHeight_Variable)
		CDXMLPutAttribute(sink_arg, kCDXProp_LineHeight, CDXCoordinateToString(m_lineHeight) );

	if (m_wrapWidth != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_WordWrapWidth, CDXCoordinateToString(m_wrapWidth) );

	if (m_labelAlignment != kCDXLabelDisplay_Auto)
		CDXMLPutAttribute(sink_arg, kCDXProp_LabelAlignment, m_labelAlignment );

	if (m_lineStarts != NULL && !m_lineStarts->empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_LineStarts, *m_lineStarts);
}

// CDXText::XMLNeedToWriteContent
//
// Override this to indicate that we need to write content between the start-tag and end-tag.
bool CDXText::XMLNeedToWriteContent() const
{
	return true;
}

// CDXText::XMLWriteContent
//
// Override this to write content between the start-tag and end-tag.
void CDXText::XMLWriteContent(XMLDataSink &sink_arg) const
{
	XMLWriteCDXString(sink_arg, GetText());
}

/**
 *  Get our text - we return the UTF-8 representation in UTF-8 mode, otherwise the legacy representation
 *
 *  @return Our text
 */
const CDXString& CDXText::GetText() const
{
#if UTF8_STD_STRING
    return m_utf8Text;
#else
    return m_legacyText;
#endif
}

/**
 *  Set our text - we set the UTF-8 representation in UTF-8 mode, otherwise the legacy representation
 *
 *  @param inString Text to set
 */
void CDXText::SetText(const CDXString& inString)
{
#if UTF8_STD_STRING
    m_utf8Text = inString;
#else
    m_legacyText = inString;
#endif
}

/**
 *  Append to our text - we append to the UTF-8 representation in UTF-8 mode, otherwise the legacy representation
 *
 *  @param inString Text to append
 */
void CDXText::AppendText(const CDXString& inString)
{
#if UTF8_STD_STRING
    m_utf8Text += inString;
#else
    m_legacyText += inString;
#endif
}

/**
 *  Copy our legacy string representation to our UTF-8 representation if (and only if) we don't already have
 *  a UTF-8 representation. Our UTF-8 representation must be in internal form.
 *
 *  @param inFontTable The font table to use
 */
void CDXText::FinalizeTextAfterRead(const CDXFontTable& inFontTable)
{
#if UTF8_STD_STRING
    CDXString convertedUTF8;

    if (m_utf8Text.empty())
    {
        // Get the encoding for our default font
        CDXCharSet defaultEncoding = inFontTable.EncodingForFontID(m_captionStyle.family);
        if (defaultEncoding == kCDXCharSetUnknown)
        {
            defaultEncoding = inFontTable.EncodingForFontID(m_labelStyle.family);
        }

        // Convert from a legacy encoding to UTF-8
        convertedUTF8 = ConvertTextEncodingWithConverter(m_legacyText, [defaultEncoding, inFontTable](const string& inString, const CDXStyle& inStyle) {
            CDXCharSet encoding = defaultEncoding;
            if (inStyle.family != UINT16(-1))
            {
                encoding = inFontTable.EncodingForFontID(inStyle.family);
            }

            return TranscodeString(inString, encoding, kCDXCharSetUTF8);
        });
    }
    else
    {
        // Historically we may have stored "internal" UTF-8 in m_utf8Text and written this out to the CDX file. We check for this
        // here, and convert any legacy characters to real UTF-8.
        convertedUTF8 = ConvertTextEncodingWithConverter(m_utf8Text, [](const string& inString, const CDXStyle&) {
            return MakeStringSafe(inString);
        });
    }

    if (convertedUTF8.str() != m_utf8Text.str())
    {
        // Our line starts are no longer valid, as the strings have changed. We could try and recalculate these, however
        // in practice client applications that care about layout (e.g. ChemDraw) will reflow the text anyway.
        DeleteAndNull(m_lineStarts);
    }

    m_utf8Text = convertedUTF8;
#endif
}

/**
 *  Copy our UTF-8 string representation to our legacy representation, using the "best" encoding
 *  we have for each font. This is less than perfect, and is for backwards compatibility only.
 *
 *  @param inFontTable The font table to update
 */
void CDXText::PrepareTextForWrite(const CDXFontTable& inFontTable)
{
#if UTF8_STD_STRING
    m_legacyText.Clear();
    
    // Get the encoding for our default font
    CDXCharSet defaultEncoding = inFontTable.EncodingForFontID(m_captionStyle.family);
    if (defaultEncoding == kCDXCharSetUnknown)
    {
        defaultEncoding = inFontTable.EncodingForFontID(m_labelStyle.family);
    }
    
    if (m_utf8Text.nstyles() == 0)
    {
        // No inline styles in the text - use our default font information
        string legacyStr = TranscodeString(m_utf8Text.str(), kCDXCharSetUTF8, defaultEncoding);
        m_legacyText.SetText(legacyStr);
    }
    else
    {
        ASSERT(m_utf8Text.stylestart(0) == 0);

        // We need to work through our styles, encoding a run at a time
        for (size_t styleIndex=0; styleIndex<m_utf8Text.nstyles(); styleIndex++)
        {
            size_t runStart = m_utf8Text.stylestart(styleIndex);
            size_t runEnd = m_utf8Text.styleend(styleIndex);
            CDXStyle style = m_utf8Text.style(styleIndex);
            
            CDXCharSet encoding = defaultEncoding;
            if (style.family != UINT16(-1))
            {
                encoding = inFontTable.EncodingForFontID(style.family);
            }
            
            string utf8Fragment(m_utf8Text.str().begin() + runStart, m_utf8Text.str().begin() + runEnd);
            string legacyFragment = TranscodeString(utf8Fragment, kCDXCharSetUTF8, encoding);

            // It's possible that we couldn't transcode the string at all
            if (legacyFragment.length() > 0)
            {
                style.startChar = m_legacyText.length();
                m_legacyText += legacyFragment;
                m_legacyText.AddStyle(style);
            }
        }
    }
#endif
}


static std::list<const CDXFontTable *> sFontTables;
void CDXPushFontTable(const CDXFontTable *fontTable)
{
	sFontTables.push_back(fontTable);
}

void CDXPopFontTable(const CDXFontTable *fontTable)
{
	if (find(sFontTables.begin(), sFontTables.end(), fontTable) != sFontTables.end())
		sFontTables.erase(find(sFontTables.begin(), sFontTables.end(), fontTable));
}

static CDXCharSet Font_IDToCharSet(int family)
{
	if (!sFontTables.empty())
	{
		try
		{
			return sFontTables.back()->LookupFont(family).m_type;
		}
		catch (...)
		{
		}
	}
	return kCDXCharSetUnknown;
}

static std::string Font_IDToName(int family)
{
	std::string	result;
	if (!sFontTables.empty())
	{
		try
		{
			result = sFontTables.back()->LookupFont(family).m_family;
		}
		catch (...)
		{
		}
	}
	return result;
}

static void XMLWriteCDXString(XMLDataSink &sink_arg, const CDXString &theText)
{
	const std::vector<CDXStyle> &v = theText.styles();
	if (v.empty())
	{
		sink_arg.os << "<" << kCDXML_string << ">";
		XMLPut(sink_arg, theText.str());
		sink_arg.os << "</" << kCDXML_string << ">";
	}
	else
	{
		if (v[0].startChar > 0)
		{
			sink_arg.os << "<" << kCDXML_string << ">";
			XMLPut(sink_arg, theText.str().substr(0, v[0].startChar));
			sink_arg.os << "</" << kCDXML_string << ">";
		}
		for (std::vector<CDXStyle>::const_iterator i = v.begin();  i != v.end();  ++i)
		{
			sink_arg.os << "<" << kCDXML_string << " "
				<< kCDXML_font << "=\"" << (i->family == 0xFFFF ? -1 : i->family) << "\" "
				<< kCDXML_size << "=\"" << i->size / 20.0 << "\"";

			if (i->color != 3)
			{
				sink_arg.os << " " << kCDXML_color << "=\"" << i->color << "\"";
                if (i->has_alpha)
                {
                    sink_arg.os << " " << kCDXML_alpha << "=\"" << i->alpha / 65535 << "\"";
                }
			}

            if (i->face != 0)
            {
                sink_arg.os << " " << kCDXML_face << "=\"" << i->face << "\"";
            }

			sink_arg.os << ">";

            size_t nChars = 0;
            if ((i + 1) == v.end())
            {
                nChars = theText.str().size();
            }
            else
            {
                nChars = ((i + 1)->startChar) - i->startChar;
            }

            if ((nChars > 0) && (i->startChar < theText.str().size()))
            {
                string utf8 = StringToUnicode(Font_IDToName(i->family), theText.str().substr(i->startChar, nChars), Font_IDToCharSet(i->family));
                XMLPut(sink_arg, utf8);
            }

			sink_arg.os << "</" << kCDXML_string << ">";

		}
	}
}
