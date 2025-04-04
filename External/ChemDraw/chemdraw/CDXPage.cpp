// CommonCS/LibCommon/Src/CDXPage.cpp
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
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <sstream>

XMLPutEnum<CDXPageDefinition>::Value sXMLPageDefinitionValuesPage[] = {
	{kCDXPageDefinition_Undefined,			"Undefined"},
	{kCDXPageDefinition_Center,				"Center"},
	{kCDXPageDefinition_TL4,				"TL4"},
	{kCDXPageDefinition_IDTerm,				"IDTerm"},
	{kCDXPageDefinition_FlushLeft,			"FlushLeft"},
	{kCDXPageDefinition_FlushRight,			"FlushRight"},
	{kCDXPageDefinition_Reaction1,			"Reaction1"},
	{kCDXPageDefinition_Reaction2,			"Reaction2"},
	{kCDXPageDefinition_MulticolumnTL4,		"MulticolumnTL4"},
	{kCDXPageDefinition_MulticolumnNonTL4,	"MulticolumnNonTL4"},
	{kCDXPageDefinition_UserDefined,		"UserDefined"}
};

XMLPutEnum<CDXPageDefinition> sXMLPageDefinitionPage(sXMLPageDefinitionValuesPage, sizeof sXMLPageDefinitionValuesPage, kCDXPageDefinition_Undefined);
void XMLPut(XMLDataSink &sink, CDXPageDefinition  v) {	sink.os << sXMLPageDefinitionPage.lookup(v); }

// *******************
// ** class CDXPage **
// *******************
//
// Specialization of CDXObject for CDXPage objects

const UINT32 CDXPage::has_boundingBox    = 0x00000002;
const UINT32 CDXPage::has_pageWidth      = 0x00000004;
const UINT32 CDXPage::has_pageHeight     = 0x00000008;
const UINT32 CDXPage::has_headerPosition = 0x00000010;
const UINT32 CDXPage::has_footerPosition = 0x00000020;
const UINT32 CDXPage::has_pageOverlap    = 0x00000040;
const UINT32 CDXPage::has_printTrimMarks = 0x00000080;
const UINT32 CDXPage::has_heightPages    = 0x00000100;
const UINT32 CDXPage::has_widthPages     = 0x00000200;
const UINT32 CDXPage::has_bgColor		 = 0x00000400;
const UINT32 CDXPage::has_boundsInParent = 0x00000800;

CDXPage::CDXPage(CDXObjectID id_arg)
	:	CDXObject			(kCDXObj_Page, id_arg)
	,	m_bgColor			(0)
	,	m_pageWidth			(0)
	,	m_pageHeight		(0)
	,	m_headerPosition	(0)
	,	m_footerPosition	(0)
	,	m_pageOverlap		(0)
	,	m_printTrimMarks	(false)
	,   m_drawingSpaceType	(kCDXDrawingSpace_Pages)
	,	m_heightPages		(0)
	,	m_widthPages		(0)
	,	m_pageDefinition	(kCDXPageDefinition_Undefined)
	,   m_flags(0)
{
}

CDXPage::~CDXPage()
{
}

CDXPage::CDXPage (const CDXPage &src)
	:	CDXObject			(src)
	,	m_boundingBox		(src.m_boundingBox)
	,	m_bgColor			(src.m_bgColor)
	,	m_pageWidth			(src.m_pageWidth)
	,	m_pageHeight		(src.m_pageHeight)
	,	m_headerPosition	(src.m_headerPosition)
	,	m_footerPosition	(src.m_footerPosition)
	,	m_pageOverlap		(src.m_pageOverlap)
	,	m_header			(src.m_header)
	,	m_footer			(src.m_footer)
	,	m_printTrimMarks	(src.m_printTrimMarks)
	,	m_drawingSpaceType	(src.m_drawingSpaceType)
	,	m_heightPages		(src.m_heightPages)
	,	m_widthPages		(src.m_widthPages)
	,	m_splitterPositions	(src.m_splitterPositions)
	,	m_pageDefinition	(src.m_pageDefinition)
	,	m_boundsInParent	(src.m_boundsInParent)
	,	m_flags				(src.m_flags)
{
}

CDXObject*	CDXPage::Clone() const
{
	return new CDXPage (*this);
}

std::string CDXPage::XMLObjectName() const
{
	return kCDXML_page;
}

void CDXPage::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_BoundingBox:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetINT32();
		CDXCoordinate left   = src_arg.GetINT32();
		CDXCoordinate bottom = src_arg.GetINT32();
		CDXCoordinate right  = src_arg.GetINT32();
		BoundingBox(CDXRectangle(top, left, bottom, right));
		break;
		}
	case kCDXProp_BackgroundColor:
		BGColor(UINT16(src_arg.GetUINT(size_arg)));
		break;
	case kCDXProp_WidthPages:
		m_widthPages = src_arg.GetUINT(size_arg);
		m_flags |= has_widthPages;
		break;
	case kCDXProp_HeightPages:
		m_heightPages = src_arg.GetUINT(size_arg);
		m_flags |= has_heightPages;
		break;
	case kCDXProp_DrawingSpaceType:
		m_drawingSpaceType = (CDXDrawingSpaceType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Width:
		m_pageWidth = src_arg.GetUINT(size_arg);
		m_flags |= has_pageWidth;
		break;
	case kCDXProp_Height:
		m_pageHeight = src_arg.GetUINT(size_arg);
		m_flags |= has_pageHeight;
		break;
	case kCDXProp_PageOverlap:
		m_pageOverlap = src_arg.GetUINT(size_arg);
		m_flags |= has_pageOverlap;
		break;
	case kCDXProp_Header:
		m_header.Read(src_arg, size_arg);
		break;
	case kCDXProp_HeaderPosition:
		m_headerPosition = src_arg.GetUINT(size_arg);
		m_flags |= has_headerPosition;
		break;
	case kCDXProp_Footer:
		m_footer.Read(src_arg, size_arg);
		break;
	case kCDXProp_FooterPosition:
		m_footerPosition = src_arg.GetUINT(size_arg);
		m_flags |= has_footerPosition;
		break;
	case kCDXProp_PrintTrimMarks:
		m_printTrimMarks = size_arg == 0 || src_arg.GetUINT(size_arg) != 0;
		m_flags |= has_printTrimMarks;
		break;
	case kCDXProp_SplitterPositions:
		CDXReadItems(m_splitterPositions, size_arg/sizeof(INT32), src_arg);
		break;
	case kCDXProp_PageDefinition:
		m_pageDefinition=(CDXPageDefinition)src_arg.GetUINT8();
		if (size_arg > 1)
			src_arg.GetString(size_arg - 1);	// Harold found a bogus CDX file from somewhere that had something strange stored in this field
		break;
	case kCDXProp_BoundsInParent:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetINT32();
		CDXCoordinate left   = src_arg.GetINT32();
		CDXCoordinate bottom = src_arg.GetINT32();
		CDXCoordinate right  = src_arg.GetINT32();
		BoundsInParent(CDXRectangle(top, left, bottom, right));
		break;
		}
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPage::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_BoundingBox:
		BoundingBox(StringToCDXRectangle(attribValue_arg));
		break;
	case kCDXProp_BackgroundColor:
		BGColor( atoi(attribValue_arg.c_str()) );
		break;
	case kCDXProp_WidthPages:
		m_widthPages = atoi(attribValue_arg.c_str());
		m_flags |= has_widthPages;
		break;
	case kCDXProp_HeightPages:
		m_heightPages = atoi(attribValue_arg.c_str());
		m_flags |= has_heightPages;
		break;
	case kCDXProp_DrawingSpaceType:
		m_drawingSpaceType = attribValue_arg == kCDXML_poster ? kCDXDrawingSpace_Poster : kCDXDrawingSpace_Pages;
		break;
	case kCDXProp_Width:
		m_pageWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_pageWidth;
		break;
	case kCDXProp_Height:
		m_pageHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_pageHeight;
		break;
	case kCDXProp_PageOverlap:
		m_pageOverlap = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_pageOverlap;
		break;
	case kCDXProp_Header:
		m_header = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_HeaderPosition:
		m_headerPosition = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_headerPosition;
		break;
	case kCDXProp_Footer:
		m_footer = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_FooterPosition:
		m_footerPosition = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_footerPosition;
		break;
	case kCDXProp_PrintTrimMarks:
		m_printTrimMarks = attribValue_arg == "yes" || attribValue_arg == "true";
		m_flags |= has_printTrimMarks;
		break;
	case kCDXProp_SplitterPositions:
	{
		std::istringstream is(attribValue_arg);
		m_splitterPositions.clear();
		while (!is.eof())
		{
			double pos;
			is >> pos;
			m_splitterPositions.push_back(CDXCoordinatefromPoints(pos));
		}
		break;
	}
	case kCDXProp_PageDefinition:
		m_pageDefinition = sXMLPageDefinitionPage.lookup(attribValue_arg);
		break;
	case kCDXProp_BoundsInParent:
		BoundsInParent(StringToCDXRectangle(attribValue_arg));
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPage::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	if (Known(has_boundingBox))
		sink_arg.PutAttribute( kCDXProp_BoundingBox, BoundingBox() );

	if (KnownBGColor())
		sink_arg.PutAttribute( kCDXProp_BackgroundColor, BGColor() );

	if (Known(has_pageWidth))
		sink_arg.PutAttribute( kCDXProp_Width, m_pageWidth );

	if (Known(has_pageHeight))
		sink_arg.PutAttribute( kCDXProp_Height, m_pageHeight );

	if (Known(has_headerPosition))
		sink_arg.PutAttribute( kCDXProp_HeaderPosition, m_headerPosition );

	if (Known(has_footerPosition))
		sink_arg.PutAttribute( kCDXProp_FooterPosition, m_footerPosition );

	if (Known(has_pageOverlap))
		sink_arg.PutAttribute( kCDXProp_PageOverlap, m_pageOverlap );

	if (Known(has_printTrimMarks) && m_printTrimMarks)
		sink_arg.PutAttribute( kCDXProp_PrintTrimMarks );

	if (Known(has_heightPages))
		sink_arg.PutAttribute( kCDXProp_HeightPages, m_heightPages );

	if (Known(has_widthPages))
		sink_arg.PutAttribute( kCDXProp_WidthPages, m_widthPages );

	if (m_header.length())
		sink_arg.PutAttributeCDXString( kCDXProp_Header, m_header );

	if (m_footer.length())
		sink_arg.PutAttributeCDXString( kCDXProp_Footer, m_footer );

	if (m_drawingSpaceType != kCDXDrawingSpace_Pages)
		sink_arg.PutAttribute( kCDXProp_DrawingSpaceType, INT8(m_drawingSpaceType) );
	
	if (!m_splitterPositions.empty())
	{
		sink_arg.PutTag( kCDXProp_SplitterPositions );
		sink_arg.Put( UINT16( 4 * m_splitterPositions.size() ) );
		Put( sink_arg, m_splitterPositions );
	}

	if (m_pageDefinition != kCDXPageDefinition_Undefined)
		sink_arg.PutAttribute( kCDXProp_PageDefinition, UINT8(m_pageDefinition));

	if (Known(has_boundsInParent))
		sink_arg.PutAttribute( kCDXProp_BoundsInParent, BoundsInParent() );

}

void CDXPage::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	if (Known(has_boundingBox))
		CDXMLPutAttribute(sink_arg, kCDXProp_BoundingBox, BoundingBox() );

	if (KnownBGColor())
		CDXMLPutAttribute(sink_arg, kCDXProp_BackgroundColor, BGColor() );

	if (Known(has_pageWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_Width, CDXCoordinateToString(m_pageWidth) );

	if (Known(has_pageHeight))
		CDXMLPutAttribute(sink_arg, kCDXProp_Height, CDXCoordinateToString(m_pageHeight) );

	if (Known(has_headerPosition))
		CDXMLPutAttribute(sink_arg, kCDXProp_HeaderPosition, CDXCoordinateToString(m_headerPosition) );

	if (Known(has_footerPosition))
		CDXMLPutAttribute(sink_arg, kCDXProp_FooterPosition, CDXCoordinateToString(m_footerPosition) );

	if (Known(has_pageOverlap))
		CDXMLPutAttribute(sink_arg, kCDXProp_PageOverlap, CDXCoordinateToString(m_pageOverlap) );

	if (Known(has_printTrimMarks) && m_printTrimMarks)
		CDXMLPutAttribute(sink_arg, kCDXProp_PrintTrimMarks );

	if (Known(has_heightPages))
		CDXMLPutAttribute(sink_arg, kCDXProp_HeightPages, m_heightPages );

	if (Known(has_widthPages))
		CDXMLPutAttribute(sink_arg, kCDXProp_WidthPages, m_widthPages );

	if (m_header.length())
		CDXMLPutAttribute(sink_arg, kCDXProp_Header, m_header.str() );

	if (m_footer.length())
		CDXMLPutAttribute(sink_arg, kCDXProp_Footer, m_footer.str() );

	if (m_drawingSpaceType != kCDXDrawingSpace_Pages)
		CDXMLPutAttribute(sink_arg, kCDXProp_DrawingSpaceType, std::string(kCDXML_poster) );

	if (!m_splitterPositions.empty())
	{
		sink_arg.os << " " << CDXMLAttributeName(kCDXProp_SplitterPositions) << "=\"";
		for (std::vector<CDXCoordinate>::const_iterator i = m_splitterPositions.begin();  i != m_splitterPositions.end();  ++i)
		{
			if (i != m_splitterPositions.begin()) sink_arg.os << " ";
			XMLPut(sink_arg, CDXCoordinateToString(*i));
		}
		sink_arg.os << "\"" << GetTextEOL();
	}
	if (m_pageDefinition != kCDXPageDefinition_Undefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_PageDefinition, m_pageDefinition);

	if (Known(has_boundsInParent))
		CDXMLPutAttribute(sink_arg, kCDXProp_BoundsInParent, BoundsInParent() );

}
