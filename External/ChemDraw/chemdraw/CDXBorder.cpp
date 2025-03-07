// CommonCS/LibCommon/Src/CDXBorder.cpp
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

#include <ostream>

#include "CDXStdObjects.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"

extern XMLPutEnum<CDXLineType> sXMLLineType;

XMLPutEnum<CDXSideType>::Value sXMLBorderSideTypeValues[] = {
	{kCDXSideType_Undefined,	"undefined"},
	{kCDXSideType_Top,			"top"},
	{kCDXSideType_Left,			"left"},
	{kCDXSideType_Bottom,		"bottom"},
	{kCDXSideType_Right,		"right"}
};

XMLPutEnum<CDXSideType> sXMLBorderSideType(sXMLBorderSideTypeValues, sizeof sXMLBorderSideTypeValues, kCDXSideType_Undefined);
void XMLPut(XMLDataSink &sink, CDXSideType  v) {	sink.os << sXMLBorderSideType.lookup(v); }

// ***********************
// **  class CDXBorder  **
// ***********************
//
// Specialization of CDXObject for CDXBorder objects

CDXBorder::CDXBorder(CDXObjectID id) 
: 	CDXObject(kCDXObj_Border, id), 
	m_side(kCDXSideType_Undefined), 
	m_lineType(CDXLineType(0)), 
	m_lineWidth(1), has_lineWidth(false), 
	m_color(0), has_color(false), m_alpha(0), has_alpha(false)
{
}

CDXBorder::CDXBorder (const CDXBorder &src) 
:	CDXObject(src), 
	m_side(src.m_side), 
	m_lineType(src.m_lineType), 
	m_lineWidth(src.m_lineWidth), has_lineWidth(src.has_lineWidth),
	m_color(src.m_color), m_alpha(src.m_alpha), has_color(src.has_color), has_alpha(src.has_alpha)
{
}


CDXObject* CDXBorder::Clone() const
{
	return new CDXBorder(*this);

}


std::string CDXBorder::XMLObjectName() const
{
	return kCDXML_border;
}

void CDXBorder::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_LineWidth:
		m_lineWidth = src_arg.GetUINT(size_arg);
		has_lineWidth = true;
		break;

	case kCDXProp_Line_Type:
		m_lineType = (CDXLineType) src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_ForegroundColor:
		m_color = UINT16(src_arg.GetUINT(size_arg));
		has_color = true;
		break;

	case kCDXProp_ForegroundAlpha:
		m_alpha = UINT16(src_arg.GetUINT(size_arg));
		has_alpha = true;
		break;

	case kCDXProp_Side:
		m_side=(CDXSideType)src_arg.GetUINT(size_arg);
		break;

	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}

}

void CDXBorder::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Line_Type:
		m_lineType = sXMLLineType.lookup(attribValue_arg);
		break;

	case kCDXProp_LineWidth:
		m_lineWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		has_lineWidth = true;
		break;

	case kCDXProp_ForegroundColor:
		m_color = atoi(attribValue_arg.c_str());
		has_color = true;
		break;

	case kCDXProp_ForegroundAlpha:
		m_alpha = atoi(attribValue_arg.c_str()) * 65535;
		has_alpha = true;
		break;

	case kCDXProp_Side:
		m_side = sXMLBorderSideType.lookup(attribValue_arg);
		break;

	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}


void CDXBorder::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	if (m_lineType != kCDXLineType_Solid)
		sink_arg.PutAttribute( kCDXProp_Line_Type, (INT8) m_lineType );

	if (has_lineWidth)
		sink_arg.PutAttribute( kCDXProp_LineWidth, m_lineWidth );
	
	if (has_color)
		sink_arg.PutAttribute(kCDXProp_ForegroundColor, (INT8) m_color);

	if (has_alpha)
		sink_arg.PutAttribute(kCDXProp_ForegroundAlpha, m_alpha);

	if (m_side != kCDXSideType_Undefined)
		sink_arg.PutAttribute(kCDXProp_Side, (INT16) m_side);
}


void CDXBorder::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	if (m_lineType != kCDXLineType_Solid)
		CDXMLPutAttribute(sink_arg, kCDXProp_Line_Type, m_lineType );

	if (has_lineWidth)
		CDXMLPutAttribute(sink_arg, kCDXProp_LineWidth, CDXCoordinateToString(m_lineWidth) );

	if (has_color)
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundColor, m_color);

	if (has_alpha)
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundAlpha, CDXValueToString(m_alpha / 65535.0, 4));

	if (m_side != kCDXSideType_Undefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_Side, m_side);
}



