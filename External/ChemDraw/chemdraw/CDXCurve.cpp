// CommonCS/LibCommon/Src/CDXCurve.cpp
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
#include "CDXMLNames.h"
#include "XMLPutEnum.h"
#include <stdlib.h> // for atoi
#include <sstream>
#include <stdexcept>

extern XMLPutEnum<CDXArrowHeadType> sXMLArrowHeadType;
extern XMLPutEnum<CDXArrowHeadPosition> sXMLArrowHeadPositionType;

// ********************
// ** class CDXCurve **
// ********************
//
// Specialization of CDXObject for CDXCurve objects

CDXCurve::CDXCurve(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_Curve, id_arg)
 , m_curveType(0)
 , m_arrowHeadType(kCDXArrowHeadType_Unspecified)
 , m_arrowHeadAtEnd(kCDXArrowHeadPosition_Unspecified)
 , m_arrowHeadAtStart(kCDXArrowHeadPosition_Unspecified)
 , m_fadePercent(1000)
 , m_headSize(0)
 , m_headCenterSize(0)
 , m_headWidth(0)
 , m_spacing(0)
 , m_closed(false)
{
}

CDXCurve::CDXCurve (const CDXCurve &src)
	:	CDXGraphicObject(src)
	,	m_curveType				(src.m_curveType)
	,	m_arrowHeadType			(src.m_arrowHeadType)
	,	m_arrowHeadAtEnd		(src.m_arrowHeadAtEnd)
	,	m_arrowHeadAtStart		(src.m_arrowHeadAtStart)
	,	m_fadePercent			(src.m_fadePercent)
	,	m_headSize				(src.m_headSize)
	,	m_headCenterSize		(src.m_headCenterSize)
	,	m_headWidth				(src.m_headWidth)
	,	m_spacing				(src.m_spacing)
	,	m_closed				(src.m_closed)
	,	m_curvePoints			(src.m_curvePoints)
	,	m_curvePoints3D			(src.m_curvePoints3D)
{
}

CDXCurve::~CDXCurve()
{
}

CDXObject*	CDXCurve::Clone() const
{
	return new CDXCurve (*this);
}

std::string CDXCurve::XMLObjectName() const
{
	return kCDXML_curve;
}

void CDXCurve::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Curve_Type:
		m_curveType = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Type:
		m_arrowHeadType = (CDXArrowHeadType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Head:
		m_arrowHeadAtEnd = (CDXArrowHeadPosition) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Tail:
		m_arrowHeadAtStart = (CDXArrowHeadPosition) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Size:
		m_headSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_CenterSize:
		m_headCenterSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Width:
		m_headWidth = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Curve_Spacing:
		m_spacing = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Closed:
		m_closed = true;
		break;
	case kCDXProp_Curve_Points:
		{
		UINT16 nPoints = src_arg.GetUINT16();
		if (size_arg != 2 + nPoints * 8)
			throw std::runtime_error("Incorrect size for Curve_Points");
		CDXReadItems(m_curvePoints, nPoints, src_arg);
		break;
		}
	case kCDXProp_Curve_Points3D:
		{
		UINT16 nPoints = src_arg.GetUINT16();
		if (size_arg != 2 + nPoints * 12)
			throw std::runtime_error("Incorrect size for Curve_Points3D");
		CDXReadItems(m_curvePoints3D, nPoints, src_arg);
		break;
		}
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXCurve::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Curve_Type:
		m_curveType = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrowhead_Type:
		m_arrowHeadType = sXMLArrowHeadType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Head:
		m_arrowHeadAtEnd = sXMLArrowHeadPositionType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Tail:
		m_arrowHeadAtStart = sXMLArrowHeadPositionType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Size:
		m_headSize = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrowhead_CenterSize:
		m_headCenterSize = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrowhead_Width:
		m_headWidth = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Curve_Spacing:
		m_spacing = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Closed:
		m_closed = attribValue_arg == "yes" || attribValue_arg == "true";
		break;
	case kCDXProp_Curve_Points:
	{
		std::istringstream is(attribValue_arg);
		double x, y;
		while (!is.eof())
		{
			is >> x >> y;
			m_curvePoints.push_back(CDXPoint2D(CDXCoordinatefromPoints(x), CDXCoordinatefromPoints(y)));
		}
		break;
	}
	case kCDXProp_Curve_Points3D:
	{
		std::istringstream is(attribValue_arg);
		double x, y, z;
		while (!is.eof())
		{
			is >> x >> y >> z;
			m_curvePoints3D.push_back(CDXPoint3D(CDXCoordinatefromPoints(x), CDXCoordinatefromPoints(y), CDXCoordinatefromPoints(z)));
		}
		break;
	}
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXCurve::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_curveType != 0)
		sink_arg.PutAttribute( kCDXProp_Curve_Type, m_curveType );

	if (m_arrowHeadType != kCDXArrowHeadType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Type, (INT16) m_arrowHeadType );

	if (m_arrowHeadAtEnd != kCDXArrowHeadPosition_Unspecified && m_arrowHeadAtEnd != kCDXArrowHeadPosition_None)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Head, (INT16) m_arrowHeadAtEnd );

	if (m_arrowHeadAtStart != kCDXArrowHeadPosition_Unspecified && m_arrowHeadAtStart != kCDXArrowHeadPosition_None)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Tail, (INT16) m_arrowHeadAtStart );

	if (m_headSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Size, m_headSize );

	if (m_headCenterSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_CenterSize, m_headCenterSize );

	if (m_headWidth != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Width, m_headWidth );

	if (m_spacing != 0)
		sink_arg.PutAttribute( kCDXProp_Curve_Spacing, m_spacing );

	if (m_fadePercent != 1000)
		sink_arg.PutAttribute( kCDXProp_FadePercent, m_fadePercent );

	if (m_closed)
		sink_arg.PutAttribute(kCDXProp_Closed);

	if (!m_curvePoints.empty())
	{
		sink_arg.PutTag( kCDXProp_Curve_Points );
		sink_arg.Put( UINT16( 2 + 8 * m_curvePoints.size() ) );
		sink_arg.Put( UINT16( m_curvePoints.size() ) );
		Put( sink_arg, m_curvePoints );
	}

	if (!m_curvePoints3D.empty())
	{
		sink_arg.PutTag( kCDXProp_Curve_Points3D );
		sink_arg.Put( UINT16( 2 + 12 * m_curvePoints3D.size() ) );
		sink_arg.Put( UINT16( m_curvePoints3D.size() ) );
		Put( sink_arg, m_curvePoints3D );
	}
}

void CDXCurve::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_curveType != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Type, m_curveType );

	if (m_arrowHeadAtEnd != kCDXArrowHeadPosition_Unspecified && m_arrowHeadAtEnd != kCDXArrowHeadPosition_None)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Head, m_arrowHeadAtEnd );

	if (m_arrowHeadAtStart != kCDXArrowHeadPosition_Unspecified && m_arrowHeadAtStart != kCDXArrowHeadPosition_None)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Tail, m_arrowHeadAtStart );

	if (m_arrowHeadType != kCDXArrowHeadType_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Type, m_arrowHeadType );

	if (m_headSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Size, m_headSize );

	if (m_headCenterSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_CenterSize, m_headCenterSize );

	if (m_headWidth != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Width, m_headWidth );

	if (m_spacing != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Spacing, m_spacing );

	if (m_fadePercent != 1000)
		CDXMLPutAttribute(sink_arg, kCDXProp_FadePercent, m_fadePercent );

	if (m_closed)
		CDXMLPutAttribute(sink_arg, kCDXProp_Closed);

	if (!m_curvePoints.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Points, m_curvePoints );

	if (!m_curvePoints3D.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Points3D, m_curvePoints3D );

}
