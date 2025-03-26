// CommonCS/LibCommon/Src/CDXArrow.cpp
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

XMLPutEnum<CDXArrowHeadType>::Value sXMLArrowHeadTypeValues[] = {
	{kCDXArrowHeadType_Unspecified,	"Unspecified"},
	{kCDXArrowHeadType_Solid,		"Solid"},
	{kCDXArrowHeadType_Hollow,		"Hollow"},
	{kCDXArrowHeadType_Angle,		"Angle"}
};

XMLPutEnum<CDXArrowHeadType> sXMLArrowHeadType(sXMLArrowHeadTypeValues, sizeof sXMLArrowHeadTypeValues, kCDXArrowHeadType_Solid);
void XMLPut(XMLDataSink &sink, CDXArrowHeadType  v) {	sink.os << sXMLArrowHeadType.lookup(v); }

XMLPutEnum<CDXArrowHeadPosition>::Value sXMLArrowHeadPositionTypeValues[] = {
	{kCDXArrowHeadPosition_Unspecified,	"Unspecified"},
	{kCDXArrowHeadPosition_None,		"None"},
	{kCDXArrowHeadPosition_Full,		"Full"},
	{kCDXArrowHeadPosition_HalfLeft,	"HalfLeft"},
	{kCDXArrowHeadPosition_HalfRight,	"HalfRight"}
};

XMLPutEnum<CDXArrowHeadPosition> sXMLArrowHeadPositionType(sXMLArrowHeadPositionTypeValues, sizeof sXMLArrowHeadPositionTypeValues, kCDXArrowHeadPosition_Unspecified);
void XMLPut(XMLDataSink &sink, CDXArrowHeadPosition  v) {	sink.os << sXMLArrowHeadPositionType.lookup(v); }

XMLPutEnum<CDXNoGoType>::Value sXMLNoGoTypeValues[] = {
	{kCDXNoGoType_Unspecified,	"Unspecified"},
	{kCDXNoGoType_None,			"None"},
	{kCDXNoGoType_Cross,		"Cross"},
	{kCDXNoGoType_Hash,			"Hash"}
};

XMLPutEnum<CDXNoGoType> sXMLNoGoType(sXMLNoGoTypeValues, sizeof sXMLNoGoTypeValues, kCDXNoGoType_Unspecified);
void XMLPut(XMLDataSink &sink, CDXNoGoType  v) {	sink.os << sXMLNoGoType.lookup(v); }

// ********************
// ** class CDXArrow **
// ********************
//
// Specialization of CDXObject for CDXArrow objects

CDXArrow::CDXArrow(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_Arrow, id_arg),
	m_arrowHeadType(kCDXArrowHeadType_Unspecified),
	m_arrowHeadHead(kCDXArrowHeadPosition_None),
	m_arrowHeadTail(kCDXArrowHeadPosition_None),
	m_headSize(0),
	m_headCenterSize(0),
	m_headWidth(0),
	m_angularSize(0),
	m_shaftSpacing(0),
	m_equilibriumRatio(100),
	m_fadePercent(1000),
	m_sourceID(-1),
	m_targetID(-1),
	m_nogo(kCDXNoGoType_Unspecified),
	m_dipole(false),
	m_has3dHead(false),
	m_has3dTail(false),
	m_hasCenter(false),
	m_hasMajorAxisEnd(false),
	m_hasMinorAxisEnd(false)
{
}

CDXArrow::CDXArrow (const CDXArrow &src)
	:	CDXGraphicObject		(src)
	,	m_arrowHeadType			(src.m_arrowHeadType)
	,	m_arrowHeadHead			(src.m_arrowHeadHead)
	,	m_arrowHeadTail			(src.m_arrowHeadTail)
	,	m_headSize				(src.m_headSize)
	,	m_headCenterSize		(src.m_headCenterSize)
	,	m_headWidth				(src.m_headWidth)
	,	m_angularSize			(src.m_angularSize)
	,	m_shaftSpacing			(src.m_shaftSpacing)
	,	m_equilibriumRatio		(src.m_equilibriumRatio)
	,	m_fadePercent			(src.m_fadePercent)
	,	m_sourceID				(src.m_sourceID)
	,	m_targetID				(src.m_targetID)
	,	m_nogo					(src.m_nogo)
	,	m_dipole				(src.m_dipole)
	,	m_3dHead				(src.m_3dHead)
	,	m_has3dHead				(src.m_has3dHead)
	,	m_3dTail				(src.m_3dTail)
	,	m_has3dTail				(src.m_has3dTail)
	,	m_center				(src.m_center)
	,	m_hasCenter				(src.m_hasCenter)
	,	m_majorAxisEnd			(src.m_majorAxisEnd)
	,	m_hasMajorAxisEnd		(src.m_hasMajorAxisEnd)
	,	m_minorAxisEnd			(src.m_minorAxisEnd)
	,	m_hasMinorAxisEnd		(src.m_hasMinorAxisEnd)
{
}

CDXArrow::~CDXArrow()
{
}

CDXObject*	CDXArrow::Clone() const
{
	return new CDXArrow (*this);
}

std::string CDXArrow::XMLObjectName() const
{
	return kCDXML_arrow;
}

void CDXArrow::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Arrowhead_Type:
		m_arrowHeadType = (CDXArrowHeadType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Head:
		m_arrowHeadHead = (CDXArrowHeadPosition) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Tail:
		m_arrowHeadTail = (CDXArrowHeadPosition) src_arg.GetUINT(size_arg);
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
	case kCDXProp_Arc_AngularSize:
		m_angularSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_ShaftSpacing:
		m_shaftSpacing = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_EquilibriumRatio:
		m_equilibriumRatio = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_SourceID:
		m_sourceID = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_TargetID:
		m_targetID = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_NoGo:
		m_nogo = (CDXNoGoType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_Dipole:
		m_dipole = true;
		break;
	case kCDXProp_3DHead:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_3dHead.x = src_arg.GetINT32();
		m_3dHead.y = src_arg.GetINT32();
		m_3dHead.z = src_arg.GetINT32();
		m_has3dHead = true;
		break;
		}
	case kCDXProp_3DTail:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_3dTail.x = src_arg.GetINT32();
		m_3dTail.y = src_arg.GetINT32();
		m_3dTail.z = src_arg.GetINT32();
		m_has3dTail = true;
		break;
		}
	case kCDXProp_3DCenter:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_center.x = src_arg.GetINT32();
		m_center.y = src_arg.GetINT32();
		m_center.z = src_arg.GetINT32();
		m_hasCenter = true;
		break;
		}
	case kCDXProp_3DMajorAxisEnd:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_majorAxisEnd.x = src_arg.GetINT32();
		m_majorAxisEnd.y = src_arg.GetINT32();
		m_majorAxisEnd.z = src_arg.GetINT32();
		m_hasMajorAxisEnd = true;
		break;
		}
	case kCDXProp_3DMinorAxisEnd:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_minorAxisEnd.x = src_arg.GetINT32();
		m_minorAxisEnd.y = src_arg.GetINT32();
		m_minorAxisEnd.z = src_arg.GetINT32();
		m_hasMinorAxisEnd = true;
		break;
		}
	case kCDXProp_FadePercent:
		{
		m_fadePercent = src_arg.GetUINT(size_arg);
		break;
		}
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXArrow::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Arrowhead_Type:
		m_arrowHeadType = sXMLArrowHeadType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Head:
		m_arrowHeadHead = sXMLArrowHeadPositionType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Tail:
		m_arrowHeadTail = sXMLArrowHeadPositionType.lookup(attribValue_arg);
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
	case kCDXProp_Arc_AngularSize:
		m_angularSize = (int)(atof(attribValue_arg.c_str()) * 10.0);
		break;
	case kCDXProp_Arrow_ShaftSpacing:
		m_shaftSpacing = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrow_EquilibriumRatio:
		m_equilibriumRatio = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrow_SourceID:
		m_sourceID = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrow_TargetID:
		m_targetID = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arrow_NoGo:
		m_nogo = sXMLNoGoType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrow_Dipole:
		m_dipole = ((attribValue_arg == "yes") || (attribValue_arg == "true"));
		break;
	case kCDXProp_3DHead:
		m_3dHead = StringToCDXPoint3D(attribValue_arg);
		m_has3dHead = true;
		break;
	case kCDXProp_3DTail:
		m_3dTail = StringToCDXPoint3D(attribValue_arg);
		m_has3dTail = true;
		break;
	case kCDXProp_3DCenter:
		m_center = StringToCDXPoint3D(attribValue_arg);
		m_hasCenter = true;
		break;
	case kCDXProp_3DMajorAxisEnd:
		m_majorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		m_hasMajorAxisEnd = true;
		break;
	case kCDXProp_3DMinorAxisEnd:
		m_minorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		m_hasMinorAxisEnd = true;
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = atoi(attribValue_arg.c_str());
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXArrow::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_arrowHeadType != kCDXArrowHeadType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Type, (INT16) m_arrowHeadType );

	if (m_arrowHeadHead != kCDXArrowHeadPosition_Unspecified && m_arrowHeadHead != kCDXArrowHeadPosition_None)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Head, (INT16) m_arrowHeadHead );

	if (m_arrowHeadTail != kCDXArrowHeadPosition_Unspecified && m_arrowHeadTail != kCDXArrowHeadPosition_None)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Tail, (INT16) m_arrowHeadTail );

	if (m_headSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Size, m_headSize );

	if (m_headCenterSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_CenterSize, m_headCenterSize );

	if (m_headWidth != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Width, m_headWidth );

	if (m_angularSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arc_AngularSize, m_angularSize );

	if (m_shaftSpacing != 0)
		sink_arg.PutAttribute( kCDXProp_Arrow_ShaftSpacing, m_shaftSpacing );

	if (m_equilibriumRatio != 100)
		sink_arg.PutAttribute( kCDXProp_Arrow_EquilibriumRatio, m_equilibriumRatio );

	if (m_sourceID != -1)
		sink_arg.PutAttribute( kCDXProp_Arrow_SourceID, m_sourceID );

	if (m_targetID != -1)
		sink_arg.PutAttribute( kCDXProp_Arrow_TargetID, m_targetID );

	if (m_nogo != kCDXNoGoType_Unspecified && m_nogo != kCDXNoGoType_None)
		sink_arg.PutAttribute(kCDXProp_Arrow_NoGo, (INT8) m_nogo);

	if (m_dipole)
		sink_arg.PutAttribute(kCDXProp_Arrow_Dipole);

	if (m_has3dHead)
		sink_arg.PutAttribute( kCDXProp_3DHead, m_3dHead );

	if (m_has3dTail)
		sink_arg.PutAttribute( kCDXProp_3DTail, m_3dTail );

	if (m_hasCenter)
		sink_arg.PutAttribute( kCDXProp_3DCenter, m_center );

	if (m_hasMajorAxisEnd)
		sink_arg.PutAttribute( kCDXProp_3DMajorAxisEnd, m_majorAxisEnd );

	if (m_hasMinorAxisEnd)
		sink_arg.PutAttribute( kCDXProp_3DMinorAxisEnd, m_minorAxisEnd );

	if (m_fadePercent != 1000)
		sink_arg.PutAttribute( kCDXProp_FadePercent, m_fadePercent );
}

void CDXArrow::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_arrowHeadHead != kCDXArrowHeadPosition_Unspecified && m_arrowHeadHead != kCDXArrowHeadPosition_None)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Head, m_arrowHeadHead );

	if (m_arrowHeadTail != kCDXArrowHeadPosition_Unspecified && m_arrowHeadTail != kCDXArrowHeadPosition_None)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Tail, m_arrowHeadTail );

	if (m_arrowHeadType != kCDXArrowHeadType_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Type, m_arrowHeadType );

	if (m_headSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Size, m_headSize );

	if (m_headCenterSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_CenterSize, m_headCenterSize );

	if (m_headWidth != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Width, m_headWidth );

	if (m_angularSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arc_AngularSize, m_angularSize / 10.0 );

	if (m_shaftSpacing != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_ShaftSpacing, m_shaftSpacing );

	if (m_equilibriumRatio != 100)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_EquilibriumRatio, m_equilibriumRatio );

	if (m_sourceID != -1)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_SourceID, m_sourceID );

	if (m_targetID != -1)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_TargetID, m_targetID );

	if (m_nogo != kCDXNoGoType_Unspecified && m_nogo != kCDXNoGoType_None)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_NoGo, m_nogo);

	if (m_dipole)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_Dipole);

	if (m_has3dHead)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DHead, m_3dHead );

	if (m_has3dTail)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DTail, m_3dTail );

	if (m_hasCenter)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DCenter, m_center );

	if (m_hasMajorAxisEnd)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMajorAxisEnd, m_majorAxisEnd );

	if (m_hasMinorAxisEnd)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMinorAxisEnd, m_minorAxisEnd );

	if (m_fadePercent != 1000)
		CDXMLPutAttribute(sink_arg, kCDXProp_FadePercent, m_fadePercent );
}
