// CommonCS/LibCommon/Src/CDXPlateBase.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

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


// ***********************
// ** class CDXPlateBase **
// ***********************
//
// Specialization of CDXObject for CDXPlateBase objects

CDXPlateBase::CDXPlateBase(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_TLCPlate, id_arg)
	,	m_showBorders(false)
	,	m_transparent(false)
	,	m_hasTransparent(false)
{
}

CDXPlateBase::CDXPlateBase (const CDXPlateBase &src)
	:	CDXGraphicObject(src)
	,	m_topleft(src.m_topleft)
	,	m_topright(src.m_topright)
	,	m_bottomright(src.m_bottomright)
	,	m_bottomleft(src.m_bottomleft)
	,	m_showBorders(src.m_showBorders)
	,	m_transparent(src.m_transparent)
	,	m_hasTransparent(src.m_hasTransparent)
{
}

void CDXPlateBase::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_TopLeft:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_topleft.y = src_arg.GetINT32();
		m_topleft.x = src_arg.GetINT32();
		break;
		}
	case kCDXProp_TopRight:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_topright.y = src_arg.GetINT32();
		m_topright.x = src_arg.GetINT32();
		break;
		}
	case kCDXProp_BottomRight:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_bottomright.y = src_arg.GetINT32();
		m_bottomright.x = src_arg.GetINT32();
		break;
		}
	case kCDXProp_BottomLeft:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_bottomleft.y = src_arg.GetINT32();
		m_bottomleft.x = src_arg.GetINT32();
		break;
		}

	case kCDXProp_ShowBorders:
		m_showBorders = true;
		break;

	case kCDXProp_Transparent:
		m_transparent = true;
		m_hasTransparent = true;
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlateBase::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_TopLeft:
		m_topleft = StringToCDXPoint2D(attribValue_arg);
		break;

	case kCDXProp_TopRight:
		m_topright = StringToCDXPoint2D(attribValue_arg);
		break;

	case kCDXProp_BottomRight:
		m_bottomright = StringToCDXPoint2D(attribValue_arg);
		break;

	case kCDXProp_BottomLeft:
		m_bottomleft = StringToCDXPoint2D(attribValue_arg);
		break;

	case kCDXProp_ShowBorders:
		m_showBorders = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_Transparent:
		m_transparent = attribValue_arg == "yes" || attribValue_arg == "true";
		m_hasTransparent = true;
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlateBase::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_TopLeft, m_topleft );
	sink_arg.PutAttribute( kCDXProp_TopRight, m_topright );
	sink_arg.PutAttribute( kCDXProp_BottomRight, m_bottomright );
	sink_arg.PutAttribute( kCDXProp_BottomLeft, m_bottomleft );

	if (m_showBorders)
		sink_arg.PutAttribute(kCDXProp_ShowBorders);

	if (m_hasTransparent && m_transparent)
		sink_arg.PutAttribute(kCDXProp_Transparent);

	CDXGraphicObject::WriteAttributesTo(sink_arg);
}

void CDXPlateBase::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute(sink_arg, kCDXProp_TopLeft, m_topleft );
	CDXMLPutAttribute(sink_arg, kCDXProp_TopRight, m_topright );
	CDXMLPutAttribute(sink_arg, kCDXProp_BottomRight, m_bottomright );
	CDXMLPutAttribute(sink_arg, kCDXProp_BottomLeft, m_bottomleft );

	if (m_showBorders)
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowBorders);

	if (m_hasTransparent && m_transparent)
		CDXMLPutAttribute(sink_arg, kCDXProp_Transparent);
	
	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}