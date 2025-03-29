// CommonCS/LibCommon/Src/CDXPlateItemBase.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2010, CambridgeSoft Corp., All Rights Reserved

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
// ** class CDXPlateItemBase  **
// ***********************
//
// Specialization of CDXObject for CDXPlateItemBase objects

CDXPlateItemBase::CDXPlateItemBase(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_TLCSpot, id_arg)
	,	m_value(0)
	,	m_width(0)
	,	m_height(0)
	,	m_displayType(0)
	,	m_showValue(false)
{
}

CDXPlateItemBase::CDXPlateItemBase (const CDXPlateItemBase &src)
	:	CDXGraphicObject(src)
	,	m_value(src.m_value)
	,	m_width(src.m_width)
	,	m_height(src.m_height)
	,	m_displayType(src.m_displayType)
	,	m_showValue(src.m_showValue)
{
}

void CDXPlateItemBase::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_Height:
		m_height = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_Curve_Type:
		m_displayType = src_arg.GetUINT(size_arg);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlateItemBase::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_Height:
		m_height = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_Curve_Type:
		m_displayType = atoi(attribValue_arg.c_str());
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlateItemBase::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_Width, m_width );
	sink_arg.PutAttribute( kCDXProp_Height, m_height );

	if (m_displayType != 0)
		sink_arg.PutAttribute( kCDXProp_Curve_Type, m_displayType );

	CDXGraphicObject::WriteAttributesTo(sink_arg);
}

void CDXPlateItemBase::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	streamsize oldPrecision=sink_arg.os.precision(8);
	sink_arg.os.precision(oldPrecision);

	CDXMLPutAttribute( sink_arg, kCDXProp_Width, m_width );
	CDXMLPutAttribute( sink_arg, kCDXProp_Height, m_height );

	if (m_displayType != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Type, m_displayType );

	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}