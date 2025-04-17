// CommonCS/LibCommon/Src/CDXTLCSpot.cpp
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
#include "CDXMLNames.h"


void CDXTLCSpot::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch(attribTag_arg)
	{
	case kCDXProp_TLC_Rf:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_value = SwapBytes(f64);
		break;

	case kCDXProp_TLC_ShowRf:
		m_showValue = true;
		break;

	case kCDXProp_TLC_Tail:
		m_tail = src_arg.GetUINT(size_arg);
		break;

	default:
		CDXPlateItemBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
	}
}

std::string CDXTLCSpot::XMLObjectName() const
{
	return kCDXML_tlcspot;
}

CDXObject*	CDXTLCSpot::Clone() const
{
	return new CDXTLCSpot(*this);
}

CDXTLCSpot::CDXTLCSpot(CDXObjectID id)
	:	CDXPlateItemBase(id)
	,	m_tail(0) 
{
}

CDXTLCSpot::CDXTLCSpot(const CDXTLCSpot &a)
	:	CDXPlateItemBase(a)
	,	m_tail(a.m_tail) 
{
}

void CDXTLCSpot::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_TLC_Rf:
		m_value = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_TLC_ShowRf:
		m_showValue = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_TLC_Tail:
		m_tail = atoi(attribValue_arg.c_str());
		break;

	default:
		CDXPlateItemBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
	}
}

void CDXTLCSpot::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	if (m_value >= 0)
	{
		sink_arg.PutAttribute(kCDXProp_TLC_Rf, FLOAT64(m_value));
	}
	
	if (m_showValue)
	{
		sink_arg.PutAttribute(kCDXProp_TLC_ShowRf);
	}

	sink_arg.PutAttribute(kCDXProp_TLC_Tail, m_tail);
	CDXPlateItemBase::WriteAttributesTo(sink_arg);
}

void CDXTLCSpot::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	if (m_value >= 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_Rf, m_value);
	}
	
	if (m_showValue)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_ShowRf);
	}

	CDXMLPutAttribute(sink_arg, kCDXProp_TLC_Tail, m_tail);
	CDXPlateItemBase::XMLWriteAttributes(sink_arg);
}