// CommonCS/LibCommon/Src/CDXGEPBand.cpp
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
#include "CDXMLNames.h"

CDXGEPBand::CDXGEPBand(CDXObjectID id)
	: CDXPlateItemBase(id)
	, m_mass(0)
{
};

CDXGEPBand::CDXGEPBand(const CDXGEPBand &a)
	: CDXPlateItemBase(a)
	, m_mass(a.m_mass)
{
};

void CDXGEPBand::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	if (attribTag_arg == kCDXProp_GEP_Value)
	{
		m_mass = src_arg.GetINT32();
	}
	else if (attribTag_arg == kCDXProp_GEP_ShowValue)
	{
		m_showValue = true;
	}
	else
	{
		CDXPlateItemBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
	}
}

std::string CDXGEPBand::XMLObjectName() const
{
	return kCDXML_gepband;
}

CDXObject*	CDXGEPBand::Clone() const
{
	return new CDXGEPBand(*this);
}

void CDXGEPBand::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	if (attribTag_arg == kCDXProp_GEP_Value)
	{
		m_mass = atoi(attribValue_arg.c_str());
	}
	else if (attribTag_arg == kCDXProp_GEP_ShowValue)
	{
		m_showValue = attribValue_arg == "yes" || attribValue_arg == "true";
	}
	else
	{
		CDXPlateItemBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
	}
}

void CDXGEPBand::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute(kCDXProp_GEP_Value, m_mass);
	
	if (m_showValue)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_ShowValue);
	}

	CDXPlateItemBase::WriteAttributesTo(sink_arg);
}

void CDXGEPBand::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute(sink_arg, kCDXProp_GEP_Value, m_mass);
	
	if (m_showValue)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_ShowValue);
	}

	CDXPlateItemBase::XMLWriteAttributes(sink_arg);
}
