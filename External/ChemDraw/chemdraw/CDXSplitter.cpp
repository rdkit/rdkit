// CommonCS/LibCommon/Src/CDXSplitter.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 2001-2004, CambridgeSoft Corp., All Rights Reserved

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
#include <sstream>

extern XMLPutEnum<CDXPageDefinition> sXMLPageDefinitionPage;

// ***********************
// ** class CDXSplitter **
// ***********************
//
// Specialization of CDXObject for CDXSplitter objects

CDXSplitter::CDXSplitter(CDXObjectID id_arg)
	:	CDXObject(kCDXObj_Splitter, id_arg)
	,	m_def(kCDXPageDefinition_Undefined)
{
}

CDXSplitter::CDXSplitter (const CDXSplitter &src)
	:	CDXObject(src)
	,	CDX_INHERIT_SRC (m_pos)
	,	CDX_INHERIT_SRC (m_def)
{
}

CDXObject*	CDXSplitter::Clone() const
{
	return new CDXSplitter(*this);
}

std::string CDXSplitter::XMLObjectName() const
{
	return kCDXML_splitter;
}

void CDXSplitter::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_2DPosition:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_pos = CDXPoint2D(x, y);
		break;
		}
	case kCDXProp_PageDefinition:
		m_def=(CDXPageDefinition)src_arg.GetUINT8();
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXSplitter::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_2DPosition:
		m_pos = StringToCDXPoint2D(attribValue_arg);
		break;
	case kCDXProp_PageDefinition:
		m_def = sXMLPageDefinitionPage.lookup(attribValue_arg);
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXSplitter::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_2DPosition, m_pos );
	if (m_def != kCDXPageDefinition_Undefined)
		sink_arg.PutAttribute( kCDXProp_PageDefinition, UINT8(m_def));

}

void CDXSplitter::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	CDXMLPutAttribute(sink_arg, kCDXProp_2DPosition, m_pos );
	if (m_def != kCDXPageDefinition_Undefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_PageDefinition, m_def);
}



