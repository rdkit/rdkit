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
#include "cs_stringUtils.h"

using namespace cs;

class CDXRLogic;
class CDXRLogicItem;

// ********************
// ** class CDXRLogic **
// ********************
//
// Specialization of CDXObject for CDXRLogic objects

CDXRLogic::CDXRLogic(CDXObjectID id_arg)
 : CDXText(id_arg)
{
}

CDXRLogic::CDXRLogic (const CDXRLogic &src)
	:	CDXText		(src)
{
}


CDXObject*	CDXRLogic::Clone() const
{
	return new CDXRLogic (*this);
}

std::string CDXRLogic::XMLObjectName() const
{
	return kCDXML_rlogic;
}

void CDXRLogic::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	CDXText::StoreAttribute(src_arg, attribTag_arg, size_arg);
}

void CDXRLogic::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	CDXText::XMLStoreAttribute(attribTag_arg, attribValue_arg);
}

void CDXRLogic::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXText::WriteAttributesTo(sink_arg);

}

void CDXRLogic::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXText::XMLWriteAttributes(sink_arg);
}


// ********************
// ** class CDXRLogicItem **
// ********************
//
// Specialization of CDXObject for CDXRLogicItem objects

CDXRLogicItem::CDXRLogicItem(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_RLogic, id_arg)
 , m_group("")
 , m_occurrence("")
 , m_itgroup("")
 , m_restH(false)
{
}

CDXRLogicItem::CDXRLogicItem (const CDXRLogicItem &src)
	:	CDXGraphicObject		(src)
	, m_group(src.m_group)
	, m_occurrence(src.m_occurrence)
	, m_itgroup(src.m_itgroup)
	, m_restH(src.m_restH)
{
}

CDXObject*	CDXRLogicItem::Clone() const
{
	return new CDXRLogicItem (*this);
}



std::string CDXRLogicItem::XMLObjectName() const
{
	return kCDXML_rlogicitem;
}

void CDXRLogicItem::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_RLogic_Group:
		m_group = src_arg.GetString(size_arg);
		break;
	case kCDXProp_RLogic_Occurrence:
		m_occurrence = src_arg.GetString(size_arg);
		break;
	case kCDXProp_RLogic_IfThenGroup:
		m_itgroup = src_arg.GetString(size_arg);
		break;
	case kCDXProp_RLogic_RestH:
		m_restH = true;
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXRLogicItem::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_RLogic_Group:
		m_group=attribValue_arg;
	break;
	case kCDXProp_RLogic_Occurrence:
		m_occurrence=attribValue_arg;
	break;
	case kCDXProp_RLogic_IfThenGroup:
		m_itgroup=attribValue_arg;
	break;
	case kCDXProp_RLogic_RestH:
		m_restH = true;
	break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXRLogicItem::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttributeString( kCDXProp_RLogic_Group, m_group );
	sink_arg.PutAttributeString( kCDXProp_RLogic_Occurrence, m_occurrence );
	sink_arg.PutAttributeString( kCDXProp_RLogic_IfThenGroup, m_itgroup );
	if (m_restH)
		sink_arg.PutAttribute( kCDXProp_RLogic_RestH );
	CDXGraphicObject::WriteAttributesTo(sink_arg);

}

void CDXRLogicItem::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute( sink_arg, kCDXProp_RLogic_Group, m_group );
	CDXMLPutAttribute( sink_arg, kCDXProp_RLogic_Occurrence, m_occurrence );
	CDXMLPutAttribute( sink_arg, kCDXProp_RLogic_IfThenGroup, m_itgroup );
	if (m_restH)
		CDXMLPutAttribute( sink_arg, kCDXProp_RLogic_RestH );
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

}

