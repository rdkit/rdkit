// CommonCS/LibCommon/Src/CDXGroupOrFragment.cpp
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
#include <sstream>

XMLPutEnum<CDXSeqType>::Value sXMLSequenceTypeValues[] = {
	{kCDXSeqType_Unknown,	"Unknown"},
    {kCDXSeqType_Peptide,	"Peptide" },
	{kCDXSeqType_Peptide1,	"Peptide1"},
	{kCDXSeqType_Peptide3,	"Peptide3"},
	{kCDXSeqType_DNA,		"DNA"},
	{kCDXSeqType_RNA,		"RNA"},
	{kCDXSeqType_Biopolymer,"Biopolymer"}
};

XMLPutEnum<CDXSeqType> sXMLSequenceType(sXMLSequenceTypeValues, sizeof sXMLSequenceTypeValues, kCDXSeqType_Unknown);
void XMLPut(XMLDataSink &sink, CDXSeqType  v) {	sink.os << sXMLSequenceType.lookup(v); }

// ******************************
// ** class CDXGroupOrFragment **
// ******************************
//
// There are some operations for which Fragments or Groups of Fragments should be permitted,
// so this superclass is used to encompass them both.

CDXGroupOrFragment::CDXGroupOrFragment(CDXTag tag, CDXObjectID id_arg)
 : CDXGraphicObject(tag, id_arg), m_integral(false)
 , m_sequenceType(kCDXSeqType_Unknown)
{
}

CDXGroupOrFragment::~CDXGroupOrFragment()
{
}

CDXGroupOrFragment::CDXGroupOrFragment (const CDXGroupOrFragment &src)
	:	CDXGraphicObject(src)
	,	CDX_INHERIT_SRC (m_connectionOrdering)
	,	CDX_INHERIT_SRC (m_sequenceType)
	,	CDX_INHERIT_SRC (m_integral)
{
}

void CDXGroupOrFragment::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Frag_ConnectionOrder:
		CDXReadItems(m_connectionOrdering, size_arg/sizeof(UINT32), src_arg);
		break;

	case kCDXProp_Frag_SequenceType:
		m_sequenceType = (CDXSeqType) src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_Group_Integral:
		m_integral = (src_arg.GetUINT8() != 0);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXGroupOrFragment::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Frag_ConnectionOrder:
	{
		std::istringstream is(attribValue_arg);
		m_connectionOrdering.clear();
		std::copy(std::istream_iterator<CDXObjectID>(is),
				  std::istream_iterator<CDXObjectID>(),
				  std::back_inserter(m_connectionOrdering));
		break;
	}

	case kCDXProp_Frag_SequenceType:
		m_sequenceType = sXMLSequenceType.lookup(attribValue_arg);
		m_integral = attribValue_arg == "yes";
		break;

	case kCDXProp_Group_Integral:
		m_integral = attribValue_arg == "yes";
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXGroupOrFragment::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (!m_connectionOrdering.empty())
	{
		sink_arg.PutTag( kCDXProp_Frag_ConnectionOrder );
		sink_arg.Put( UINT16( 4 * m_connectionOrdering.size() ) );
		Put( sink_arg, m_connectionOrdering );
	}
    
    if (m_sequenceType != kCDXSeqType_Unknown)
    {
        sink_arg.PutAttribute(kCDXProp_Frag_SequenceType, INT8(m_sequenceType));
    }

    if (m_integral)
    {
        sink_arg.PutAttribute(kCDXProp_Group_Integral, UINT8(m_integral));
    }
}

void CDXGroupOrFragment::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (!m_connectionOrdering.empty())
    {
		CDXMLPutAttribute(sink_arg, kCDXProp_Frag_ConnectionOrder, m_connectionOrdering);
    }

    if (m_sequenceType != kCDXSeqType_Unknown)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Frag_SequenceType, m_sequenceType);
    }

    if (m_integral)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Group_Integral, (bool)m_integral);
    }
}

void CDXGroupOrFragment::AddAttachmentNode (const CDXNode *pNode)
{
	m_connectionOrdering.push_back (pNode->GetObjectID());
}

CDXNode *CDXGroupOrFragment::GetAttachmentNode(size_t attachmentIndex)
{
	// Return the attachmentIndex'th external connection point node
	// Might return NULL in a variety of error cases
	
	if (m_connectionOrdering.size() == 0)
	{
		// We don't have an ordering, so go searching
		CDXObjectsRange nodes = ContainedObjects(kCDXObj_Node);
		for (CDXObjectsByTag::const_iterator i = nodes.begin();  i != nodes.end();  ++i)
		{
			CDXNode *n = dynamic_cast<CDXNode *>(GetObject(i));
			if (n != NULL && n->m_nodeType == kCDXNodeType_ExternalConnectionPoint)
			{
				if (attachmentIndex == 0)
					return n;
				--attachmentIndex;
			}
		}
		return 0;
	}
	else if (attachmentIndex >= m_connectionOrdering.size())
		return 0;
	else
		return dynamic_cast<CDXNode *>(FindByID(m_connectionOrdering[attachmentIndex]));
}
