// CommonCS/LibCommon/Src/CDXBracketedGroup.cpp
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

extern XMLPutEnum<CDXBracketUsage> sXMLBracketUsage;
extern XMLPutEnum<CDXPolymerRepeatPattern> sXMLPolymerRepeatPattern;
extern XMLPutEnum<CDXPolymerFlipType> sXMLPolymerFlipType;

// *****************************
// ** class CDXBracketedGroup **
// *****************************
//
// Specialization of CDXObject for CDXBracketedGroup objects

CDXBracketedGroup::CDXBracketedGroup(CDXObjectID id_arg)
	: CDXObject	(kCDXObj_BracketedGroup, id_arg)
	, m_usage			(kCDXBracketUsage_Unspecified)
	, m_repeatPattern	(kCDXPolymerRepeatPattern_HeadToTail)
	, m_flipType		(kCDXPolymerFlipType_Unspecified)
	, m_repeatCount		(2)
	, m_componentOrder	(0)
{
}

CDXBracketedGroup::CDXBracketedGroup (const CDXBracketedGroup &src)
	:	CDXObject	(src)
	,	m_bracketedObjects		(src.m_bracketedObjects)
	,	m_usage					(src.m_usage)
	,	m_repeatPattern			(src.m_repeatPattern)
	,	m_flipType				(src.m_flipType)
	,	m_repeatCount			(src.m_repeatCount)
	,	m_componentOrder		(src.m_componentOrder)
	,	m_SRULabel				(src.m_SRULabel)
{
}

CDXBracketedGroup::~CDXBracketedGroup()
{
}

CDXObject*	CDXBracketedGroup::Clone() const
{
	return new CDXBracketedGroup (*this);
}

std::string CDXBracketedGroup::XMLObjectName() const
{
	return kCDXML_bracketedgroup;
}

void CDXBracketedGroup::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_BracketedObjects:
	{
		CDXReadItems(m_bracketedObjects, size_arg/sizeof(UINT32), src_arg);
		break;
	}
	case kCDXProp_Bracket_Usage:
		m_usage = (CDXBracketUsage) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Polymer_RepeatPattern:
		m_repeatPattern = (CDXPolymerRepeatPattern) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Polymer_FlipType:
		m_flipType = (CDXPolymerFlipType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_RepeatCount:
		{
			FLOAT64 f64;
			src_arg.GetBytes((char *)&f64, sizeof(f64));
			m_repeatCount = SwapBytes(f64);
		}
		break;
	case kCDXProp_Bracket_ComponentOrder:
		m_componentOrder = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_SRULabel:
		{
			CDXString cdxs;
			cdxs.Read(src_arg, size_arg);
			m_SRULabel = cdxs.str();
		}
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXBracketedGroup::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
    switch (attribTag_arg)
    {
        case kCDXProp_BracketedObjects:
        {
            std::istringstream is(value_arg);
            m_bracketedObjects.clear();
            std::copy(std::istream_iterator<CDXObjectID>(is),
                std::istream_iterator<CDXObjectID>(),
                std::back_inserter(m_bracketedObjects));
            break;
        }

        case kCDXProp_Bracket_Usage:
            m_usage = sXMLBracketUsage.lookup(value_arg);
            break;

        case kCDXProp_Polymer_RepeatPattern:
            m_repeatPattern = sXMLPolymerRepeatPattern.lookup(value_arg);
            break;

        case kCDXProp_Polymer_FlipType:
            m_flipType = sXMLPolymerFlipType.lookup(value_arg);
            break;

        case kCDXProp_Bracket_RepeatCount:
            m_repeatCount = atof(value_arg.c_str());
            break;

        case kCDXProp_Bracket_ComponentOrder:
            m_componentOrder = atoi(value_arg.c_str());
            break;

        case kCDXProp_Bracket_SRULabel:
            m_SRULabel = UnicodeToString(value_arg);
            break;

        default:
            CDXObject::XMLStoreAttribute(attribTag_arg, value_arg);
            break;
    }
}

void CDXBracketedGroup::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	if (!m_bracketedObjects.empty())
	{
		sink_arg.PutTag( kCDXProp_BracketedObjects );
		sink_arg.Put( UINT16( 4 * m_bracketedObjects.size() ) );
		Put( sink_arg, m_bracketedObjects );
	}

    if (m_usage != kCDXBracketUsage_Unspecified)
    {
        sink_arg.PutAttribute(kCDXProp_Bracket_Usage, INT16(m_usage));
    }

    if (m_repeatPattern != kCDXPolymerRepeatPattern_HeadToTail)
    {
        sink_arg.PutAttribute(kCDXProp_Polymer_RepeatPattern, INT16(m_repeatPattern));
    }

    if (m_flipType != kCDXPolymerFlipType_Unspecified)
    {
        sink_arg.PutAttribute(kCDXProp_Polymer_FlipType, INT16(m_flipType));
    }

    if (m_usage == kCDXBracketUsage_MultipleGroup)
    {
        sink_arg.PutAttribute(kCDXProp_Bracket_RepeatCount, m_repeatCount);
    }

    if (m_componentOrder != 0)
    {
        sink_arg.PutAttribute(kCDXProp_Bracket_ComponentOrder, INT16(m_componentOrder));
    }

    if (!m_SRULabel.empty())
    {
        sink_arg.PutAttributeCDXString(kCDXProp_Bracket_SRULabel, CDXString(m_SRULabel));
    }
}

void CDXBracketedGroup::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

    if (!m_bracketedObjects.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_BracketedObjects, m_bracketedObjects);
    }

    if (m_usage != kCDXBracketUsage_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_Usage, m_usage);
    }

    if (m_repeatPattern != kCDXPolymerRepeatPattern_HeadToTail)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Polymer_RepeatPattern, m_repeatPattern);
    }

    if (m_flipType != kCDXPolymerFlipType_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Polymer_FlipType, m_flipType);
    }

    if (m_usage == kCDXBracketUsage_MultipleGroup)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_RepeatCount, m_repeatCount);
    }

    if (m_componentOrder != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_ComponentOrder, m_componentOrder);
    }

    if (!m_SRULabel.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_SRULabel, StringToUnicode("", m_SRULabel));
    }
}

// ********************************
// ** class CDXBracketAttachment **
// ********************************
//
// Specialization of CDXObject for CDXBracketAttachment objects

CDXBracketAttachment::CDXBracketAttachment(CDXObjectID id_arg)
	: CDXObject	(kCDXObj_BracketAttachment, id_arg)
	, m_graphicID		(0)
{
}

CDXBracketAttachment::CDXBracketAttachment (const CDXBracketAttachment &src)
	:	CDXObject	(src)
	,	m_graphicID		(src.m_graphicID)
{
}

CDXBracketAttachment::~CDXBracketAttachment()
{
}

CDXObject*	CDXBracketAttachment::Clone() const
{
	return new CDXBracketAttachment (*this);
}

std::string CDXBracketAttachment::XMLObjectName() const
{
	return kCDXML_bracketattachment;
}

void CDXBracketAttachment::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bracket_GraphicID:
		m_graphicID = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXBracketAttachment::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bracket_GraphicID:
		m_graphicID = atoi(value_arg.c_str());
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, value_arg);
		break;
	}
}


void CDXBracketAttachment::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_Bracket_GraphicID,			m_graphicID );
}

void CDXBracketAttachment::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_GraphicID,		m_graphicID);
}

// ***************************
// ** class CDXCrossingBond **
// ***************************
//
// Specialization of CDXObject for CDXCrossingBond objects

CDXCrossingBond::CDXCrossingBond(CDXObjectID id_arg)
	: CDXObject	(kCDXObj_CrossingBond, id_arg)
	, m_bondID			(0)
	, m_innerAtomID		(0)
{
}

CDXCrossingBond::CDXCrossingBond (const CDXCrossingBond &src)
	:	CDXObject	(src)
	,	m_bondID		(src.m_bondID)
	,	m_innerAtomID	(src.m_innerAtomID)
{
}

CDXCrossingBond::~CDXCrossingBond()
{
}

CDXObject*	CDXCrossingBond::Clone() const
{
	return new CDXCrossingBond (*this);
}

std::string CDXCrossingBond::XMLObjectName() const
{
	return kCDXML_crossingbond;
}

void CDXCrossingBond::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bracket_BondID:
		m_bondID = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_InnerAtomID:
		m_innerAtomID = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXCrossingBond::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bracket_BondID:
		m_bondID = atoi(value_arg.c_str());
		break;
	case kCDXProp_Bracket_InnerAtomID:
		m_innerAtomID = atoi(value_arg.c_str());
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, value_arg);
		break;
	}
}


void CDXCrossingBond::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_Bracket_BondID,			m_bondID );
	sink_arg.PutAttribute( kCDXProp_Bracket_InnerAtomID,	m_innerAtomID );
}

void CDXCrossingBond::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_BondID,		m_bondID);
	CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_InnerAtomID,	m_innerAtomID);
}

