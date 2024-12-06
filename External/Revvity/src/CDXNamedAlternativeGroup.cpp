// CommonCS/LibCommon/Src/CDXNamedAlternativeGroup.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <ostream>

// ************************************
// ** class CDXNamedAlternativeGroup **
// ************************************
//
// Specialization of CDXObject for CDXNamedAlternativeGroup objects

const UINT16 CDXNamedAlternativeGroup::has_textFrame  = 0x0001;
const UINT16 CDXNamedAlternativeGroup::has_groupFrame = 0x0002;

CDXNamedAlternativeGroup::CDXNamedAlternativeGroup(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_NamedAlternativeGroup, id_arg)
 , m_altGroupFlags(0)
 , m_valence(0)
{
}

CDXNamedAlternativeGroup::CDXNamedAlternativeGroup (const CDXNamedAlternativeGroup &src)
	:	CDXGraphicObject(src)
	,	m_textFrame(src.m_textFrame)
	,	m_groupFrame(src.m_groupFrame)
	,	m_altGroupFlags(src.m_altGroupFlags)
	,	m_valence(src.m_valence)
{
}

CDXNamedAlternativeGroup::~CDXNamedAlternativeGroup()
{
}

CDXObject*	CDXNamedAlternativeGroup::Clone() const
{
	return new CDXNamedAlternativeGroup (*this);
}

std::string CDXNamedAlternativeGroup::XMLObjectName() const
{
	return kCDXML_altgroup;
}

void CDXNamedAlternativeGroup::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_NamedAlternativeGroup_TextFrame:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetUINT32();
		CDXCoordinate left   = src_arg.GetUINT32();
		CDXCoordinate bottom = src_arg.GetUINT32();
		CDXCoordinate right  = src_arg.GetUINT32();
		TextFrame(CDXRectangle(top, left, bottom, right));
		break;
		}
		break;
	case kCDXProp_NamedAlternativeGroup_GroupFrame:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetUINT32();
		CDXCoordinate left   = src_arg.GetUINT32();
		CDXCoordinate bottom = src_arg.GetUINT32();
		CDXCoordinate right  = src_arg.GetUINT32();
		GroupFrame(CDXRectangle(top, left, bottom, right));
		break;
		}
		break;
	case kCDXProp_NamedAlternativeGroup_Valence:
		m_valence = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXNamedAlternativeGroup::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_NamedAlternativeGroup_TextFrame:
		TextFrame(StringToCDXRectangle(attribValue_arg));
		break;
	case kCDXProp_NamedAlternativeGroup_GroupFrame:
		GroupFrame(StringToCDXRectangle(attribValue_arg));
		break;
	case kCDXProp_NamedAlternativeGroup_Valence:
		m_valence = atoi(attribValue_arg.c_str());
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXNamedAlternativeGroup::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (KnownTextFrame())
		sink_arg.PutAttribute( kCDXProp_NamedAlternativeGroup_TextFrame, m_textFrame );

	if (KnownGroupFrame())
		sink_arg.PutAttribute( kCDXProp_NamedAlternativeGroup_GroupFrame, m_groupFrame );

	if (m_valence != 0)
		sink_arg.PutAttribute( kCDXProp_NamedAlternativeGroup_Valence, m_valence );
}

void CDXNamedAlternativeGroup::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (KnownTextFrame())
		CDXMLPutAttribute(sink_arg, kCDXProp_NamedAlternativeGroup_TextFrame, m_textFrame );

	if (KnownGroupFrame())
		CDXMLPutAttribute(sink_arg, kCDXProp_NamedAlternativeGroup_GroupFrame, m_groupFrame );

	if (m_valence != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_NamedAlternativeGroup_Valence, m_valence );

}

