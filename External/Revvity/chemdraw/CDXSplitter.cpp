// CommonCS/LibCommon/Src/CDXSplitter.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 2001-2004, CambridgeSoft Corp., All Rights Reserved

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



