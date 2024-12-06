// CommonCS/LibCommon/Src/CDXTemplateGrid.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <ostream>

// ***************************
// ** class CDXTemplateGrid **
// ***************************
//
// Specialization of CDXObject for CDXTemplateGrid objects

CDXTemplateGrid::CDXTemplateGrid(CDXObjectID id_arg)
	: CDXObject(kCDXObj_TemplateGrid, id_arg)
	, m_paneHeight(0)
	, m_numRows(0)
	, m_numColumns(0)
	, m_extent(0, 0)
{
}

CDXTemplateGrid::CDXTemplateGrid (const CDXTemplateGrid &src)
	: CDXObject		(src)
	, m_paneHeight	(src.m_paneHeight)
	, m_numRows		(src.m_numRows)
	, m_numColumns	(src.m_numColumns)
	, m_extent		(src.m_extent)
{
}

CDXTemplateGrid::~CDXTemplateGrid()
{
}

CDXObject*	CDXTemplateGrid::Clone() const
{
	return new CDXTemplateGrid (*this);
}

std::string CDXTemplateGrid::XMLObjectName() const
{
	return kCDXML_templategrid;
}

void CDXTemplateGrid::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Template_PaneHeight:
		m_paneHeight = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Template_NumRows:
		m_numRows = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Template_NumColumns:
		m_numColumns = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_2DExtent:
	{	if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_extent = CDXPoint2D(x, y);
		break;
	}
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXTemplateGrid::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Template_PaneHeight:
		m_paneHeight = atoi(value_arg.c_str());
		break;
	case kCDXProp_Template_NumRows:
		m_numRows = atoi(value_arg.c_str());
		break;
	case kCDXProp_Template_NumColumns:
		m_numColumns = atoi(value_arg.c_str());
		break;
	case kCDXProp_2DExtent:
		m_extent = StringToCDXPoint2D(value_arg);
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, value_arg);
		break;
	}
}

void CDXTemplateGrid::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);
	sink_arg.PutAttribute( kCDXProp_Template_PaneHeight, m_paneHeight );
	sink_arg.PutAttribute( kCDXProp_Template_NumRows, m_numRows );
	sink_arg.PutAttribute( kCDXProp_Template_NumColumns, m_numColumns );
	sink_arg.PutAttribute( kCDXProp_2DExtent, m_extent );
}

void CDXTemplateGrid::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_PaneHeight, m_paneHeight );
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_NumRows, m_numRows );
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_NumColumns, m_numColumns );
	CDXMLPutAttribute( sink_arg, kCDXProp_2DExtent, m_extent );
}
