// CommonCS/LibCommon/Src/CDXTable.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 2001-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <sstream>

// ***********************
// **  class CDXTable   **
// ***********************
//
// Specialization of CDXObject for CDXTable objects

CDXTable::CDXTable(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_Table, id_arg)
{
}

CDXTable::CDXTable (const CDXTable &src)
	:	CDXGraphicObject(src)
{
}

CDXObject*	CDXTable::Clone() const
{
	return new CDXTable(*this);
}

std::string CDXTable::XMLObjectName() const
{
	return kCDXML_table;
}

void CDXTable::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
//	switch (attribTag_arg)
	{
//	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
//		break;
	}
}

void CDXTable::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
//	switch (attribTag_arg)
	{
//	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
//		break;
	}
}

void CDXTable::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

}

void CDXTable::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}



