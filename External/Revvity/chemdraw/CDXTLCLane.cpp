// CommonCS/LibCommon/Src/CDXTLCLane.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <sstream>

// ***********************
// ** class CDXPlateLane  **
// ***********************
//
// Specialization of CDXObject for CDXPlateLane objects

CDXPlateLane::CDXPlateLane(CDXTag tag, CDXObjectID id)
	:	CDXGraphicObject(tag, id)
{
}

CDXPlateLane::CDXPlateLane (const CDXPlateLane &src)
	:	CDXGraphicObject(src)
{
}

CDXTLCLane::CDXTLCLane(CDXObjectID id)
	: CDXPlateLane(kCDXObj_GEPLane, id)
{
}

CDXTLCLane::CDXTLCLane(const CDXTLCLane &a)
	: CDXPlateLane(a)
{
}

CDXObject*	CDXTLCLane::Clone() const
{
	return new CDXTLCLane(*this);
}

std::string CDXTLCLane::XMLObjectName() const
{
	return kCDXML_tlclane;
}

// CDXGEPLane class

CDXGEPLane::CDXGEPLane(CDXObjectID id)
	: CDXPlateLane(kCDXObj_GEPLane, id)
{
}

CDXGEPLane::CDXGEPLane(const CDXGEPLane &a)
	: CDXPlateLane(a)
{
}

CDXObject*	CDXGEPLane::Clone() const
{
	return new CDXGEPLane(*this);
}

std::string CDXGEPLane::XMLObjectName() const
{
	return kCDXML_geplane;
}

void CDXGEPLane::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	CDXPlateLane::XMLStoreAttribute(attribTag_arg, attribValue_arg);
}

void CDXGEPLane::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXPlateLane::XMLWriteAttributes(sink_arg);
}

void CDXGEPLane::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXPlateLane::WriteAttributesTo(sink_arg);
}

void CDXGEPLane::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	CDXPlateLane::StoreAttribute(src_arg, attribTag_arg, size_arg);
}