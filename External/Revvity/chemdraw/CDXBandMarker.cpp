// LibCommon\Src\CDXBandMarker.cpp
// Copyright 2010, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <sstream>
#include "cs_stringUtils.h"

using namespace cs;

CDXBandMarker::CDXBandMarker(CDXObjectID id_arg)
 : CDXObjectTag(id_arg)
{
}

CDXBandMarker::CDXBandMarker (const CDXBandMarker &src)
	:	CDXObjectTag		(src)
{
}


CDXObject*	CDXBandMarker::Clone() const
{
	return new CDXBandMarker (*this);
}

std::string CDXBandMarker::XMLObjectName() const
{
	return kCDXML_marker;
}

void CDXBandMarker::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	CDXObjectTag::StoreAttribute(src_arg, attribTag_arg, size_arg);
}

void CDXBandMarker::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	CDXObjectTag::XMLStoreAttribute(attribTag_arg, attribValue_arg);
}

void CDXBandMarker::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObjectTag::WriteAttributesTo(sink_arg);
}

void CDXBandMarker::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObjectTag::XMLWriteAttributes(sink_arg);
}