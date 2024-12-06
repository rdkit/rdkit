// CommonCS/LibCommon/Src/CDXGroup.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"

// ********************
// ** class CDXGroup **
// ********************
//
// Specialization of CDXObject for CDXGroup objects

CDXGroup::CDXGroup(CDXObjectID id_arg)
 : CDXGroupOrFragment(kCDXObj_Group, id_arg)
{
}

CDXGroup::~CDXGroup()
{
}

CDXObject*	CDXGroup::Clone() const
{
	return new CDXGroup (*this);
}

CDXGroup::CDXGroup(const CDXGroup &src)
	:	CDXGroupOrFragment (src)
{
}

CDXGroup::CDXGroup(const CDXFragment &src)
	:	CDXGroupOrFragment (src)
{
}

std::string CDXGroup::XMLObjectName() const
{
	return kCDXML_group;
}
