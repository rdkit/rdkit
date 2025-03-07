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

class CDXPlasmidMarker;
class CDXPlasmidRegion;
class CDXPlasmidMap;

// ********************
// ** class CDXPlasmidMarker **
// ********************
//
// Specialization of CDXObject for CDXPlasmidMarker objects

CDXPlasmidMarker::CDXPlasmidMarker(CDXObjectID id_arg)
 : CDXObjectTag(id_arg)
 , m_markerOffset(0)
 , m_markerAngle(0)
{
}

CDXPlasmidMarker::CDXPlasmidMarker (const CDXPlasmidMarker &src)
	:	CDXObjectTag		(src)
	, m_markerOffset(src.m_markerOffset)
	, m_markerAngle(src.m_markerAngle)
{
}


CDXObject*	CDXPlasmidMarker::Clone() const
{
	return new CDXPlasmidMarker (*this);
}

std::string CDXPlasmidMarker::XMLObjectName() const
{
	return kCDXML_plasmidmarker;
}

void CDXPlasmidMarker::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_MarkerOffset:
		m_markerOffset = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_PlasmidMap_MarkerAngle:
		m_markerAngle = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXObjectTag::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlasmidMarker::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_MarkerOffset:
		m_markerOffset = StrToNum(attribValue_arg);
		break;

	case kCDXProp_PlasmidMap_MarkerAngle:
		m_markerAngle = StrToNum(attribValue_arg);
		break;

	default:
		CDXObjectTag::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlasmidMarker::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_MarkerOffset, m_markerOffset );
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_MarkerAngle, m_markerAngle );
	CDXObjectTag::WriteAttributesTo(sink_arg);

}

void CDXPlasmidMarker::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_MarkerOffset, m_markerOffset );
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_MarkerAngle, m_markerAngle );
	CDXObjectTag::XMLWriteAttributes(sink_arg);

}


// ********************
// ** class CDXPlasmidRegion **
// ********************
//
// Specialization of CDXObject for CDXPlasmidRegion objects

CDXPlasmidRegion::CDXPlasmidRegion(CDXObjectID id_arg)
 : CDXArrow(id_arg)
 , m_regionStart(0)
 , m_regionEnd(0)
 , m_regionOffset(0)
{
}

CDXPlasmidRegion::CDXPlasmidRegion (const CDXPlasmidRegion &src)
	:	CDXArrow		(src)
	, m_regionStart(src.m_regionStart)
	, m_regionEnd(src.m_regionEnd)
	, m_regionOffset(src.m_regionOffset)
{
}


CDXObject*	CDXPlasmidRegion::Clone() const
{
	return new CDXPlasmidRegion (*this);
}

std::string CDXPlasmidRegion::XMLObjectName() const
{
	return kCDXML_plasmidregion;
}

void CDXPlasmidRegion::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_RegionStart:
		m_regionStart = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_PlasmidMap_RegionEnd:
		m_regionEnd = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_PlasmidMap_RegionOffset:
		m_regionOffset = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXArrow::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlasmidRegion::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_RegionStart:
		m_regionStart = StrToNum(attribValue_arg);
		break;
	case kCDXProp_PlasmidMap_RegionEnd:
		m_regionEnd = StrToNum(attribValue_arg);
		break;
	case kCDXProp_PlasmidMap_RegionOffset:
		m_regionOffset = StrToNum(attribValue_arg);
		break;

	default:
		CDXArrow::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlasmidRegion::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_RegionStart, m_regionStart );
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_RegionEnd, m_regionEnd );
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_RegionOffset, m_regionOffset );
	CDXArrow::WriteAttributesTo(sink_arg);

}

void CDXPlasmidRegion::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_RegionStart, m_regionStart );
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_RegionEnd, m_regionEnd );
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_RegionOffset, m_regionOffset );
	CDXArrow::XMLWriteAttributes(sink_arg);

}



// ********************
// ** class CDXPlasmidMap **
// ********************
//
// Specialization of CDXObject for CDXPlasmidMap objects

CDXPlasmidMap::CDXPlasmidMap(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_PlasmidMap, id_arg)
 , m_numberBasePairs(0)
 , m_ringRadius(0)
{
}

CDXPlasmidMap::CDXPlasmidMap (const CDXPlasmidMap &src)
	:	CDXGraphicObject		(src)
	, m_numberBasePairs(src.m_numberBasePairs)
	, m_ringRadius(src.m_ringRadius)
{
}


CDXObject*	CDXPlasmidMap::Clone() const
{
	return new CDXPlasmidMap (*this);
}

std::string CDXPlasmidMap::XMLObjectName() const
{
	return kCDXML_plasmidmap;
}

void CDXPlasmidMap::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_NumberBasePairs:
		m_numberBasePairs = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_PlasmidMap_RingRadius:
		m_ringRadius = src_arg.GetUINT(size_arg);
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlasmidMap::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PlasmidMap_NumberBasePairs:
		m_numberBasePairs = StrToNum(attribValue_arg);
		break;
	case kCDXProp_PlasmidMap_RingRadius:
		m_ringRadius = StrToNum(attribValue_arg);
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlasmidMap::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_NumberBasePairs, m_numberBasePairs );
	sink_arg.PutAttribute( kCDXProp_PlasmidMap_RingRadius, m_ringRadius );
	CDXGraphicObject::WriteAttributesTo(sink_arg);

}

void CDXPlasmidMap::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_NumberBasePairs, m_numberBasePairs );
	CDXMLPutAttribute( sink_arg, kCDXProp_PlasmidMap_RingRadius, m_ringRadius );
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

}
