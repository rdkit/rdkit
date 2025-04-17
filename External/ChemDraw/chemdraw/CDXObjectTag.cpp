// CommonCS/LibCommon/Src/CDXObjectTag.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2007, CambridgeSoft Corp., All Rights Reserved

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
#include "cs_assert.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <stdio.h>	// for sprintf
#include <sstream>

XMLPutEnum<CDXObjectTagType>::Value sXMLObjectTagTypeValues[] = {
	{kCDXObjectTagType_Undefined,		"Unknown"},
	{kCDXObjectTagType_Double,			"Double"},
	{kCDXObjectTagType_Int32,			"Long"},        // Note that "Long" in the file is a signed 32-bit int
	{kCDXObjectTagType_String,			"String"}
};

XMLPutEnum<CDXObjectTagType> sXMLObjectTagType(sXMLObjectTagTypeValues, sizeof sXMLObjectTagTypeValues, kCDXObjectTagType_Undefined);
void XMLPut(XMLDataSink &sink, CDXObjectTagType  v) {	sink.os << sXMLObjectTagType.lookup(v); }


XMLPutEnum<CDXPositioningType>::Value sXMLPositioningTypeValues[] = {
	{kCDXPositioningType_Auto,		"auto"},
	{kCDXPositioningType_Angle,		"angle"},
	{kCDXPositioningType_Offset,	"offset"},
	{kCDXPositioningType_Absolute,	"absolute"}
};

XMLPutEnum<CDXPositioningType> sXMLPositioningType(sXMLPositioningTypeValues, sizeof sXMLPositioningTypeValues, kCDXPositioningType_Auto);
void XMLPut(XMLDataSink &sink, CDXPositioningType  v) {	sink.os << sXMLPositioningType.lookup(v); }


// ************************
// ** class CDXObjectTag **
// ************************
//
// Specialization of CDXObject for CDXObjectTag objects


CDXObject*	CDXObjectTag::Clone() const
{
	return new CDXObjectTag (*this);
}


std::string CDXObjectTag::XMLObjectName() const
{
	return kCDXML_objecttag;
}


CDXObjectTag::CDXObjectTag(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_ObjectTag, id_arg)
	,	m_Type(kCDXObjectTagType_Undefined)
	,	m_Name(kCDXTagType_Unknown)
	,	m_Tracking(false)
	,	m_Persistent(true)
	,	m_positioning(kCDXPositioningType_Auto)
	,	m_positioningAngle(0.0)
	,	m_positioningOffset(CDXPoint2D(0,0))
{
    m_Int32Val = 0;
    m_DoubleVal = 0;
}

CDXObjectTag::CDXObjectTag (const CDXObjectTag &src)
	:	CDXGraphicObject(src)
	,	CDX_INHERIT_SRC (m_Type)
	,	CDX_INHERIT_SRC (m_Name)
	,	CDX_INHERIT_SRC (m_Tracking)
	,	CDX_INHERIT_SRC (m_Persistent)
	,	CDX_INHERIT_SRC (m_positioning)
	,	CDX_INHERIT_SRC (m_positioningAngle)
	,	CDX_INHERIT_SRC (m_positioningOffset)
{
    m_Int32Val = 0;
    m_DoubleVal = 0;
    m_StringVal.clear();

    switch (m_Type)
    {
    case kCDXObjectTagType_Double: m_DoubleVal = src.m_DoubleVal; break;
    case kCDXObjectTagType_Int32:  m_Int32Val = src.m_Int32Val;   break;
    case kCDXObjectTagType_String: m_StringVal = src.m_StringVal; break;
    case kCDXObjectTagType_Undefined:                             break;
    }
}

#if TARGET_OS_WIN32
#pragma warning( push )
#pragma warning( disable : 4800 )	// conversion-to-bool inefficiency warning
#endif

// Take a CDX input stream and create a CDX structure
void CDXObjectTag::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg) {

	case kCDXProp_ObjectTag_Type:
		m_Type=(CDXObjectTagType)src_arg.GetUINT16();
		break;

	case kCDXProp_ObjectTag_Tracking:
		m_Tracking=src_arg.GetUINT8();
		break;

	case kCDXProp_ObjectTag_Persistent:
		m_Persistent=src_arg.GetUINT8();
		break;

	case kCDXProp_ObjectTag_Value:
		switch (m_Type) {

		case kCDXObjectTagType_Double:
		{
			FLOAT64 f64;
			src_arg.GetBytes((char*)&f64, sizeof f64);
			m_DoubleVal = SwapBytes(f64);
		}
		break;

		case kCDXObjectTagType_Int32:
			m_Int32Val=src_arg.GetUINT32();
			break;

#if 0	// can delete this passage.   Temp stub only -heh 7/15/07
		case 4 /* formerly, kCDXObjectTagType_Annotation */:
			ASSERT (false);	// CDX should be updated
			// A few ChemFinder records contain this tag
#endif

		case kCDXObjectTagType_String:
		{
			CDXString cdxs;
			cdxs.Read(src_arg, size_arg);
			m_StringVal = cdxs.str();
		}
		break;

		default:  // Ignore any ObjectTag values with unknown or undefined type.
			ASSERT (false);	// shouldn't occur
			CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
		}
	break;	

	case kCDXProp_Name:
		m_Name = src_arg.GetString(size_arg);
		break;
	case kCDXProp_Positioning:
		m_positioning = (CDXPositioningType) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_PositioningAngle:
		m_positioningAngle = src_arg.GetINT32();
		break;
	case kCDXProp_PositioningOffset:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_positioningOffset = CDXPoint2D(x, y);
		break;
		}
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}


// Read CDXML and create a CDX structure	
void CDXObjectTag::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_ObjectTag_Type:
		m_Type = sXMLObjectTagType.lookup(attribValue_arg);
	break;

	case kCDXProp_ObjectTag_Tracking:
		m_Tracking = attribValue_arg == "yes";
	break;

	case kCDXProp_ObjectTag_Persistent:
		m_Persistent = attribValue_arg == "yes";
	break;

	case kCDXProp_Name:
		m_Name=attribValue_arg;
	break;

	case kCDXProp_ObjectTag_Value:
		switch (m_Type) {

		case kCDXObjectTagType_Double:
			m_DoubleVal = atof(attribValue_arg.c_str());
		break;

		case kCDXObjectTagType_Int32:
			m_Int32Val = atoi(attribValue_arg.c_str());
		break;

		case kCDXObjectTagType_String:
            m_StringVal = UnicodeToString(attribValue_arg);
		break;

		default:
		break;
		}
	break;

	case kCDXProp_Positioning:
		m_positioning = sXMLPositioningType.lookup(attribValue_arg);
		break;
	case kCDXProp_PositioningAngle:
		m_positioningAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
		break;
	case kCDXProp_PositioningOffset:
		m_positioningOffset = StringToCDXPoint2D(attribValue_arg);
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}

}

#if TARGET_OS_WIN32
#pragma warning( pop )
#endif

// Take a CDX structure and create a CDX output stream.
void CDXObjectTag::WriteAttributesTo(CDXDataSink &sink_arg) const {
	
	sink_arg.PutAttribute( kCDXProp_ObjectTag_Type, UINT16(m_Type));

	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_positioning != kCDXPositioningType_Auto)
		sink_arg.PutAttribute( kCDXProp_Positioning, (UINT8) m_positioning );

	if (m_positioningAngle != 0)
		sink_arg.PutAttribute( kCDXProp_PositioningAngle, m_positioningAngle );

	if (m_positioningOffset.x != 0 || m_positioningOffset.y != 0)
		sink_arg.PutAttribute( kCDXProp_PositioningOffset, m_positioningOffset );

	if (m_Tracking)
		sink_arg.PutAttribute( kCDXProp_ObjectTag_Tracking, UINT8(m_Tracking));

	if (!m_Persistent)
		sink_arg.PutAttribute( kCDXProp_ObjectTag_Persistent, UINT8(m_Persistent));

	sink_arg.PutAttributeString(kCDXProp_Name, m_Name);

	switch (m_Type) {

	case kCDXObjectTagType_Double:
		{
		FLOAT64 scratch8=SwapBytes(m_DoubleVal);
		sink_arg.PutAttribute(kCDXProp_ObjectTag_Value, FLOAT64(scratch8));
		break;
		}

	case kCDXObjectTagType_Int32:
		sink_arg.PutAttribute(kCDXProp_ObjectTag_Value, m_Int32Val);
		break;

	case kCDXObjectTagType_String:
        if (!m_StringVal.empty())
        {
            sink_arg.PutAttributeCDXString(kCDXProp_ObjectTag_Value, CDXString(m_StringVal));
        }
		break;

	default:
		break;
	}
}
	
	
// Take a CDX structure and write a CDXML output stream
void CDXObjectTag::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Type, m_Type);

	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_positioning != kCDXPositioningType_Auto)
		CDXMLPutAttribute(sink_arg, kCDXProp_Positioning, m_positioning );

	if (m_positioningAngle != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_PositioningAngle, CdxAngleToDegrees(m_positioningAngle) );

	if (m_positioningOffset.x != 0 || m_positioningOffset.y != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_PositioningOffset, m_positioningOffset );

	CDXMLPutAttribute(sink_arg, kCDXProp_Name, StringToUnicode("", m_Name));

	if (m_Tracking)
		CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Tracking, (bool)m_Tracking);

	if (!m_Persistent)
		CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Persistent, (bool)m_Persistent);

	switch (m_Type) {

	case kCDXObjectTagType_Double:
		CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Value, m_DoubleVal);
		break;

	case kCDXObjectTagType_Int32:
		CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Value, m_Int32Val);
		break;

	case kCDXObjectTagType_String:
        if (!m_StringVal.empty())
        {
            CDXMLPutAttribute(sink_arg, kCDXProp_ObjectTag_Value, StringToUnicode("", m_StringVal));
        }
		break;

	default:
		break;
	}
}
