// CommonCS/LibCommon/Src/CDXGEPPlate.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

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
#include "CDXMLNames.h"

CDXGEPPlate::CDXGEPPlate(CDXObjectID id)
	: CDXPlateBase(id)
	,	m_scaleUnitID(0)
	,	m_minRange(-1)
	,	m_maxRange(-1)
	,	m_showScale(false)
	,	m_labelsAngle(0)
	,	m_axisWidth(0)
{
}

CDXGEPPlate::CDXGEPPlate(const CDXGEPPlate &a)
	:	CDXPlateBase(a)
	,	m_showScale(a.m_showScale) 
	,	m_labelsAngle(a.m_labelsAngle)
	,	m_scaleUnitID(a.m_scaleUnitID)
	,	m_minRange(a.m_minRange)
	,	m_maxRange(a.m_maxRange)
	,	m_axisWidth(a.m_axisWidth)
{
}

CDXObject*	CDXGEPPlate::Clone() const
{
	return new CDXGEPPlate(*this);
}

std::string CDXGEPPlate::XMLObjectName() const
{
	return kCDXML_gepplate;
}

void CDXGEPPlate::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXPlateBase::XMLWriteAttributes(sink_arg);
	if (m_scaleUnitID != 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_ScaleUnit, m_scaleUnitID);
	}

	if (m_minRange < m_maxRange)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_StartRange, m_minRange);
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_EndRange, m_maxRange);
	}

	if (m_showScale)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_ShowScale);
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_AxisWidth, m_axisWidth);
	}

	if (m_labelsAngle != 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_LaneLabelsAngle, m_labelsAngle);
	}
}

void CDXGEPPlate::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)	
	{
	case kCDXProp_GEP_ShowScale:
		m_showScale = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_GEP_ScaleUnit:
		m_scaleUnitID = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_GEP_StartRange:
		m_minRange = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_GEP_EndRange:
		m_maxRange = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_GEP_LaneLabelsAngle:
		m_labelsAngle = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_GEP_AxisWidth:
		m_axisWidth = atof(attribValue_arg.c_str());
		break;

	default: 
		CDXPlateBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXGEPPlate::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)	
	{
	case kCDXProp_GEP_ShowScale:
		m_showScale = true; // If it appears, then it's true
		break;

	case kCDXProp_GEP_ScaleUnit:
		m_scaleUnitID = src_arg.GetINT32();
		break;

	case kCDXProp_GEP_StartRange:
		m_minRange = src_arg.GetINT32();
		break;

	case kCDXProp_GEP_EndRange:
		m_maxRange = src_arg.GetINT32();
		break;

	case kCDXProp_GEP_LaneLabelsAngle:
		m_labelsAngle = src_arg.GetINT32();
		break;

	case kCDXProp_GEP_AxisWidth:
		m_axisWidth = src_arg.GetFLOAT64();
		break;

	default: 
		CDXPlateBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXGEPPlate::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXPlateBase::WriteAttributesTo(sink_arg);
	if (m_scaleUnitID != 0)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_ScaleUnit, (INT32)m_scaleUnitID);
	}

	if (m_minRange < m_maxRange)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_StartRange, (INT32)m_minRange);
		sink_arg.PutAttribute(kCDXProp_GEP_EndRange, (INT32)m_maxRange);
	}

	if (m_showScale)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_ShowScale);		
		sink_arg.PutAttribute(kCDXProp_GEP_AxisWidth, m_axisWidth);
	}

	if (m_labelsAngle != 0)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_LaneLabelsAngle, (INT32)m_labelsAngle);
	}
}