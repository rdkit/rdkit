// CommonCS/LibCommon/Src/CDXTLCPlate.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"


CDXObject*	CDXTLCPlate::Clone() const
{
	return new CDXTLCPlate(*this);
}

std::string CDXTLCPlate::XMLObjectName() const
{
	return kCDXML_tlcplate;
}

CDXTLCPlate::CDXTLCPlate(CDXObjectID id)
	:	CDXPlateBase(id) 
	,	m_originFraction(-1)
	,	m_solventFrontFraction(-1)
	,	m_showOrigin(false)
	,	m_showSolventFront(false)
	,	m_showSideTicks(false)
{
}

CDXTLCPlate::CDXTLCPlate(const CDXTLCPlate &src) 
	:	CDXPlateBase(src)
	,	m_originFraction(src.m_originFraction)
	,	m_solventFrontFraction(src.m_solventFrontFraction)
	,	m_showOrigin(src.m_showOrigin)
	,	m_showSolventFront(src.m_showSolventFront)
	,	m_showSideTicks(src.m_showSideTicks)
{
}

void CDXTLCPlate::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch (attribTag_arg)
	{
	case kCDXProp_TLC_OriginFraction:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_originFraction = SwapBytes(f64);
		break;

	case kCDXProp_TLC_SolventFrontFraction:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_solventFrontFraction = SwapBytes(f64);
		break;

	case kCDXProp_TLC_ShowOrigin:
		m_showOrigin = true;
		break;

	case kCDXProp_TLC_ShowSolventFront:
		m_showSolventFront = true;
		break;

	case kCDXProp_TLC_ShowSideTicks:
		m_showSideTicks = true;
		break;

	default:
		CDXPlateBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
	}
}

void CDXTLCPlate::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_TLC_OriginFraction:
		m_originFraction = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_TLC_SolventFrontFraction:
		m_solventFrontFraction = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_TLC_ShowOrigin:
		m_showOrigin = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_TLC_ShowSolventFront:
		m_showSolventFront = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_TLC_ShowSideTicks:
		m_showSideTicks = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	default:
		CDXPlateBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
	}
}

void CDXTLCPlate::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	if (m_originFraction >= 0)
		sink_arg.PutAttribute( kCDXProp_TLC_OriginFraction, FLOAT64(m_originFraction));

	if (m_solventFrontFraction >= 0)
		sink_arg.PutAttribute( kCDXProp_TLC_SolventFrontFraction, FLOAT64(m_solventFrontFraction));

	if (m_showOrigin)
		sink_arg.PutAttribute(kCDXProp_TLC_ShowOrigin);

	if (m_showSolventFront)
		sink_arg.PutAttribute(kCDXProp_TLC_ShowSolventFront);

	if (m_showSideTicks)
		sink_arg.PutAttribute(kCDXProp_TLC_ShowSideTicks);

	CDXPlateBase::WriteAttributesTo(sink_arg);
}

void CDXTLCPlate::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	streamsize oldPrecision=sink_arg.os.precision(8);
	if (m_originFraction >= 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_OriginFraction, m_originFraction);

	if (m_solventFrontFraction >= 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_SolventFrontFraction, m_solventFrontFraction);
	sink_arg.os.precision(oldPrecision);

	if (m_showOrigin)
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_ShowOrigin);

	if (m_showSolventFront)
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_ShowSolventFront);

	if (m_showSideTicks)
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_ShowSideTicks);

	CDXPlateBase::XMLWriteAttributes(sink_arg);
}