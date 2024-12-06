// CommonCS/LibCommon/Src/CDXTLCSpot.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"


void CDXTLCSpot::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch(attribTag_arg)
	{
	case kCDXProp_TLC_Rf:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_value = SwapBytes(f64);
		break;

	case kCDXProp_TLC_ShowRf:
		m_showValue = true;
		break;

	case kCDXProp_TLC_Tail:
		m_tail = src_arg.GetUINT(size_arg);
		break;

	default:
		CDXPlateItemBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
	}
}

std::string CDXTLCSpot::XMLObjectName() const
{
	return kCDXML_tlcspot;
}

CDXObject*	CDXTLCSpot::Clone() const
{
	return new CDXTLCSpot(*this);
}

CDXTLCSpot::CDXTLCSpot(CDXObjectID id)
	:	CDXPlateItemBase(id)
	,	m_tail(0) 
{
}

CDXTLCSpot::CDXTLCSpot(const CDXTLCSpot &a)
	:	CDXPlateItemBase(a)
	,	m_tail(a.m_tail) 
{
}

void CDXTLCSpot::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_TLC_Rf:
		m_value = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_TLC_ShowRf:
		m_showValue = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_TLC_Tail:
		m_tail = atoi(attribValue_arg.c_str());
		break;

	default:
		CDXPlateItemBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
	}
}

void CDXTLCSpot::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	if (m_value >= 0)
	{
		sink_arg.PutAttribute(kCDXProp_TLC_Rf, FLOAT64(m_value));
	}
	
	if (m_showValue)
	{
		sink_arg.PutAttribute(kCDXProp_TLC_ShowRf);
	}

	sink_arg.PutAttribute(kCDXProp_TLC_Tail, m_tail);
	CDXPlateItemBase::WriteAttributesTo(sink_arg);
}

void CDXTLCSpot::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	if (m_value >= 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_Rf, m_value);
	}
	
	if (m_showValue)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_TLC_ShowRf);
	}

	CDXMLPutAttribute(sink_arg, kCDXProp_TLC_Tail, m_tail);
	CDXPlateItemBase::XMLWriteAttributes(sink_arg);
}