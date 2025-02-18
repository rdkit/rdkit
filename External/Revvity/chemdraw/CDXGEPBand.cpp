// CommonCS/LibCommon/Src/CDXGEPBand.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2010, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXMLNames.h"

CDXGEPBand::CDXGEPBand(CDXObjectID id)
	: CDXPlateItemBase(id)
	, m_mass(0)
{
};

CDXGEPBand::CDXGEPBand(const CDXGEPBand &a)
	: CDXPlateItemBase(a)
	, m_mass(a.m_mass)
{
};

void CDXGEPBand::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	if (attribTag_arg == kCDXProp_GEP_Value)
	{
		m_mass = src_arg.GetINT32();
	}
	else if (attribTag_arg == kCDXProp_GEP_ShowValue)
	{
		m_showValue = true;
	}
	else
	{
		CDXPlateItemBase::StoreAttribute(src_arg, attribTag_arg, size_arg);
	}
}

std::string CDXGEPBand::XMLObjectName() const
{
	return kCDXML_gepband;
}

CDXObject*	CDXGEPBand::Clone() const
{
	return new CDXGEPBand(*this);
}

void CDXGEPBand::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	if (attribTag_arg == kCDXProp_GEP_Value)
	{
		m_mass = atoi(attribValue_arg.c_str());
	}
	else if (attribTag_arg == kCDXProp_GEP_ShowValue)
	{
		m_showValue = attribValue_arg == "yes" || attribValue_arg == "true";
	}
	else
	{
		CDXPlateItemBase::XMLStoreAttribute(attribTag_arg, attribValue_arg);
	}
}

void CDXGEPBand::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute(kCDXProp_GEP_Value, m_mass);
	
	if (m_showValue)
	{
		sink_arg.PutAttribute(kCDXProp_GEP_ShowValue);
	}

	CDXPlateItemBase::WriteAttributesTo(sink_arg);
}

void CDXGEPBand::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute(sink_arg, kCDXProp_GEP_Value, m_mass);
	
	if (m_showValue)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_GEP_ShowValue);
	}

	CDXPlateItemBase::XMLWriteAttributes(sink_arg);
}
