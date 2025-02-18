// CommonCS/LibCommon/Src/CDXPlateItemBase.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2010, CambridgeSoft Corp., All Rights Reserved

#include "CDXStdObjects.h"

// ***********************
// ** class CDXPlateItemBase  **
// ***********************
//
// Specialization of CDXObject for CDXPlateItemBase objects

CDXPlateItemBase::CDXPlateItemBase(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_TLCSpot, id_arg)
	,	m_value(0)
	,	m_width(0)
	,	m_height(0)
	,	m_displayType(0)
	,	m_showValue(false)
{
}

CDXPlateItemBase::CDXPlateItemBase (const CDXPlateItemBase &src)
	:	CDXGraphicObject(src)
	,	m_value(src.m_value)
	,	m_width(src.m_width)
	,	m_height(src.m_height)
	,	m_displayType(src.m_displayType)
	,	m_showValue(src.m_showValue)
{
}

void CDXPlateItemBase::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_Height:
		m_height = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_Curve_Type:
		m_displayType = src_arg.GetUINT(size_arg);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXPlateItemBase::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_Height:
		m_height = atoi(attribValue_arg.c_str());
		break;

	case kCDXProp_Curve_Type:
		m_displayType = atoi(attribValue_arg.c_str());
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXPlateItemBase::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_Width, m_width );
	sink_arg.PutAttribute( kCDXProp_Height, m_height );

	if (m_displayType != 0)
		sink_arg.PutAttribute( kCDXProp_Curve_Type, m_displayType );

	CDXGraphicObject::WriteAttributesTo(sink_arg);
}

void CDXPlateItemBase::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	streamsize oldPrecision=sink_arg.os.precision(8);
	sink_arg.os.precision(oldPrecision);

	CDXMLPutAttribute( sink_arg, kCDXProp_Width, m_width );
	CDXMLPutAttribute( sink_arg, kCDXProp_Height, m_height );

	if (m_displayType != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Curve_Type, m_displayType );

	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}