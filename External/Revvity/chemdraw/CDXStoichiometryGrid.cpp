// CommonCS/LibCommon/Src/CDXStoichiometryGrid.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 2001-2004, CambridgeSoft Corp., All Rights Reserved

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
#include <sstream>
#include "cs_stringUtils.h"

using namespace cs;

class CDXStoichiometryGrid;
class CDXSGComponent;
class CDXSGDatum;

// ***********************
// **  class CDXStoichiometryGrid   **
// ***********************
//
// Specialization of CDXObject for CDXStoichiometryGrid objects

CDXStoichiometryGrid::CDXStoichiometryGrid(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_StoichiometryGrid, id_arg)
{
}

CDXStoichiometryGrid::CDXStoichiometryGrid (const CDXStoichiometryGrid &src)
	:	CDXGraphicObject(src)
{
}

CDXObject*	CDXStoichiometryGrid::Clone() const
{
	return new CDXStoichiometryGrid(*this);
}

std::string CDXStoichiometryGrid::XMLObjectName() const
{
	return kCDXML_stoichiometrygrid;
}

void CDXStoichiometryGrid::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
}

void CDXStoichiometryGrid::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
}

void CDXStoichiometryGrid::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

}

void CDXStoichiometryGrid::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}



// **************************
// **  class CDXSGComponent **
// **************************
//
// Specialization of CDXObject for CDXSGComponent objects

CDXSGComponent::CDXSGComponent(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_SGComponent, id_arg)
	,	m_width(0)
	,	m_isReactant(false)
	,	m_isHeader(false)
	,   m_referenceID(-1)
{
}

CDXSGComponent::CDXSGComponent (const CDXSGComponent &src)
	:	CDXGraphicObject(src)
	,	m_width(src.m_width)
	,	m_isReactant(src.m_isReactant)
	,	m_isHeader(src.m_isHeader)
	,   m_referenceID(src.m_referenceID)
{
}

CDXObject*	CDXSGComponent::Clone() const
{
	return new CDXSGComponent(*this);
}

std::string CDXSGComponent::XMLObjectName() const
{
	return kCDXML_sgcomponent;
}

void CDXSGComponent::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_SG_ComponentReferenceID:
		m_referenceID = src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_SG_ComponentIsReactant:
		m_isReactant = true;
		break;

	case kCDXProp_SG_ComponentIsHeader:
		m_isHeader = true;
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXSGComponent::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Width:
		m_width = StrToNum(attribValue_arg);
		break;

	case kCDXProp_SG_ComponentReferenceID:
		m_referenceID = StrToNum(attribValue_arg);
		break;

	case kCDXProp_SG_ComponentIsReactant:
		m_isReactant = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_SG_ComponentIsHeader:
		m_isHeader = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXSGComponent::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	sink_arg.PutAttribute( kCDXProp_Width, m_width );
	if (m_referenceID != -1)
		sink_arg.PutAttribute( kCDXProp_SG_ComponentReferenceID, m_referenceID );
	if (m_isReactant)
		sink_arg.PutAttribute( kCDXProp_SG_ComponentIsReactant );
	if (m_isHeader)
		sink_arg.PutAttribute( kCDXProp_SG_ComponentIsHeader );

	CDXGraphicObject::WriteAttributesTo(sink_arg);

}

void CDXSGComponent::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXMLPutAttribute( sink_arg, kCDXProp_Width, m_width );
	if (m_referenceID != -1)
		CDXMLPutAttribute( sink_arg, kCDXProp_SG_ComponentReferenceID, m_referenceID );
	if (m_isReactant)
		CDXMLPutAttribute( sink_arg, kCDXProp_SG_ComponentIsReactant );
	if (m_isHeader)
		CDXMLPutAttribute( sink_arg, kCDXProp_SG_ComponentIsHeader );

	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}





// ***********************
// **  class CDXSGDatum **
// ***********************
//
// Specialization of CDXObject for CDXSGDatum objects

CDXSGDatum::CDXSGDatum(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_SGDatum, id_arg)
	,	m_type(kCDXSGPropertyType_Unknown)
	,	m_datatype(kCDXSGDataType_Unspecified)
	,	m_propertystring("")
	,	m_propertyvalue(0)
	,	m_isreadonly(false)
	,	m_isedited(false)
	,	m_ishidden(false)
{
}

CDXSGDatum::CDXSGDatum (const CDXSGDatum &src)
	:	CDXGraphicObject(src)
	,	m_type(src.m_type)
	,	m_datatype(src.m_datatype)
	,	m_propertystring(src.m_propertystring)
	,	m_propertyvalue(src.m_propertyvalue)
	,	m_isreadonly(src.m_isreadonly)
	,	m_isedited(src.m_isedited)
	,	m_ishidden(src.m_ishidden)
{
}

CDXObject*	CDXSGDatum::Clone() const
{
	return new CDXSGDatum(*this);
}

std::string CDXSGDatum::XMLObjectName() const
{
	return kCDXML_sgdatum;
}

void CDXSGDatum::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch (attribTag_arg)
	{
	case kCDXProp_SG_DataType:
		m_datatype = (CDXSGDataType)src_arg.GetUINT16();
		break;

	case kCDXProp_SG_PropertyType:
		m_type = (CDXSGPropertyType)src_arg.GetUINT16();
		break;

	case kCDXProp_SG_DataValue:	
		if (m_datatype != kCDXSGDataType_String)
		{
			src_arg.GetBytes((char *)&f64, sizeof(f64));
			m_propertyvalue = SwapBytes(f64);
		}
		else
		{
			m_propertystring = src_arg.GetString(size_arg);
		}
		break;

	case kCDXProp_IsReadOnly:
		m_isreadonly = true;
		break;

	case kCDXProp_IsEdited:
		m_isedited = true;
		break;

	case kCDXProp_IsHidden:
		m_ishidden = true;
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXSGDatum::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_SG_DataType:
		m_datatype = (CDXSGDataType)StrToNum(attribValue_arg);
		break;

	case kCDXProp_SG_PropertyType:
		m_type = (CDXSGPropertyType)StrToNum(attribValue_arg);
		break;

	case kCDXProp_SG_DataValue:	
		if (m_datatype != kCDXSGDataType_String)
			m_propertyvalue = StrToDub(attribValue_arg);
		else
			m_propertystring = attribValue_arg;
		break;

	case kCDXProp_IsReadOnly:
		m_isreadonly = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_IsEdited:
		m_isedited = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_IsHidden:
		m_ishidden = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXSGDatum::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	if (m_datatype != kCDXSGDataType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_SG_DataType, (UINT16)m_datatype );

	if (m_datatype != kCDXSGDataType_String)
		sink_arg.PutAttribute( kCDXProp_SG_DataValue, FLOAT64(m_propertyvalue));
	else
		sink_arg.PutAttributeString( kCDXProp_SG_DataValue, m_propertystring);

	if (m_type != kCDXSGPropertyType_Unknown)
		sink_arg.PutAttribute( kCDXProp_SG_PropertyType, (UINT16)m_type );

	if (m_isreadonly)
		sink_arg.PutAttribute( kCDXProp_IsReadOnly );

	if (m_isedited)
		sink_arg.PutAttribute( kCDXProp_IsEdited );

	if (m_ishidden)
		sink_arg.PutAttribute( kCDXProp_IsHidden );

	CDXGraphicObject::WriteAttributesTo(sink_arg);


}

void CDXSGDatum::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	if (m_datatype != kCDXSGDataType_Unspecified)
		CDXMLPutAttribute( sink_arg, kCDXProp_SG_DataType, m_datatype );

	if (m_datatype != kCDXSGDataType_String)
		CDXMLPutAttribute(sink_arg, kCDXProp_SG_DataValue, m_propertyvalue);
	else
		CDXMLPutAttribute(sink_arg, kCDXProp_SG_DataValue, m_propertystring);

	if (m_type != kCDXSGPropertyType_Unknown)
		CDXMLPutAttribute( sink_arg, kCDXProp_SG_PropertyType, m_type );

	if (m_isreadonly)
		CDXMLPutAttribute( sink_arg, kCDXProp_IsReadOnly, m_isreadonly );

	if (m_isedited)
		CDXMLPutAttribute( sink_arg, kCDXProp_IsEdited, m_isedited );

	if (m_ishidden)
		CDXMLPutAttribute( sink_arg, kCDXProp_IsHidden, m_ishidden );

	CDXGraphicObject::XMLWriteAttributes(sink_arg);
}



