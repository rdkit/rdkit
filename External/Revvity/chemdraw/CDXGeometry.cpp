// CommonCS/LibCommon/Src/CDXGeometry.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

#include <ostream>

#pragma warning( disable : 4786 ) // identifier was truncated to 'number' characters in the debug information

#include "CDXStdObjects.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <sstream>
#include <stdexcept>

XMLPutEnum<CDXGeometricFeature>::Value sXMLGeometricFeatureValues[] = {
	{kCDXGeometricFeature_Undefined,						"Undefined"},
	{kCDXGeometricFeature_PointFromPointPointDistance,		"PointFromPointPointDistance"},
	{kCDXGeometricFeature_PointFromPointPointPercentage,	"PointFromPointPointPercentage"},
	{kCDXGeometricFeature_PointFromPointNormalDistance,		"PointFromPointNormalDistance"},
	{kCDXGeometricFeature_LineFromPoints,					"LineFromPoints"},
	{kCDXGeometricFeature_PlaneFromPoints,					"PlaneFromPoints"},
	{kCDXGeometricFeature_PlaneFromPointLine,				"PlaneFromPointLine"},
	{kCDXGeometricFeature_CentroidFromPoints,				"CentroidFromPoints"},
	{kCDXGeometricFeature_NormalFromPointPlane,				"NormalFromPointPlane"}
};

XMLPutEnum<CDXGeometricFeature> sXMLGeometricFeature(sXMLGeometricFeatureValues, sizeof sXMLGeometricFeatureValues, kCDXGeometricFeature_Undefined);
void XMLPut(XMLDataSink &sink, CDXGeometricFeature  v) {	sink.os << sXMLGeometricFeature.lookup(v); }

// *************************
// **  class CDXGeometry  **
// *************************
//
// Specialization of CDXObject for CDXGeometry objects

CDXGeometry::CDXGeometry(CDXObjectID id) 
: 	CDXGraphicObject(kCDXObj_Geometry, id), 
	m_geometricFeature(kCDXGeometricFeature_Undefined), 
	m_relationValue(0),
	m_pointIsDirected(false)
{
}

CDXGeometry::CDXGeometry (const CDXGeometry &src) 
:	CDXGraphicObject(src), 
	m_geometricFeature(src.m_geometricFeature), 
	m_relationValue(src.m_relationValue),
	m_pointIsDirected(src.m_pointIsDirected),
	m_basisObjects(src.m_basisObjects)
{
}


CDXObject* CDXGeometry::Clone() const
{
	return new CDXGeometry(*this);
}


std::string CDXGeometry::XMLObjectName() const
{
	return kCDXML_geometry;
}

void CDXGeometry::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch (attribTag_arg)
	{
	case kCDXProp_GeometricFeature:
		m_geometricFeature = (CDXGeometricFeature) src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_RelationValue:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_relationValue = SwapBytes(f64);
		break;

	case kCDXProp_PointIsDirected:
		m_pointIsDirected = size_arg == 0 || src_arg.GetUINT(size_arg) != 0;
		break;

	case kCDXProp_BasisObjects:
        m_basisObjects = ReadObjectIDList(src_arg, size_arg);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}

}

void CDXGeometry::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_GeometricFeature:
		m_geometricFeature = sXMLGeometricFeature.lookup(attribValue_arg);
		break;

	case kCDXProp_RelationValue:
		m_relationValue = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_PointIsDirected:
		m_pointIsDirected = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_BasisObjects:
		{
			std::istringstream is(attribValue_arg);
			m_basisObjects.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_basisObjects));
		}
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}


void CDXGeometry::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_geometricFeature != kCDXGeometricFeature_Undefined)
		sink_arg.PutAttribute( kCDXProp_GeometricFeature, (INT8) m_geometricFeature );
	
	if (m_relationValue != 0)
		sink_arg.PutAttribute( kCDXProp_RelationValue, FLOAT64(m_relationValue));

	if (m_pointIsDirected)
		sink_arg.PutAttribute( kCDXProp_PointIsDirected );

    sink_arg.PutAttributeForObjectIDList(kCDXProp_BasisObjects, m_basisObjects);
}


void CDXGeometry::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_geometricFeature != kCDXGeometricFeature_Undefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_GeometricFeature, m_geometricFeature );

	if (m_relationValue != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_RelationValue, m_relationValue);

	if (m_pointIsDirected)
		CDXMLPutAttribute(sink_arg, kCDXProp_PointIsDirected );

    CDXMLPutAttributeForObjectIDList(sink_arg, kCDXProp_BasisObjects, m_basisObjects);
}
