// CommonCS/LibCommon/Src/CDXBioDraw.cpp
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
#include <stdlib.h> // for atoi
#include <sstream>

/*
Explanation of the attributes of the CDXBioShape object

CDXBioShapeType m_type;
	it is an integer value describing the general look of a BioShape object.
	Supported values are:
		
	kCDXBioShape1SubstrateEnzyme
	kCDXBioShape2SubstrateEnzyme
	kCDXBioShapeReceptor
	kCDXBioShapeGProteinAlpha
	kCDXBioShapeGProteinBeta
	kCDXBioShapeGProteinGamma
	kCDXBioShapeImmunoglobin
	kCDXBioShapeIonChannel
	kCDXBioShapeEndoplasmicReticulum
	kCDXBioShapeGolgi
	kCDXBioShapeMitochondrion
	kCDXBioShapeCloud

	kCDXBioShapeMembraneLine
	kCDXBioShapeMembraneArc
	kCDXBioShapeMembraneEllipse
	kCDXBioShapeDNA
	kCDXBioShapeHelixProtein
	kCDXBioShapeGolgi1Element
	kCDXBioShapeGolgi2Element

	The first group is transformed to the target coordinate system
	only by using the Position() and m_majorAxisEnd parameter. the second group is 
	transformed by using Position(), m_majorAxisEnd and m_minorAxisEnd parameters.

CDXPoint2D Position();
	It is the TopLeft point of the bounding rectangle for the first group, and
	the first pivot point for the second group

CDXPoint2D m_majorAxisEnd;
	It is the BottomRight point of the bounding rectangle for the first group, and
	the second pivot point for the second group

CDXPoint2D m_minorAxisEnd;
	It is used only for the first group of bioshapes. It is the TopRight point of the
	bounding rectangle at creation, and the if the object is rotated, than the rotated
	TopRight point.

vector<double> m_params;
	the number and meaning of the parameters depends on the value of m_type attribute:

kCDXBioShape1SubstrateEnzyme:
	0. height
	1. width
	The proportion of the two parameters define the look of the shape

kCDXBioShapeReceptor:
	0. radius of bottom part
	1. radius of middle part
	2. height of bottom part
	3. height of the neck
	4. height of top part up to the lowest point of surface
	5. height of the top part from the lowest point of surface to the end

kCDXBioShapeGProteinGamma:
	0. height of upper part
	1. height of lower part
	2. width

kCDXBioShapeHelixProtein:
	0. diameter of the cylinder
	1. height of the cylinder
	2. distance in between cylinders
	3. diameter of the pipe connecting cylinders

kCDXBioShapeMembraneArc:
	0. size of the circle
	1. a param of ellipse
	2. b param of ellipse
	3. start angle
	3. end angle

kCDXBioShapeMembraneLine:
	0. size of the circle

kCDXBioShapeMembraneEllipse:
	0. size of the circle
	1. a param of ellipse
	2. b param of ellipse
	3. start angle
	3. end angle

kCDXBioShapeGolgi1Element:
	0. height
	1. width of the shape
	2. length of the shape

kCDXBioShapeGolgi2Element:
	0. height
	1. width of the shape
	2. length of the shape

*/

XMLPutEnum<CDXBioShapeType>::Value sXMLBioShapeTypeValues[] = {
 	{kCDXBioShapeUndefined,				"Undefined"},				
 	{kCDXBioShape1SubstrateEnzyme,		"1SubstrateEnzyme"},			
 	{kCDXBioShape2SubstrateEnzyme,		"2SubstrateEnzyme"},			
 	{kCDXBioShapeReceptor,				"Receptor"},				
 	{kCDXBioShapeGProteinAlpha,			"GProteinAlpha"},				
 	{kCDXBioShapeGProteinBeta,			"GProteinBeta"},				
 	{kCDXBioShapeGProteinGamma,			"GProteinGamma"},				
 	{kCDXBioShapeImmunoglobin,			"Immunoglobin"},			
 	{kCDXBioShapeIonChannel,			"IonChannel"},			
 	{kCDXBioShapeEndoplasmicReticulum,	"EndoplasmicReticulum"},		
 	{kCDXBioShapeGolgi,					"Golgi"},					
 	{kCDXBioShapeMembraneLine,			"MembraneLine"},			
 	{kCDXBioShapeMembraneArc,			"MembraneArc"},				
 	{kCDXBioShapeMembraneEllipse,		"MembraneEllipse"},			
	{kCDXBioShapeMembraneMicelle,		"MembraneMicelle"},			
 	{kCDXBioShapeDNA,					"DNA"},					
 	{kCDXBioShapeHelixProtein,			"HelixProtein"},				
	{kCDXBioShapeMitochondrion,			"Mitochondrion"},				
	{kCDXBioShapeCloud,					"Cloud"},			
	{kCDXBioShapetRNA,					"tRNA"},		
	{kCDXBioShapeRibosomeA,				"RibosomeA"},			
	{kCDXBioShapeRibosomeB,				"RibosomeB"}			
};

XMLPutEnum<CDXBioShapeType> sXMLBioShapeType(sXMLBioShapeTypeValues, sizeof sXMLBioShapeTypeValues, kCDXBioShapeUndefined);
void XMLPut(XMLDataSink &sink, CDXBioShapeType  v) {	sink.os << sXMLBioShapeType.lookup(v); }

CDXBioShape::CDXBioShape(CDXObjectID id_arg)
	: CDXGraphicObject(kCDXObj_BioShape, id_arg)
	, m_type(kCDXBioShapeUndefined)
	, m_fadePercent(1000)
	, m_1SubstrateEnzyme_ReceptorSize(0)
	, m_receptor_NeckWidth(0)
	, m_gprotein_UpperHeight(0)
	, m_membrane_ElementSize(0)
	, m_membrane_StartAngle(0)
	, m_membrane_EndAngle(0)
	, m_dna_WaveHeight(0)
	, m_dna_WaveLength(0)
	, m_dna_WaveWidth(0)
	, m_dna_Offset(0)
	, m_helixProtein_CylinderWidth(0)
	, m_helixProtein_CylinderHeight(0)
	, m_helixProtein_CylinderDistance(0)
	, m_helixProtein_PipeWidth(0)
	, m_helixProtein_Extra(0)
	, m_flags(0)
{
}

CDXBioShape::CDXBioShape(const CDXBioShape &src)
	: CDXGraphicObject(src)	
	, m_type(src.m_type)
	, m_majorAxisEnd(src.m_majorAxisEnd)
	, m_minorAxisEnd(src.m_minorAxisEnd)
	, m_fadePercent(src.m_fadePercent)
	, m_1SubstrateEnzyme_ReceptorSize(src.m_1SubstrateEnzyme_ReceptorSize)
	, m_receptor_NeckWidth(src.m_receptor_NeckWidth)
	, m_gprotein_UpperHeight(src.m_gprotein_UpperHeight)
	, m_membrane_ElementSize(src.m_membrane_ElementSize)
	, m_membrane_StartAngle(src.m_membrane_StartAngle)
	, m_membrane_EndAngle(src.m_membrane_EndAngle)
	, m_dna_WaveHeight(src.m_dna_WaveHeight)
	, m_dna_WaveLength(src.m_dna_WaveLength)
	, m_dna_WaveWidth(src.m_dna_WaveWidth)
	, m_dna_Offset(src.m_dna_Offset)
	, m_helixProtein_CylinderWidth(src.m_helixProtein_CylinderWidth)
	, m_helixProtein_CylinderHeight(src.m_helixProtein_CylinderHeight)
	, m_helixProtein_CylinderDistance(src.m_helixProtein_CylinderDistance)
	, m_helixProtein_PipeWidth(src.m_helixProtein_PipeWidth)
	, m_helixProtein_Extra(src.m_helixProtein_Extra)
	, m_flags(src.m_flags)

{
}

void CDXBioShape::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_BioShape_Type:
		m_type = (CDXBioShapeType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_3DMajorAxisEnd:
		m_majorAxisEnd.x = src_arg.GetUINT32();
		m_majorAxisEnd.y = src_arg.GetUINT32();
		m_majorAxisEnd.z = src_arg.GetUINT32();
		Known(has_majorAxisEnd, true);
		break;
	case kCDXProp_3DMinorAxisEnd:
		m_minorAxisEnd.x = src_arg.GetUINT32();
		m_minorAxisEnd.y = src_arg.GetUINT32();
		m_minorAxisEnd.z = src_arg.GetUINT32();
		Known(has_minorAxisEnd, true);
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = src_arg.GetUINT(size_arg);
		Known(CDXBioShape::has_fadePercent, true);
		break;
	case kCDXProp_1SubstrateEnzyme_ReceptorSize:
		m_1SubstrateEnzyme_ReceptorSize = src_arg.GetUINT(size_arg);
		Known(has_1SubstrateEnzyme_ReceptorSize, true);
		break;
	case kCDXProp_Receptor_NeckWidth:
		m_receptor_NeckWidth = src_arg.GetUINT(size_arg);
		Known(has_receptor_NeckWidth, true);
		break;
	case kCDXProp_Gprotein_UpperHeight:
		m_gprotein_UpperHeight = src_arg.GetUINT(size_arg);
		Known(has_gprotein_UpperHeight, true);
		break;
	case kCDXProp_Membrane_ElementSize:
		m_membrane_ElementSize = src_arg.GetUINT(size_arg);
		Known(has_membrane_ElementSize, true);
		break;
	case kCDXProp_Membrane_StartAngle:
		m_membrane_StartAngle = src_arg.GetUINT(size_arg);
		Known(has_membrane_StartAngle, true);
		break;
	case kCDXProp_Membrane_EndAngle:
		m_membrane_EndAngle = src_arg.GetUINT(size_arg);
		Known(has_membrane_EndAngle, true);
		break;
	case kCDXProp_DNA_WaveHeight:
		m_dna_WaveHeight = src_arg.GetUINT(size_arg);
		Known(has_dna_WaveHeight, true);
		break;
	case kCDXProp_DNA_WaveLength:
		m_dna_WaveLength = src_arg.GetUINT(size_arg);
		Known(has_dna_WaveLength, true);
		break;
	case kCDXProp_DNA_WaveWidth:
		m_dna_WaveWidth = src_arg.GetUINT(size_arg);
		Known(has_dna_WaveWidth, true);
		break;
	case kCDXProp_DNA_Offset:
		m_dna_Offset = src_arg.GetUINT(size_arg);
		Known(has_dna_Offset, true);
		break;
	case kCDXProp_HelixProtein_CylinderWidth:
		m_helixProtein_CylinderWidth = src_arg.GetUINT(size_arg);
		Known(has_helixProtein_CylinderWidth, true);
		break;
	case kCDXProp_HelixProtein_CylinderHeight:
		m_helixProtein_CylinderHeight = src_arg.GetUINT(size_arg);
		Known(has_helixProtein_CylinderHeight, true);
		break;
	case kCDXProp_HelixProtein_CylinderDistance:
		m_helixProtein_CylinderDistance = src_arg.GetUINT(size_arg);
		Known(has_helixProtein_CylinderDistance, true);
		break;
	case kCDXProp_HelixProtein_PipeWidth:
		m_helixProtein_PipeWidth = src_arg.GetUINT(size_arg);
		Known(has_helixProtein_PipeWidth, true);
		break;
	case kCDXProp_HelixProtein_Extra:
		m_helixProtein_Extra = src_arg.GetUINT(size_arg);
		Known(has_helixProtein_Extra, true);
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}
void CDXBioShape::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_BioShape_Type:
		m_type = sXMLBioShapeType.lookup(attribValue_arg);
		break;
	case kCDXProp_3DMajorAxisEnd:
		m_majorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		Known(has_majorAxisEnd, true);
		break;
	case kCDXProp_3DMinorAxisEnd:
		m_minorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		Known(has_minorAxisEnd, true);
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = atoi(attribValue_arg.c_str());
		Known(has_fadePercent, true);
		break;
	case kCDXProp_1SubstrateEnzyme_ReceptorSize:
		m_1SubstrateEnzyme_ReceptorSize = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_1SubstrateEnzyme_ReceptorSize, true);
		break;
	case kCDXProp_Receptor_NeckWidth:
		m_receptor_NeckWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_receptor_NeckWidth, true);
		break;
	case kCDXProp_Gprotein_UpperHeight:
		m_gprotein_UpperHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_gprotein_UpperHeight, true);
		break;
	case kCDXProp_Membrane_ElementSize:
		m_membrane_ElementSize = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_membrane_ElementSize, true);
		break;
	case kCDXProp_Membrane_StartAngle:
		m_membrane_StartAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
		Known(has_membrane_StartAngle, true);
		break;
	case kCDXProp_Membrane_EndAngle:
		m_membrane_EndAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
		Known(has_membrane_EndAngle, true);
		break;
	case kCDXProp_DNA_WaveHeight:
		m_dna_WaveHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_dna_WaveHeight, true);
		break;
	case kCDXProp_DNA_WaveLength:
		m_dna_WaveLength = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_dna_WaveLength, true);
		break;
	case kCDXProp_DNA_WaveWidth:
		m_dna_WaveWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_dna_WaveWidth, true);
		break;
	case kCDXProp_DNA_Offset:
		m_dna_Offset = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_dna_Offset, true);
		break;
	case kCDXProp_HelixProtein_CylinderWidth:
		m_helixProtein_CylinderWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_helixProtein_CylinderWidth, true);
		break;
	case kCDXProp_HelixProtein_CylinderHeight:
		m_helixProtein_CylinderHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_helixProtein_CylinderHeight, true);
		break;
	case kCDXProp_HelixProtein_CylinderDistance:
		m_helixProtein_CylinderDistance = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_helixProtein_CylinderDistance, true);
		break;
	case kCDXProp_HelixProtein_PipeWidth:
		m_helixProtein_PipeWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_helixProtein_PipeWidth, true);
		break;
	case kCDXProp_HelixProtein_Extra:
		m_helixProtein_Extra = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		Known(has_helixProtein_Extra, true);
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXBioShape::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_type != kCDXBioShapeUndefined)
		sink_arg.PutAttribute(kCDXProp_BioShape_Type, (INT16)m_type);

	if (Known(has_majorAxisEnd))
		sink_arg.PutAttribute(kCDXProp_3DMajorAxisEnd, m_majorAxisEnd);

	if (Known(has_minorAxisEnd))
		sink_arg.PutAttribute(kCDXProp_3DMinorAxisEnd, m_minorAxisEnd);

	if (Known(has_fadePercent))
		sink_arg.PutAttribute( kCDXProp_FadePercent, m_fadePercent );

	if (Known(has_1SubstrateEnzyme_ReceptorSize))
		sink_arg.PutAttribute(kCDXProp_1SubstrateEnzyme_ReceptorSize, m_1SubstrateEnzyme_ReceptorSize);

	if (Known(has_receptor_NeckWidth))
		sink_arg.PutAttribute(kCDXProp_Receptor_NeckWidth, m_receptor_NeckWidth);

	if (Known(has_gprotein_UpperHeight))
		sink_arg.PutAttribute(kCDXProp_Gprotein_UpperHeight, m_gprotein_UpperHeight);

	if (Known(has_membrane_ElementSize))
		sink_arg.PutAttribute(kCDXProp_Membrane_ElementSize, m_membrane_ElementSize);

	if (Known(has_membrane_StartAngle))
		sink_arg.PutAttribute(kCDXProp_Membrane_StartAngle, m_membrane_StartAngle);

	if (Known(has_membrane_EndAngle))
		sink_arg.PutAttribute(kCDXProp_Membrane_EndAngle, m_membrane_EndAngle);

	if (Known(has_dna_WaveHeight))
		sink_arg.PutAttribute(kCDXProp_DNA_WaveHeight, m_dna_WaveHeight);

	if (Known(has_dna_WaveLength))
		sink_arg.PutAttribute(kCDXProp_DNA_WaveLength, m_dna_WaveLength);

	if (Known(has_dna_WaveWidth))
		sink_arg.PutAttribute(kCDXProp_DNA_WaveWidth, m_dna_WaveWidth);

	if (Known(has_dna_Offset))
		sink_arg.PutAttribute(kCDXProp_DNA_Offset, m_dna_Offset);

	if (Known(has_helixProtein_CylinderWidth))
		sink_arg.PutAttribute(kCDXProp_HelixProtein_CylinderWidth, m_helixProtein_CylinderWidth);

	if (Known(has_helixProtein_CylinderHeight))
		sink_arg.PutAttribute(kCDXProp_HelixProtein_CylinderHeight, m_helixProtein_CylinderHeight);

	if (Known(has_helixProtein_CylinderDistance))
		sink_arg.PutAttribute(kCDXProp_HelixProtein_CylinderDistance, m_helixProtein_CylinderDistance);

	if (Known(has_helixProtein_PipeWidth))
		sink_arg.PutAttribute(kCDXProp_HelixProtein_PipeWidth, m_helixProtein_PipeWidth);

	if (Known(has_helixProtein_Extra))
		sink_arg.PutAttribute(kCDXProp_HelixProtein_Extra, m_helixProtein_Extra);

}

void CDXBioShape::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_type != kCDXBioShapeUndefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_BioShape_Type, m_type);

	if (Known(has_majorAxisEnd))
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMajorAxisEnd, m_majorAxisEnd);

	if (Known(has_minorAxisEnd))
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMinorAxisEnd, m_minorAxisEnd);

	if (Known(has_fadePercent))
		CDXMLPutAttribute(sink_arg, kCDXProp_FadePercent, m_fadePercent );

	if (Known(has_1SubstrateEnzyme_ReceptorSize))
		CDXMLPutAttribute(sink_arg, kCDXProp_1SubstrateEnzyme_ReceptorSize, CDXCoordinateToString(m_1SubstrateEnzyme_ReceptorSize));

	if (Known(has_receptor_NeckWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_Receptor_NeckWidth, CDXCoordinateToString(m_receptor_NeckWidth));

	if (Known(has_gprotein_UpperHeight))
		CDXMLPutAttribute(sink_arg, kCDXProp_Gprotein_UpperHeight, CDXCoordinateToString(m_gprotein_UpperHeight));

	if (Known(has_membrane_ElementSize))
		CDXMLPutAttribute(sink_arg, kCDXProp_Membrane_ElementSize, CDXCoordinateToString(m_membrane_ElementSize));

	if (Known(has_membrane_StartAngle))
		CDXMLPutAttribute(sink_arg, kCDXProp_Membrane_StartAngle, CdxAngleToDegrees(m_membrane_StartAngle));

	if (Known(has_membrane_EndAngle))
		CDXMLPutAttribute(sink_arg, kCDXProp_Membrane_EndAngle, CdxAngleToDegrees(m_membrane_EndAngle));

	if (Known(has_dna_WaveHeight))
		CDXMLPutAttribute(sink_arg, kCDXProp_DNA_WaveHeight, CDXCoordinateToString(m_dna_WaveHeight));

	if (Known(has_dna_WaveLength))
		CDXMLPutAttribute(sink_arg, kCDXProp_DNA_WaveLength, CDXCoordinateToString(m_dna_WaveLength));

	if (Known(has_dna_WaveWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_DNA_WaveWidth, CDXCoordinateToString(m_dna_WaveWidth));

	if (Known(has_dna_Offset))
		CDXMLPutAttribute(sink_arg, kCDXProp_DNA_Offset, CDXCoordinateToString(m_dna_Offset));

	if (Known(has_helixProtein_CylinderWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_HelixProtein_CylinderWidth, CDXCoordinateToString(m_helixProtein_CylinderWidth));

	if (Known(has_helixProtein_CylinderHeight))
		CDXMLPutAttribute(sink_arg, kCDXProp_HelixProtein_CylinderHeight, CDXCoordinateToString(m_helixProtein_CylinderHeight));

	if (Known(has_helixProtein_CylinderDistance))
		CDXMLPutAttribute(sink_arg, kCDXProp_HelixProtein_CylinderDistance, CDXCoordinateToString(m_helixProtein_CylinderDistance));

	if (Known(has_helixProtein_PipeWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_HelixProtein_PipeWidth, CDXCoordinateToString(m_helixProtein_PipeWidth));

	if (Known(has_helixProtein_Extra))
		CDXMLPutAttribute(sink_arg, kCDXProp_HelixProtein_Extra, CDXCoordinateToString(m_helixProtein_Extra));
}

std::string CDXBioShape::XMLObjectName() const
{
	return kCDXML_bioshape;
}

CDXObject* CDXBioShape::Clone() const
{
	return new CDXBioShape (*this);
}

CDXBioShape::~CDXBioShape()
{
}
