// CommonCS/LibCommon/Src/CDXGraphic.cpp
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
#include "CDXUtils.h"

XMLPutEnum<CDXGraphicType>::Value sXMLGraphicTypeValues[] = {
	{kCDXGraphicType_Undefined,	"Undefined"},
	{kCDXGraphicType_Line,		"Line"},
	{kCDXGraphicType_Arc,		"Arc"},
	{kCDXGraphicType_Rectangle,	"Rectangle"},
	{kCDXGraphicType_Oval,		"Oval"},
	{kCDXGraphicType_Orbital,	"Orbital"},
	{kCDXGraphicType_Bracket,	"Bracket"},
	{kCDXGraphicType_Symbol,	"Symbol"}
};

XMLPutEnum<CDXGraphicType> sXMLGraphicType(sXMLGraphicTypeValues, sizeof sXMLGraphicTypeValues, kCDXGraphicType_Undefined);
void XMLPut(XMLDataSink &sink, CDXGraphicType  v) {	sink.os << sXMLGraphicType.lookup(v); }

XMLPutEnum<CDXBracketType>::Value sXMLBracketTypeValues[] = {
	{kCDXBracketType_RoundPair,		"RoundPair"},
	{kCDXBracketType_SquarePair,	"SquarePair"},
	{kCDXBracketType_CurlyPair,		"CurlyPair"},
	{kCDXBracketType_Square,		"Square"},
	{kCDXBracketType_Curly,			"Curly"},
	{kCDXBracketType_Round,			"Round"}
};

XMLPutEnum<CDXBracketType> sXMLBracketType(sXMLBracketTypeValues, sizeof sXMLBracketTypeValues, kCDXBracketType_RoundPair);
void XMLPut(XMLDataSink &sink, CDXBracketType  v) {	sink.os << sXMLBracketType.lookup(v); }

XMLPutEnum<CDXRectangleType>::Value sXMLRectangleTypeValues[] = {
	{kCDXRectangleType_Plain,		"Plain"},
	{kCDXRectangleType_RoundEdge,	"RoundEdge"},
	{kCDXRectangleType_Shadow,		"Shadow"},
	{kCDXRectangleType_Shaded,		"Shaded"},
	{kCDXRectangleType_Filled,		"Filled"},
	{kCDXRectangleType_Dashed,		"Dashed"},
	{kCDXRectangleType_Bold,		"Bold"}
};

XMLPutEnum<CDXRectangleType> sXMLRectangleType(sXMLRectangleTypeValues, sizeof sXMLRectangleTypeValues, kCDXRectangleType_Plain);
void XMLPut(XMLDataSink &sink, CDXRectangleType  v) { sink.os << sXMLRectangleType.asTokens(v); }

XMLPutEnum<CDXOvalType>::Value sXMLOvalTypeValues[] = {
	{kCDXOvalType_Circle,	"Circle"},
	{kCDXOvalType_Shaded,	"Shaded"},
	{kCDXOvalType_Filled,	"Filled"},
	{kCDXOvalType_Dashed,	"Dashed"},
	{kCDXOvalType_Bold,		"Bold"},
	{kCDXOvalType_Shadowed,	"Shadowed"}
};

XMLPutEnum<CDXOvalType> sXMLOvalType(sXMLOvalTypeValues, sizeof sXMLOvalTypeValues, kCDXOvalType_Circle);
void XMLPut(XMLDataSink &sink, CDXOvalType  v)
{
	UINT16 intVal = UINT16(v);
	UINT16 i = 1;
	bool first = true;
	// Go throught the bits one at a time, removing them from the bitmask, until it's empty
	while (intVal != 0)
	{
		if ((intVal & i) != 0)
		{
			if (!first)
				sink.os << " ";
			first = false;
			sink.os << sXMLOvalType.lookup(CDXOvalType(i));
			intVal &= ~i;
		}
		i <<= 1;
	}
}

XMLPutEnum<CDXSymbolType>::Value sXMLSymbolTypeValues[] = {
	{kCDXSymbolType_LonePair,		"LonePair"},
	{kCDXSymbolType_Electron,		"Electron"},
	{kCDXSymbolType_RadicalCation,	"RadicalCation"},
	{kCDXSymbolType_RadicalAnion,	"RadicalAnion"},
	{kCDXSymbolType_CirclePlus,		"CirclePlus"},
	{kCDXSymbolType_CircleMinus,	"CircleMinus"},
	{kCDXSymbolType_Dagger,			"Dagger"},
	{kCDXSymbolType_DoubleDagger,	"DoubleDagger"},
	{kCDXSymbolType_Plus,			"Plus"},
	{kCDXSymbolType_Minus,			"Minus"},
	{kCDXSymbolType_Racemic,		"Racemic"},
	{kCDXSymbolType_Absolute,		"Absolute"},
	{kCDXSymbolType_Relative,		"Relative"},
	{kCDXSymbolType_LonePairBar,	"LonePairBar"},
};

XMLPutEnum<CDXSymbolType> sXMLSymbolType(sXMLSymbolTypeValues, sizeof sXMLSymbolTypeValues, kCDXSymbolType_LonePair);
void XMLPut(XMLDataSink &sink, CDXSymbolType  v) {	sink.os << sXMLSymbolType.lookup(v); }

//TODO:
// These are bit encoded - need to figure the valid combinations
XMLPutEnum<CDXArrowType>::Value sXMLArrowTypeValues[] = {
	{kCDXArrowType_NoHead,			"NoHead"},
	{kCDXArrowType_HalfHead,		"HalfHead"},
	{kCDXArrowType_FullHead,		"FullHead"},
	{kCDXArrowType_Resonance,		"Resonance"},
	{kCDXArrowType_Equilibrium,		"Equilibrium"},
	{kCDXArrowType_Hollow,			"Hollow"},
	{kCDXArrowType_RetroSynthetic,	"RetroSynthetic"}
};

XMLPutEnum<CDXArrowType> sXMLArrowType(sXMLArrowTypeValues, sizeof sXMLArrowTypeValues, kCDXArrowType_NoHead);
void XMLPut(XMLDataSink &sink, CDXArrowType  v) {	sink.os << sXMLArrowType.lookup(v); }

XMLPutEnum<CDXOrbitalType>::Value sXMLOrbitalTypeValues[] = {
	{kCDXOrbitalType_s,				"s"},
	{kCDXOrbitalType_oval,			"oval"},
	{kCDXOrbitalType_lobe,			"lobe"},
	{kCDXOrbitalType_p,				"p"},
	{kCDXOrbitalType_hybridPlus,	"hybridPlus"},
	{kCDXOrbitalType_hybridMinus,	"hybridMinus"},
	{kCDXOrbitalType_dz2Plus,		"dz2Plus"},
	{kCDXOrbitalType_dz2Minus,		"dz2Minus"},
	{kCDXOrbitalType_dxy,			"dxy"},
	{kCDXOrbitalType_sShaded,		"sShaded"},
	{kCDXOrbitalType_ovalShaded,		"ovalShaded"},
	{kCDXOrbitalType_lobeShaded,		"lobeShaded"},
	{kCDXOrbitalType_pShaded,			"pShaded"},
	{kCDXOrbitalType_sFilled,			"sFilled"},
	{kCDXOrbitalType_ovalFilled,		"ovalFilled"},
	{kCDXOrbitalType_lobeFilled,		"lobeFilled"},
	{kCDXOrbitalType_pFilled,			"pFilled"},
	{kCDXOrbitalType_hybridPlusFilled,	"hybridPlusFilled"},
	{kCDXOrbitalType_hybridMinusFilled,	"hybridMinusFilled"},
	{kCDXOrbitalType_dz2PlusFilled,		"dz2PlusFilled"},
	{kCDXOrbitalType_dz2MinusFilled,	"dz2MinusFilled"},
	{kCDXOrbitalType_dxyFilled,			"dxyFilled"}
};

XMLPutEnum<CDXOrbitalType> sXMLOrbitalType(sXMLOrbitalTypeValues, sizeof sXMLOrbitalTypeValues, kCDXOrbitalType_s);
void XMLPut(XMLDataSink &sink, CDXOrbitalType  v) {	sink.os << sXMLOrbitalType.lookup(v); }

XMLPutEnum<CDXFrameType>::Value sXMLFrameTypeValues[] = {
	{kCDXFrameType_Unspecified,	"Unspecified"},
	{kCDXFrameType_None,		"None"},
	{kCDXFrameType_Plain,		"Plain"},
	{kCDXFrameType_RoundEdge,	"RoundEdge"},
	{kCDXFrameType_Shadow,		"Shadow"},
	{kCDXFrameType_Shaded,		"Shaded"},
	{kCDXFrameType_Filled,		"Filled"},
	{kCDXFrameType_Dashed,		"Dashed"},
	{kCDXFrameType_Bold,		"Bold"}
};

XMLPutEnum<CDXFrameType> sXMLFrameType(sXMLFrameTypeValues, sizeof sXMLFrameTypeValues, kCDXFrameType_Unspecified);
void XMLPut(XMLDataSink &sink, CDXFrameType  v) { sink.os << sXMLFrameType.asTokens(v); }

XMLPutEnum<CDXBracketUsage>::Value sXMLBracketUsageValues[] = {
	{kCDXBracketUsage_Unspecified,			"Unspecified"},			
	{kCDXBracketUsage_Anypolymer,			"Anypolymer"},			
	{kCDXBracketUsage_Component,			"Component"},			
	{kCDXBracketUsage_Copolymer,			"Copolymer"},			
	{kCDXBracketUsage_CopolymerAlternating,	"CopolymerAlternating"},	
	{kCDXBracketUsage_CopolymerBlock,		"CopolymerBlock"},		
	{kCDXBracketUsage_CopolymerRandom,		"CopolymerRandom"},		
	{kCDXBracketUsage_Crosslink,			"Crosslink"},			
	{kCDXBracketUsage_Generic,				"Generic"},				
	{kCDXBracketUsage_Graft,				"Graft"},				
	{kCDXBracketUsage_Mer,					"Mer"},					
	{kCDXBracketUsage_MixtureOrdered,		"MixtureOrdered"},		
	{kCDXBracketUsage_MixtureUnordered,		"MixtureUnordered"},		
	{kCDXBracketUsage_Modification,			"Modification"},			
	{kCDXBracketUsage_Monomer,				"Monomer"},				
	{kCDXBracketUsage_MultipleGroup,		"MultipleGroup"},		
	{kCDXBracketUsage_SRU,					"SRU"}					
};

XMLPutEnum<CDXBracketUsage> sXMLBracketUsage(sXMLBracketUsageValues, sizeof sXMLBracketUsageValues, kCDXBracketUsage_Unspecified);
void XMLPut(XMLDataSink &sink, CDXBracketUsage  v) {	sink.os << sXMLBracketUsage.lookup(v); }

XMLPutEnum<CDXPolymerRepeatPattern>::Value sXMLPolymerRepeatPatternValues[] = {
	{kCDXPolymerRepeatPattern_HeadToTail,	"HeadToTail"},			
	{kCDXPolymerRepeatPattern_HeadToHead,	"HeadToHead"},			
	{kCDXPolymerRepeatPattern_EitherUnknown,"EitherUnknown"}					
};

XMLPutEnum<CDXPolymerRepeatPattern> sXMLPolymerRepeatPattern(sXMLPolymerRepeatPatternValues, sizeof sXMLPolymerRepeatPatternValues, kCDXPolymerRepeatPattern_HeadToTail);
void XMLPut(XMLDataSink &sink, CDXPolymerRepeatPattern  v) {	sink.os << sXMLPolymerRepeatPattern.lookup(v); }

XMLPutEnum<CDXPolymerFlipType>::Value sXMLPolymerFlipTypeValues[] = {
	{kCDXPolymerFlipType_Unspecified,		"Unspecified"},			
	{kCDXPolymerFlipType_NoFlip,			"NoFlip"},			
	{kCDXPolymerFlipType_Flip,				"Flip"}					
};

XMLPutEnum<CDXPolymerFlipType> sXMLPolymerFlipType(sXMLPolymerFlipTypeValues, sizeof sXMLPolymerFlipTypeValues, kCDXPolymerFlipType_Unspecified);
void XMLPut(XMLDataSink &sink, CDXPolymerFlipType  v) {	sink.os << sXMLPolymerFlipType.lookup(v); }

// **********************
// ** class CDXGraphic **
// **********************
//
// Specialization of CDXObject for CDXGraphic objects

CDXGraphic::CDXGraphic(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_Graphic, id_arg)
	, m_graphicType(kCDXGraphicType_Undefined)
	, m_arrowType(kCDXArrowType_NoHead)
	, m_bracketType(kCDXBracketType_RoundPair)
	, m_rectangleType(kCDXRectangleType_Plain)
	, m_ovalType(kCDXOvalType_Plain)
	, m_symbolType(kCDXSymbolType_LonePair)
	, m_orbitalType(kCDXOrbitalType_s)
	, m_frameType(kCDXFrameType_Unspecified)
	, m_fadePercent(1000)
	, m_headSize(0)
	, m_angularSize(0)
	, m_bracketLipSize(0)
	, m_bracketUsage(kCDXBracketUsage_Unspecified)
	, m_polymerRepeatPattern(kCDXPolymerRepeatPattern_HeadToTail)
	, m_polymerFlipType(kCDXPolymerFlipType_Unspecified)
	, m_has3dHead(false)
	, m_has3dTail(false)
	, m_shadowSize(0)
	, m_hasCenter(false)
	, m_hasMajorAxisEnd(false)
	, m_hasMinorAxisEnd(false)
	, m_cornerRadius(0)
{
}

CDXGraphic::CDXGraphic (const CDXGraphic &src)
	:	CDXGraphicObject		(src)
	,	m_graphicType			(src.m_graphicType)
	,	m_arrowType				(src.m_arrowType)
	,	m_bracketType			(src.m_bracketType)
	,	m_rectangleType			(src.m_rectangleType)
	,	m_ovalType				(src.m_ovalType)
	,	m_symbolType			(src.m_symbolType)
	,	m_orbitalType			(src.m_orbitalType)
	,	m_frameType				(src.m_frameType)
	,	m_headSize				(src.m_headSize)
	,	m_angularSize			(src.m_angularSize)
	,	m_bracketLipSize		(src.m_bracketLipSize)
	,	m_bracketUsage			(src.m_bracketUsage)
	,	m_polymerRepeatPattern	(src.m_polymerRepeatPattern)
	,	m_polymerFlipType		(src.m_polymerFlipType)
	,	m_fadePercent			(src.m_fadePercent)
	,	m_3dHead				(src.m_3dHead)
	,	m_has3dHead				(src.m_has3dHead)
	,	m_3dTail				(src.m_3dTail)
	,	m_has3dTail				(src.m_has3dTail)
	,	m_shadowSize			(src.m_shadowSize)
	,	m_center				(src.m_center)
	,	m_hasCenter				(src.m_hasCenter)
	,	m_majorAxisEnd			(src.m_majorAxisEnd)
	,	m_hasMajorAxisEnd		(src.m_hasMajorAxisEnd)
	,	m_minorAxisEnd			(src.m_minorAxisEnd)
	,	m_hasMinorAxisEnd		(src.m_hasMinorAxisEnd)
	,	m_cornerRadius			(src.m_cornerRadius)
{
}

CDXGraphic::~CDXGraphic()
{
}

CDXObject*	CDXGraphic::Clone() const
{
	return new CDXGraphic (*this);
}

std::string CDXGraphic::XMLObjectName() const
{
	return kCDXML_graphic;
}

void CDXGraphic::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Graphic_Type:
		m_graphicType = (CDXGraphicType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrow_Type:
		m_arrowType = (CDXArrowType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Rectangle_Type:
		m_rectangleType = (CDXRectangleType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Oval_Type:
		m_ovalType = (CDXOvalType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Orbital_Type:
		m_orbitalType = (CDXOrbitalType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Frame_Type:
		m_frameType = (CDXFrameType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_Type:
		m_bracketType = (CDXBracketType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Symbol_Type:
		m_symbolType = (CDXSymbolType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arrowhead_Size:
		m_headSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Arc_AngularSize:
		m_angularSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_CornerRadius:
		m_cornerRadius = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_LipSize:
		m_bracketLipSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bracket_Usage:
		m_bracketUsage = (CDXBracketUsage) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Polymer_RepeatPattern:
		m_polymerRepeatPattern = (CDXPolymerRepeatPattern) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Polymer_FlipType:
		m_polymerFlipType = (CDXPolymerFlipType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_3DHead:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_3dHead.x = src_arg.GetINT32();
		m_3dHead.y = src_arg.GetINT32();
		m_3dHead.z = src_arg.GetINT32();
		m_has3dHead = true;
		break;
		}
	case kCDXProp_3DTail:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_3dTail.x = src_arg.GetINT32();
		m_3dTail.y = src_arg.GetINT32();
		m_3dTail.z = src_arg.GetINT32();
		m_has3dTail = true;
		break;
		}
	case kCDXProp_ShadowSize:
		m_shadowSize = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_3DCenter:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_center.x = src_arg.GetINT32();
		m_center.y = src_arg.GetINT32();
		m_center.z = src_arg.GetINT32();
		m_hasCenter = true;
		break;
		}
	case kCDXProp_3DMajorAxisEnd:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_majorAxisEnd.x = src_arg.GetINT32();
		m_majorAxisEnd.y = src_arg.GetINT32();
		m_majorAxisEnd.z = src_arg.GetINT32();
		m_hasMajorAxisEnd = true;
		break;
		}
	case kCDXProp_3DMinorAxisEnd:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_minorAxisEnd.x = src_arg.GetINT32();
		m_minorAxisEnd.y = src_arg.GetINT32();
		m_minorAxisEnd.z = src_arg.GetINT32();
		m_hasMinorAxisEnd = true;
		break;
		}
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

CDXRectangleType &operator |= (CDXRectangleType &lhs, const CDXRectangleType &rhs);
CDXOvalType &operator |= (CDXOvalType &lhs, const CDXOvalType &rhs);

CDXRectangleType &operator |= (CDXRectangleType &lhs, const CDXRectangleType &rhs)
{
	return lhs = CDXRectangleType(UINT16(lhs) | UINT16(rhs));
}

CDXOvalType &operator |= (CDXOvalType &lhs, const CDXOvalType &rhs)
{
	return lhs = CDXOvalType(UINT16(lhs) | UINT16(rhs));
}

void CDXGraphic::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Graphic_Type:
		m_graphicType = sXMLGraphicType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrow_Type:
		m_arrowType = sXMLArrowType.lookup(attribValue_arg);
		break;
    case kCDXProp_Rectangle_Type:
		m_rectangleType = sXMLRectangleType.lookupTokens(attribValue_arg);
		break;
	case kCDXProp_Oval_Type:
		m_ovalType = sXMLOvalType.lookupTokens(attribValue_arg);
		break;
	case kCDXProp_Orbital_Type:
		m_orbitalType = sXMLOrbitalType.lookup(attribValue_arg);
		break;
	case kCDXProp_Frame_Type:
		m_frameType = sXMLFrameType.lookup(attribValue_arg);
		break;
	case kCDXProp_Bracket_Type:
		m_bracketType = sXMLBracketType.lookup(attribValue_arg);
		break;
	case kCDXProp_Symbol_Type:
		m_symbolType = sXMLSymbolType.lookup(attribValue_arg);
		break;
	case kCDXProp_Arrowhead_Size:
		m_headSize = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Arc_AngularSize:
		m_angularSize = (int)(atof(attribValue_arg.c_str()) * 10.0);
		break;
	case kCDXProp_CornerRadius:
		m_cornerRadius = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_FadePercent:
		m_fadePercent = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Bracket_LipSize:
		m_bracketLipSize = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Bracket_Usage:
		m_bracketUsage = sXMLBracketUsage.lookup(attribValue_arg);
		break;
	case kCDXProp_Polymer_RepeatPattern:
		m_polymerRepeatPattern = sXMLPolymerRepeatPattern.lookup(attribValue_arg);
		break;
	case kCDXProp_Polymer_FlipType:
		m_polymerFlipType = sXMLPolymerFlipType.lookup(attribValue_arg);
		break;
	case kCDXProp_3DHead:
		m_3dHead = StringToCDXPoint3D(attribValue_arg);
		m_has3dHead = true;
		break;
	case kCDXProp_3DTail:
		m_3dTail = StringToCDXPoint3D(attribValue_arg);
		m_has3dTail = true;
		break;
	case kCDXProp_ShadowSize:
		m_shadowSize = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_3DCenter:
		m_center = StringToCDXPoint3D(attribValue_arg);
		m_hasCenter = true;
		break;
	case kCDXProp_3DMajorAxisEnd:
		m_majorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		m_hasMajorAxisEnd = true;
		break;
	case kCDXProp_3DMinorAxisEnd:
		m_minorAxisEnd = StringToCDXPoint3D(attribValue_arg);
		m_hasMinorAxisEnd = true;
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXGraphic::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_Graphic_Type, (INT16) m_graphicType );

	if (m_arrowType != kCDXArrowType_NoHead)
		sink_arg.PutAttribute( kCDXProp_Arrow_Type, (INT16) m_arrowType );

	if (m_graphicType == kCDXGraphicType_Rectangle)
		sink_arg.PutAttribute( kCDXProp_Rectangle_Type, (INT16) m_rectangleType );

	if (m_graphicType == kCDXGraphicType_Oval ||
		(m_graphicType == kCDXGraphicType_Orbital && m_ovalType != kCDXOvalType_Plain))
	{
		sink_arg.PutAttribute( kCDXProp_Oval_Type, (INT16) m_ovalType );
	}

	if (m_graphicType == kCDXGraphicType_Orbital)
		sink_arg.PutAttribute( kCDXProp_Orbital_Type, (INT16) m_orbitalType );

	if (m_graphicType == kCDXGraphicType_Bracket)
		sink_arg.PutAttribute( kCDXProp_Bracket_Type, (INT16) m_bracketType );

	if (m_graphicType == kCDXGraphicType_Symbol)
		sink_arg.PutAttribute( kCDXProp_Symbol_Type, (INT16) m_symbolType );

	if (m_frameType != kCDXFrameType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Frame_Type, (INT16) m_frameType );

	if (m_headSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arrowhead_Size, m_headSize );

	if (m_angularSize != 0)
		sink_arg.PutAttribute( kCDXProp_Arc_AngularSize, m_angularSize );

	if (m_cornerRadius != 0)
		sink_arg.PutAttribute( kCDXProp_CornerRadius, m_cornerRadius );

	if (m_fadePercent != 1000)
		sink_arg.PutAttribute( kCDXProp_FadePercent, m_fadePercent );

	if (m_bracketLipSize != 0)
		sink_arg.PutAttribute( kCDXProp_Bracket_LipSize, m_bracketLipSize );

	if (m_bracketUsage != kCDXBracketUsage_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Bracket_Usage, (INT8) m_bracketUsage );

	if (m_polymerRepeatPattern != kCDXPolymerRepeatPattern_HeadToTail)
		sink_arg.PutAttribute( kCDXProp_Polymer_RepeatPattern, (INT8) m_polymerRepeatPattern );

	if (m_polymerFlipType != kCDXPolymerFlipType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Polymer_FlipType, (INT8) m_polymerFlipType );

	if (m_has3dHead)
		sink_arg.PutAttribute( kCDXProp_3DHead, m_3dHead );

	if (m_has3dTail)
		sink_arg.PutAttribute( kCDXProp_3DTail, m_3dTail );

	if (m_shadowSize != 0)
		sink_arg.PutAttribute( kCDXProp_ShadowSize, m_shadowSize );

	if (m_hasCenter)
		sink_arg.PutAttribute( kCDXProp_3DCenter, m_center );

	if (m_hasMajorAxisEnd)
		sink_arg.PutAttribute( kCDXProp_3DMajorAxisEnd, m_majorAxisEnd );

	if (m_hasMinorAxisEnd)
		sink_arg.PutAttribute( kCDXProp_3DMinorAxisEnd, m_minorAxisEnd );
}

void CDXGraphic::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	CDXMLPutAttribute(sink_arg, kCDXProp_Graphic_Type, m_graphicType );

	if (m_arrowType != kCDXArrowType_NoHead)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrow_Type, m_arrowType );

	if (m_graphicType == kCDXGraphicType_Rectangle)
		CDXMLPutAttribute(sink_arg, kCDXProp_Rectangle_Type, m_rectangleType );

	if (m_graphicType == kCDXGraphicType_Oval ||
		(m_graphicType == kCDXGraphicType_Orbital && m_ovalType != kCDXOvalType_Plain))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Oval_Type, m_ovalType );
	}

	if (m_graphicType == kCDXGraphicType_Orbital)
		CDXMLPutAttribute(sink_arg, kCDXProp_Orbital_Type, m_orbitalType );

	if (m_graphicType == kCDXGraphicType_Bracket)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_Type, m_bracketType );

	if (m_graphicType == kCDXGraphicType_Symbol)
		CDXMLPutAttribute(sink_arg, kCDXProp_Symbol_Type, m_symbolType );

	if (m_frameType != kCDXFrameType_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Frame_Type, m_frameType );

	if (m_headSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arrowhead_Size, m_headSize );

	if (m_angularSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Arc_AngularSize, m_angularSize / 10.0 );

	if (m_cornerRadius != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_CornerRadius, m_cornerRadius );

	if (m_fadePercent != 1000)
		CDXMLPutAttribute(sink_arg, kCDXProp_FadePercent, m_fadePercent );

	if (m_bracketLipSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_LipSize, m_bracketLipSize );

	if (m_bracketUsage != kCDXBracketUsage_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bracket_Usage, m_bracketUsage );

	if (m_polymerRepeatPattern != kCDXPolymerRepeatPattern_HeadToTail)
		CDXMLPutAttribute(sink_arg, kCDXProp_Polymer_RepeatPattern, m_polymerRepeatPattern );

	if (m_polymerFlipType != kCDXPolymerFlipType_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Polymer_FlipType, m_polymerFlipType );

	if (m_has3dHead)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DHead, m_3dHead );

	if (m_has3dTail)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DTail, m_3dTail );

	if (m_shadowSize != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_ShadowSize, m_shadowSize );

	if (m_hasCenter)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DCenter, m_center );

	if (m_hasMajorAxisEnd)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMajorAxisEnd, m_majorAxisEnd );

	if (m_hasMinorAxisEnd)
		CDXMLPutAttribute(sink_arg, kCDXProp_3DMinorAxisEnd, m_minorAxisEnd );
}

std::vector<CDXObjectID> CDXGraphic::GetOtherBracketsInGroup() const
{
	std::vector<CDXObjectID> matchingBrackets;
	const CDXBracketedGroup* group = GetBracketedGroupForBracket();
	if (group)
	{
		CDXObjectsRange attachments = group->ContainedObjects(kCDXObj_BracketAttachment);
		for (CDXObjectsByTag::const_iterator itAttach = attachments.begin(); itAttach != attachments.end(); ++itAttach)
		{
			CDXBracketAttachment* attachment = dynamic_cast<CDXBracketAttachment*>(GetObject(itAttach));
			if (attachment && (attachment->m_graphicID != GetObjectID()))
			{
				// Add all except for the specified bracket to the list
				matchingBrackets.push_back(attachment->m_graphicID);
			}
		}
	}

	return matchingBrackets;
}

const CDXBracketedGroup* CDXGraphic::GetBracketedGroupForBracket() const
{
	const CDXPage* cdxPage = GetCDXPageForObject(*this);
	if (cdxPage)
	{
		CDXObjectsRange groups = cdxPage->ContainedObjects(kCDXObj_BracketedGroup);
		for (CDXObjectsByTag::const_iterator it = groups.begin(); it != groups.end(); ++it)
		{
			CDXBracketedGroup* group = dynamic_cast<CDXBracketedGroup*>(GetObject(it));
			if (group)
			{
				CDXObjectsRange attachments = group->ContainedObjects(kCDXObj_BracketAttachment);
				for (CDXObjectsByTag::const_iterator itAttach = attachments.begin(); itAttach != attachments.end(); ++itAttach)
				{
					CDXBracketAttachment* attachment = dynamic_cast<CDXBracketAttachment*>(GetObject(itAttach));
					if (attachment && (attachment->m_graphicID == GetObjectID()))
					{
						return group;
					}
				}
			}
		}
	}

	return NULL;
}

