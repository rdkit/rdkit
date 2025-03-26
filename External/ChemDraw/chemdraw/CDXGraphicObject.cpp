// CommonCS/LibCommon/Src/CDXGraphicObject.cpp
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
#include "CDXUnicode.h"
#include "CDXMLNames.h"
#include "XMLPutEnum.h"
#include <stdlib.h>	// for atoi and friends
#include <ostream>

extern XMLPutEnum<CDXTextJustification> sXMLTextJustificationText;

XMLPutEnum<CDXAminoAcidTermini>::Value sXMLAminoAcidTerminiValues[] = {
	{kCDXAminoAcidTerminiHOH,		"HOH"},
	{kCDXAminoAcidTerminiNH2COOH,	"NH2COOH"}
};

XMLPutEnum<CDXAminoAcidTermini> sXMLAminoAcidTerminiText(sXMLAminoAcidTerminiValues, sizeof sXMLAminoAcidTerminiValues, kCDXAminoAcidTerminiHOH);
void XMLPut(XMLDataSink &sink, CDXAminoAcidTermini  v) {	sink.os << sXMLAminoAcidTerminiText.asTokens(v); ; }

XMLPutEnum<CDXLineType>::Value sXMLLineTypeValues[] = {
	{kCDXLineType_Solid,	"Solid"},
	{kCDXLineType_Dashed,	"Dashed"},
	{kCDXLineType_Bold,		"Bold"},
	{kCDXLineType_Wavy,		"Wavy"}
};

XMLPutEnum<CDXLineType> sXMLLineType(sXMLLineTypeValues, sizeof sXMLLineTypeValues, kCDXLineType_Solid);
void XMLPut(XMLDataSink &sink, CDXLineType  v) {	sink.os << sXMLLineType.asTokens(v); ; }

XMLPutEnum<CDXFillType>::Value sXMLFillTypeValues[] = {
	{kCDXFillType_Unspecified,	"Unspecified"},
	{kCDXFillType_None,			"None"},
	{kCDXFillType_Solid,		"Solid"},
	{kCDXFillType_Shaded,		"Shaded"}
};

XMLPutEnum<CDXFillType> sXMLFillType(sXMLFillTypeValues, sizeof sXMLFillTypeValues, kCDXFillType_Unspecified);
void XMLPut(XMLDataSink &sink, CDXFillType  v) {	sink.os << sXMLFillType.lookup(v); }

// ****************************
// ** class CDXGraphicObject **
// ****************************
//
// Specialization of CDXObject for CDXGraphicObject objects

CDXGraphicObject::CDXGraphicObject(CDXTag tag_arg, CDXObjectID id_arg)
    : CDXObject(tag_arg, id_arg)
    , m_supersededBy(0)
    , m_zOrder(0)
    , m_fgColor(0)
    , m_bgColor(0)
    , highlightColor(0)
    , m_fgAlpha(0)
    , m_bgAlpha(0)
    , m_visible(true)
    , m_propertiesRepresented(nullptr)
    , m_chemicalProperties(nullptr)
    , m_ignoreErrors(false)
    , m_errorDescription(nullptr)
    , m_fractionalWidths(false)
    , m_interpretChemically(false)
    , m_showAtomQuery(false)
    , m_showAtomStereo(false)
    , m_showAtomEnhancedStereo(false)
    , m_showAtomNumber(false)
    , m_showResidueID(false)
    , m_showBondQuery(false)
    , m_showBondRxn(false)
    , m_showBondStereo(false)
    , m_showTerminalCarbons(false)
    , m_showNonTerminalCarbons(false)
    , m_hideImplicitHydrogens(false)
    , m_hashSpacing(0)
    , m_marginWidth(0)
    , m_lineWidth(0)
    , m_boldWidth(0)
    , m_bondLength(0)
    , m_bondSpacing(0)
    , m_bondSpacingAbs(0)
    , m_bondSpacingType(kCDXBondSpacingRelative)
    , m_chainAngle(0)
    , m_labelLineHeight(0)
    , m_captionLineHeight(0)
    , m_labelJustification(kCDXTextJustification_Auto)
    , m_captionJustification(kCDXTextJustification_Auto)
    , m_aminoAcidTermini(kCDXAminoAcidTerminiHOH)
    , m_showSequenceTermini(true)
    , m_showSequenceBonds(true)
    , m_showSequenceUnlinkedBranches(false)
    , m_residueWrapCount(40)
    , m_residueBlockCount(40)
    , m_lineType(kCDXLineType_Solid)
    , m_fillType(kCDXFillType_Unspecified)
    , m_flags(0)
    , m_flags2(0)
    , m_rxnAutonumberStyle(kCDXRxnAutonumberStyle_Roman)
    , m_rxnAutonumberStart(1)
    , m_rxnAutonumberConditions(false)
    , m_rxnAutonumberFormat("(#)")
    , showPerspective(false)
{
}

CDXGraphicObject::~CDXGraphicObject()
{
	DeleteAndNull(m_propertiesRepresented);
	DeleteAndNull(m_chemicalProperties);
	DeleteAndNull(m_errorDescription);
}

CDXGraphicObject::CDXGraphicObject (const CDXGraphicObject &src)
	:	CDXObject			(src)
	,	m_supersededBy		(src.m_supersededBy)
	,	m_2dPosition		(src.m_2dPosition)
	,	m_3dPosition		(src.m_3dPosition)
	,	m_boundingBox		(src.m_boundingBox)
	,	m_zOrder			(src.m_zOrder)
	,	m_fgColor			(src.m_fgColor)
	,	m_bgColor			(src.m_bgColor)
    ,   highlightColor(src.highlightColor)
	,	m_fgAlpha			(src.m_fgAlpha)
	,	m_bgAlpha			(src.m_bgAlpha)
	,	m_visible			(src.m_visible)
	,	m_propertiesRepresented (src.m_propertiesRepresented == NULL ? NULL : new std::vector<CDXPropRep>(*src.m_propertiesRepresented))
	,	m_chemicalProperties	(src.m_chemicalProperties == NULL ? NULL : new std::vector<CDXChemProp>(*src.m_chemicalProperties))
	,	m_ignoreErrors		(src.m_ignoreErrors)
	,	m_errorDescription	(src.m_errorDescription == NULL ? NULL : new CDXString(*src.m_errorDescription))
	,	m_fractionalWidths	(src.m_fractionalWidths)
	,	m_interpretChemically	(src.m_interpretChemically)
	,	m_showAtomQuery		(src.m_showAtomQuery)
	,	m_showAtomStereo	(src.m_showAtomStereo)
	,	m_showAtomEnhancedStereo	(src.m_showAtomEnhancedStereo)
	,	m_showAtomNumber	(src.m_showAtomNumber)
	,	m_showResidueID		(src.m_showResidueID)
	,	m_showBondQuery		(src.m_showBondQuery)
	,	m_showBondRxn		(src.m_showBondRxn)
	,	m_showBondStereo	(src.m_showBondStereo)
	,	m_showTerminalCarbons	(src.m_showTerminalCarbons)
	,	m_showNonTerminalCarbons(src.m_showNonTerminalCarbons)
	,	m_hideImplicitHydrogens	(src.m_hideImplicitHydrogens)
	,	m_labelStyle		(src.m_labelStyle)
	,	m_captionStyle		(src.m_captionStyle)
	,	m_hashSpacing		(src.m_hashSpacing)
	,	m_marginWidth		(src.m_marginWidth)
	,	m_lineWidth			(src.m_lineWidth)
	,	m_boldWidth			(src.m_boldWidth)
	,	m_bondLength		(src.m_bondLength)
	,	m_bondSpacing		(src.m_bondSpacing)
	,	m_bondSpacingAbs	(src.m_bondSpacingAbs)
	,	m_bondSpacingType	(src.m_bondSpacingType)
	,	m_chainAngle		(src.m_chainAngle)
	,	m_labelLineHeight	(src.m_labelLineHeight)
	,	m_captionLineHeight	(src.m_captionLineHeight)
	,	m_labelJustification(src.m_labelJustification)
	,	m_captionJustification	(src.m_captionJustification)
	,   m_aminoAcidTermini		(src.m_aminoAcidTermini)
	,	m_showSequenceTermini	(src.m_showSequenceTermini)
	,   m_showSequenceUnlinkedBranches    (src.m_showSequenceUnlinkedBranches)
	,	m_showSequenceBonds	(src.m_showSequenceBonds)
	,	m_lineType			(src.m_lineType)
	,	m_fillType			(src.m_fillType)
	,	m_flags				(src.m_flags)
	,	m_flags2			(src.m_flags2)
	,	m_residueWrapCount	(src.m_residueWrapCount)
	,	m_residueBlockCount	(src.m_residueBlockCount)
    ,   m_rxnAutonumberStyle(src.m_rxnAutonumberStyle)
    ,   m_rxnAutonumberConditions(src.m_rxnAutonumberConditions)
    ,   m_rxnAutonumberStart(src.m_rxnAutonumberStart)
    ,   m_rxnAutonumberFormat(src.m_rxnAutonumberFormat)
    ,   showPerspective(src.showPerspective)
{
}

void CDXGraphicObject::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_SupersededBy:
		m_supersededBy = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_ZOrder:
		ZOrder (src_arg.GetUINT (size_arg));
		break;
	case kCDXProp_2DPosition:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		Position(CDXPoint2D(x, y));
		break;
		}
	case kCDXProp_3DPosition:
		{
		if (size_arg != 12)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate x = src_arg.GetINT32();
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate z = src_arg.GetINT32();
		Position3D(CDXPoint3D(x, y, z));
		break;
		}
	case kCDXProp_BoundingBox:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetINT32();
		CDXCoordinate left   = src_arg.GetINT32();
		CDXCoordinate bottom = src_arg.GetINT32();
		CDXCoordinate right  = src_arg.GetINT32();
		BoundingBox(CDXRectangle(top, left, bottom, right));
		break;
		}
	case kCDXProp_ForegroundColor:
		{
		FGColor(UINT16(src_arg.GetUINT(size_arg)));
		break;
		}
	case kCDXProp_BackgroundColor:
		{
		BGColor(UINT16(src_arg.GetUINT(size_arg)));
		break;
		}
	case kCDXProp_ForegroundAlpha:
		{
		m_fgAlpha = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_fgAlpha;
		break;
		}
	case kCDXProp_BackgroundAlpha:
		{
		m_bgAlpha = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_bgAlpha;
		break;
		}
	case kCDXProp_HighlightColor:
		SetHighlightColor(CDXColorIndex(src_arg.GetUINT(size_arg)));
		break;
            
    case kCDXProp_StructurePerspective:
        SetShowPerspective((size_arg != 0) && (src_arg.GetINT(size_arg) != 0));
        break;
            
	case kCDXProp_Visible:
		m_visible = (size_arg != 0) && src_arg.GetINT(size_arg) != 0;
		break;
	case kCDXProp_RepresentsProperty:
		{
		if ((size_arg % 6) != 0)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		if (size_arg > 0)
		{
			for (ptrdiff_t i = size_arg;  i > 0;  i -= 6)
			{
				CDXObjectID objectID = src_arg.GetINT32();
				CDXTag propTag = src_arg.GetINT16();
				PropertyRepresented(CDXPropRep(objectID, propTag));
			}
		}
		break;
		}
	case kCDXProp_IgnoreWarnings:
		m_ignoreErrors = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		break;
	case kCDXProp_ChemicalWarning:
		if (size_arg > 0)
		{
			DeleteAndNull(m_errorDescription);
			m_errorDescription = new CDXString;
			m_errorDescription->Read(src_arg, size_arg);
		}
		break;
	case kCDXProp_ChainAngle:
		m_chainAngle = src_arg.GetUINT(size_arg);
		if (size_arg == 2) // some old files stored this in 10ths of a degree
			m_chainAngle = (m_chainAngle * 65536) / 10;
		m_flags |= has_chainAngle;
		break;
	case kCDXProp_BondSpacing:
		m_bondSpacing = src_arg.GetUINT(size_arg);
		m_flags |= has_bondSpacing;
		break;
	case kCDXProp_BondSpacingAbs:
		m_bondSpacingAbs = src_arg.GetUINT(size_arg);
		m_flags |= has_bondSpacingAbs;
		break;
	case kCDXProp_BondLength:
		m_bondLength = src_arg.GetUINT(size_arg);
		m_flags |= has_bondLength;
		break;
	case kCDXProp_BoldWidth:
		m_boldWidth = src_arg.GetUINT(size_arg);
		m_flags |= has_boldWidth;
		break;
	case kCDXProp_LineWidth:
		m_lineWidth = src_arg.GetUINT(size_arg);
		m_flags |= has_lineWidth;
		break;
	case kCDXProp_MarginWidth:
		m_marginWidth = src_arg.GetUINT(size_arg);
		m_flags |= has_marginWidth;
		break;
	case kCDXProp_HashSpacing:
		m_hashSpacing = src_arg.GetUINT(size_arg);
		m_flags |= has_hashSpacing;
		break;
	case kCDXProp_LabelStyle:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_labelStyle.family = src_arg.GetUINT16();
		m_labelStyle.face   = src_arg.GetUINT16();
		m_labelStyle.size   = src_arg.GetUINT16();
		m_labelStyle.color  = src_arg.GetUINT16();
		m_flags |= has_labelStyle;
		break;
		}
	case kCDXProp_CaptionStyle:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_captionStyle.family = src_arg.GetUINT16();
		m_captionStyle.face   = src_arg.GetUINT16();
		m_captionStyle.size   = src_arg.GetUINT16();
		m_captionStyle.color  = src_arg.GetUINT16();
		m_flags |= has_captionStyle;
		break;
		}
	case kCDXProp_FractionalWidths:
		m_fractionalWidths = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_fractionalWidths;
		break;
	case kCDXProp_InterpretChemically:
		m_interpretChemically = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_interpretChemically;
		break;
	case kCDXProp_Atom_ShowQuery:
		m_showAtomQuery = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showAtomQuery;
		break;
	case kCDXProp_Atom_ShowStereo:
		m_showAtomStereo = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showAtomStereo;
		break;
	case kCDXProp_Atom_ShowEnhancedStereo:
		m_showAtomEnhancedStereo = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showAtomEnhancedStereo;
		break;
	case kCDXProp_Atom_ShowAtomNumber:
		m_showAtomNumber = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showAtomNumber;
		break;
	case kCDXProp_Atom_ShowResidueID:
		m_showResidueID = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_showResidueID;
		break;
	case kCDXProp_Bond_ShowQuery:
		m_showBondQuery = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showBondQuery;
		break;
	case kCDXProp_Bond_ShowRxn:
		m_showBondRxn = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showBondRxn;
		break;
	case kCDXProp_Bond_ShowStereo:
		m_showBondStereo = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showBondStereo;
		break;
	case kCDXProp_Atom_ShowTerminalCarbonLabels:
		m_showTerminalCarbons = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showTerminalCarbons;
		break;
	case kCDXProp_Atom_ShowNonTerminalCarbonLabels:
		m_showNonTerminalCarbons = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_showNonTerminalCarbons;
		break;
	case kCDXProp_Atom_HideImplicitHydrogens:
		m_hideImplicitHydrogens = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_hideImplicitHydrogens;
		break;
	case kCDXProp_LabelJustification:
		m_labelJustification = (CDXTextJustification) src_arg.GetINT(size_arg);
        m_flags |= has_labelJustification;
		break;
	case kCDXProp_CaptionJustification:
		m_captionJustification = (CDXTextJustification) src_arg.GetINT(size_arg);
        m_flags |= has_captionJustification;
		break;
	case kCDXProp_AminoAcidTermini:
		m_aminoAcidTermini = (CDXAminoAcidTermini) src_arg.GetINT(size_arg);
		m_flags2 |= has_aminoAcidTermini;
		break;
	case kCDXProp_ShowSequenceTermini:
		m_showSequenceTermini = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_showSequenceTermini;
		break;
	case kCDXProp_ShowSequenceBonds:
		m_showSequenceBonds = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_showSequenceBonds;
		break;
	case kCDXProp_ShowSequenceUnlinkedBranches:
		m_showSequenceUnlinkedBranches = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_showSequenceUnlinkedBranches;
		break;
	case kCDXProp_ResidueWrapCount:
		m_residueWrapCount = src_arg.GetINT(size_arg);
		m_flags2 |= has_residueWrapCount;
		break;
	case kCDXProp_ResidueBlockCount:
		m_residueBlockCount = src_arg.GetINT(size_arg);
		m_flags2 |= has_residueBlockCount;
		break;
	case kCDXProp_LabelLineHeight:
		m_labelLineHeight = src_arg.GetUINT(size_arg);
        if ((size_arg == 2) &&
            (m_labelLineHeight != kCDXLineHeight_Automatic) &&
            (m_labelLineHeight != kCDXLineHeight_Variable))
        {
            // pre-6.0 files have this in 1440'ths of an inch
            m_labelLineHeight = (m_labelLineHeight * 0x10000) / 20; 
        }

        m_flags |= has_labelLineHeight;
		break;
	case kCDXProp_CaptionLineHeight:
		m_captionLineHeight = src_arg.GetUINT(size_arg);
        if ((size_arg == 2) && 
            (m_captionLineHeight != kCDXLineHeight_Automatic) && 
            (m_captionLineHeight != kCDXLineHeight_Variable))
        {
            // pre-6.0 files have this in 1440'ths of an inch
            m_captionLineHeight = (m_captionLineHeight * 0x10000) / 20; 
        }

        m_flags |= has_captionLineHeight;
		break;
	case kCDXProp_Line_Type:
		m_lineType = (CDXLineType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Fill_Type:
		m_fillType = (CDXFillType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_ChemicalPropertyName:
    case kCDXProp_ChemicalPropertyID:
    case kCDXProp_ChemicalPropertyFragmentLabel:
	case kCDXProp_ChemicalPropertyFormula:
	case kCDXProp_ChemicalPropertyExactMass:
	case kCDXProp_ChemicalPropertyMolWeight:
	case kCDXProp_ChemicalPropertyMOverZ:
	case kCDXProp_ChemicalPropertyAnalysis:
	case kCDXProp_ChemicalPropertyBoilingPoint:
	case kCDXProp_ChemicalPropertyMeltingPoint:
	case kCDXProp_ChemicalPropertyCriticalTemp:
	case kCDXProp_ChemicalPropertyCriticalPressure:
	case kCDXProp_ChemicalPropertyCriticalVolume:
	case kCDXProp_ChemicalPropertyGibbsEnergy:
	case kCDXProp_ChemicalPropertyLogP:
	case kCDXProp_ChemicalPropertyMR:
	case kCDXProp_ChemicalPropertyHenrysLaw:
	case kCDXProp_ChemicalPropertyHeatOfForm:
	case kCDXProp_ChemicalPropertytPSA:
	case kCDXProp_ChemicalPropertyCLogP:
	case kCDXProp_ChemicalPropertyCMR:
	case kCDXProp_ChemicalPropertyLogS:
	case kCDXProp_ChemicalPropertyPKa:
		{
			CDXString str;
			str.Read(src_arg, size_arg);
			ChemicalProperty(CDXChemProp(attribTag_arg, str));
		}
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

// Figure out the bounds of this object and all its children
CDXRectangle *CDXGraphicObject::BoundsIncludingChildren() const
{
	CDXRectangle *r = NULL;
	if (KnownBoundingBox())
	{
		r = new CDXRectangle(BoundingBox());
	}
	else if (KnownPosition())
	{
		r = new CDXRectangle;
		r->left = r->right = Position().x;
		r->top = r->bottom = Position().y;
	}

	std::unique_ptr<CDXRectangle> subBounds(CDXObject::BoundsIncludingChildren());
	if (subBounds.get() != NULL)
	{
		if (r == NULL)
		{
			r = subBounds.release();
		}
		else
		{
			r->left = ((r->left < subBounds->left) ? r->left : subBounds->left);
			r->top = ((r->top < subBounds->top) ? r->top : subBounds->top);
			r->right = ((r->right > subBounds->right) ? r->right : subBounds->right);
			r->bottom = ((r->bottom > subBounds->bottom) ? r->bottom : subBounds->bottom);
		}
	}
	return r;
}

CDXLineType &operator |= (CDXLineType &lhs, const CDXLineType &rhs)
{
	return lhs = CDXLineType(UINT16(lhs) | UINT16(rhs));
}

void CDXGraphicObject::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_SupersededBy:
		m_supersededBy = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_ZOrder:
		ZOrder( atoi(attribValue_arg.c_str()) );
		break;
	case kCDXProp_2DPosition:
		Position(StringToCDXPoint2D(attribValue_arg));
		break;
	case kCDXProp_3DPosition:
		Position3D(StringToCDXPoint3D(attribValue_arg));
		break;
	case kCDXProp_BoundingBox:
		BoundingBox(StringToCDXRectangle(attribValue_arg));
		break;
	case kCDXProp_ForegroundColor:
		FGColor( atoi(attribValue_arg.c_str()) );
		break;
	case kCDXProp_BackgroundColor:
		BGColor( atoi(attribValue_arg.c_str()) );
		break;
	case kCDXProp_ForegroundAlpha:
		m_fgAlpha = atoi(attribValue_arg.c_str()) * 65535;
		m_flags2 |= has_fgAlpha;
		break;
	case kCDXProp_BackgroundAlpha:
		m_bgAlpha = atoi(attribValue_arg.c_str()) * 65535;
		m_flags2 |= has_bgAlpha;
		break;
            
    case kCDXProp_HighlightColor:
         SetHighlightColor(atoi(attribValue_arg.c_str()));
         break;
            
    case kCDXProp_StructurePerspective:
        SetShowPerspective((attribValue_arg == "yes") || (attribValue_arg == "true"));
        break;
            
	case kCDXProp_Visible:
		m_visible = (attribValue_arg == "yes") || (attribValue_arg == "true");
		break;
	case kCDXProp_IgnoreWarnings:
		m_ignoreErrors = (attribValue_arg == "yes") || (attribValue_arg == "true");
		break;
	case kCDXProp_ChemicalWarning:
		DeleteAndNull(m_errorDescription);
		m_errorDescription = new CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_ChainAngle:
		m_chainAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
		m_flags |= has_chainAngle;
		break;
	case kCDXProp_BondSpacing:
		m_bondSpacing = (int)(10 * atof(attribValue_arg.c_str()));
		m_flags |= has_bondSpacing;
		break;
	case kCDXProp_BondSpacingAbs:
		m_bondSpacingAbs = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_bondSpacingAbs;
		break;
	case kCDXProp_BondLength:
		m_bondLength = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_bondLength;
		break;
	case kCDXProp_BoldWidth:
		m_boldWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_boldWidth;
		break;
	case kCDXProp_LineWidth:
		m_lineWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_lineWidth;
		break;
	case kCDXProp_MarginWidth:
		m_marginWidth = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_marginWidth;
		break;
	case kCDXProp_HashSpacing:
		m_hashSpacing = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_hashSpacing;
		break;
	case kCDXProp_LabelStyleFont:
		m_labelStyle.family = atoi(attribValue_arg.c_str());
		m_flags |= has_labelStyle;
		break;
	case kCDXProp_LabelStyleFace:
		m_labelStyle.face = atoi(attribValue_arg.c_str());
		m_flags |= has_labelStyle;
		break;
	case kCDXProp_LabelStyleSize:
		m_labelStyle.size = (int)(20.0 * atof(attribValue_arg.c_str()));
		m_flags |= has_labelStyle;
		break;
	case kCDXProp_LabelStyleColor:
		m_labelStyle.color = atoi(attribValue_arg.c_str());
		m_flags |= has_labelStyle;
		break;
	case kCDXProp_CaptionStyleFont:
		m_captionStyle.family = atoi(attribValue_arg.c_str());
		m_flags |= has_captionStyle;
		break;
	case kCDXProp_CaptionStyleFace:
		m_captionStyle.face = atoi(attribValue_arg.c_str());
		m_flags |= has_captionStyle;
		break;
	case kCDXProp_CaptionStyleSize:
		m_captionStyle.size = (int)(20.0 * atof(attribValue_arg.c_str()));
		m_flags |= has_captionStyle;
		break;
	case kCDXProp_CaptionStyleColor:
		m_captionStyle.color = atoi(attribValue_arg.c_str());
		m_flags |= has_captionStyle;
		break;
	case kCDXProp_LabelJustification:
		m_labelJustification = sXMLTextJustificationText.lookup(attribValue_arg);
		m_flags |= has_labelJustification;
		break;
	case kCDXProp_CaptionJustification:
		m_captionJustification = sXMLTextJustificationText.lookup(attribValue_arg);
		m_flags |= has_captionJustification;
		break;
	case kCDXProp_AminoAcidTermini:
		m_aminoAcidTermini = sXMLAminoAcidTerminiText.lookup(attribValue_arg);
		m_flags2 |= has_aminoAcidTermini;
		break;
	case kCDXProp_ShowSequenceTermini:
		m_showSequenceTermini = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_showSequenceTermini;
		break;
	case kCDXProp_ShowSequenceBonds:
		m_showSequenceBonds = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_showSequenceBonds;
		break;
	case kCDXProp_ShowSequenceUnlinkedBranches:
		m_showSequenceUnlinkedBranches = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_showSequenceUnlinkedBranches;
		break;
	case kCDXProp_ResidueWrapCount:
		m_residueWrapCount = atoi(attribValue_arg.c_str());
		m_flags2 |= has_residueWrapCount;
		break;
	case kCDXProp_ResidueBlockCount:
		m_residueBlockCount = atoi(attribValue_arg.c_str());
		m_flags2 |= has_residueBlockCount;
		break;
	case kCDXProp_LabelLineHeight:
		if (attribValue_arg == kCDXML_auto)
			m_labelLineHeight = kCDXLineHeight_Automatic;
		else if (attribValue_arg == kCDXML_variable)
			m_labelLineHeight = kCDXLineHeight_Variable;
		else
			m_labelLineHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_labelLineHeight;
		break;
	case kCDXProp_CaptionLineHeight:
		if (attribValue_arg == kCDXML_auto)
			m_captionLineHeight = kCDXLineHeight_Automatic;
		else if (attribValue_arg == kCDXML_variable)
			m_captionLineHeight = kCDXLineHeight_Variable;
		else
			m_captionLineHeight = CDXCoordinatefromPoints(atof(attribValue_arg.c_str()));
		m_flags |= has_captionLineHeight;
		break;
	case kCDXProp_FractionalWidths:
		m_fractionalWidths = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_fractionalWidths;
		break;
	case kCDXProp_InterpretChemically:
		m_interpretChemically = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_interpretChemically;
		break;
	case kCDXProp_Atom_ShowQuery:
		m_showAtomQuery = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showAtomQuery;
		break;
	case kCDXProp_Atom_ShowStereo:
		m_showAtomStereo = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showAtomStereo;
		break;
	case kCDXProp_Atom_ShowEnhancedStereo:
		m_showAtomEnhancedStereo = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showAtomEnhancedStereo;
		break;
	case kCDXProp_Atom_ShowAtomNumber:
		m_showAtomNumber = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showAtomNumber;
		break;
	case kCDXProp_Atom_ShowResidueID:
		m_showResidueID = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_showResidueID;
		break;
	case kCDXProp_Bond_ShowQuery:
		m_showBondQuery = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showBondQuery;
		break;
	case kCDXProp_Bond_ShowRxn:
		m_showBondRxn = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showBondRxn;
		break;
	case kCDXProp_Bond_ShowStereo:
		m_showBondStereo = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showBondStereo;
		break;
	case kCDXProp_Atom_ShowTerminalCarbonLabels:
		m_showTerminalCarbons = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showTerminalCarbons;
		break;
	case kCDXProp_Atom_ShowNonTerminalCarbonLabels:
		m_showNonTerminalCarbons = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_showNonTerminalCarbons;
		break;
	case kCDXProp_Atom_HideImplicitHydrogens:
		m_hideImplicitHydrogens = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_hideImplicitHydrogens;
		break;
	case kCDXProp_Line_Type:
		m_lineType = sXMLLineType.lookupTokens(attribValue_arg);
		break;
	case kCDXProp_Fill_Type:
		m_fillType = sXMLFillType.lookup(attribValue_arg);
		break;
    case kCDXProp_ChemicalPropertyName:
    case kCDXProp_ChemicalPropertyID:
    case kCDXProp_ChemicalPropertyFragmentLabel:
	case kCDXProp_ChemicalPropertyFormula:
	case kCDXProp_ChemicalPropertyExactMass:
	case kCDXProp_ChemicalPropertyMolWeight:
	case kCDXProp_ChemicalPropertyMOverZ:
	case kCDXProp_ChemicalPropertyAnalysis:
	case kCDXProp_ChemicalPropertyBoilingPoint:
	case kCDXProp_ChemicalPropertyMeltingPoint:
	case kCDXProp_ChemicalPropertyCriticalTemp:
	case kCDXProp_ChemicalPropertyCriticalPressure:
	case kCDXProp_ChemicalPropertyCriticalVolume:
	case kCDXProp_ChemicalPropertyGibbsEnergy:
	case kCDXProp_ChemicalPropertyLogP:
	case kCDXProp_ChemicalPropertyMR:
	case kCDXProp_ChemicalPropertyHenrysLaw:
	case kCDXProp_ChemicalPropertyHeatOfForm:
	case kCDXProp_ChemicalPropertytPSA:
	case kCDXProp_ChemicalPropertyCLogP:
	case kCDXProp_ChemicalPropertyCMR:
	case kCDXProp_ChemicalPropertyLogS:
	case kCDXProp_ChemicalPropertyPKa:
		{
			CDXString str;
			str = CDXString(UnicodeToString(attribValue_arg));
			ChemicalProperty(CDXChemProp(attribTag_arg, str));
			break;
		}
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXGraphicObject::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	if (m_supersededBy != 0)
		sink_arg.PutAttribute( kCDXProp_SupersededBy, m_supersededBy );

	if (KnownPosition())
		sink_arg.PutAttribute( kCDXProp_2DPosition, Position() );

	if (KnownPosition3D())
		sink_arg.PutAttribute( kCDXProp_3DPosition, Position3D() );

	if (KnownBoundingBox())
		sink_arg.PutAttribute( kCDXProp_BoundingBox, BoundingBox() );

	if (KnownZOrder())
		sink_arg.PutAttribute( kCDXProp_ZOrder, (INT16)ZOrder() );

	if (KnownFGColor())
		sink_arg.PutAttribute( kCDXProp_ForegroundColor, FGColor() );

	if (HasHighlightColor())
	{
		sink_arg.PutAttribute(kCDXProp_HighlightColor, GetHighlightColor());
	}

	if (KnownBGColor())
		sink_arg.PutAttribute( kCDXProp_BackgroundColor, BGColor() );

	if (Known(has_fgAlpha))
		sink_arg.PutAttribute(kCDXProp_ForegroundAlpha, m_fgAlpha);

	if (Known(has_bgAlpha))
		sink_arg.PutAttribute(kCDXProp_BackgroundAlpha, m_bgAlpha);

	if (!m_visible)
		sink_arg.PutAttribute(kCDXProp_Visible, (INT8)m_visible);
    
    if (showPerspective)
    {
        sink_arg.PutAttribute(kCDXProp_StructurePerspective, (INT8)showPerspective);
    }

	if (m_propertiesRepresented != NULL && !m_propertiesRepresented->empty())
		sink_arg.PutAttribute( kCDXProp_RepresentsProperty, *m_propertiesRepresented );

	if (m_chemicalProperties != NULL && !m_chemicalProperties->empty())
		sink_arg.PutAttribute( *m_chemicalProperties );

	if (m_ignoreErrors)
		sink_arg.PutAttribute(kCDXProp_IgnoreWarnings);

	if (m_errorDescription != NULL && !m_errorDescription->empty())
		sink_arg.PutAttributeCDXString(kCDXProp_ChemicalWarning, *m_errorDescription);

	if (Known(has_fractionalWidths))
		sink_arg.PutAttribute( kCDXProp_FractionalWidths, (INT8)m_fractionalWidths );
	
	if (Known(has_interpretChemically))
		sink_arg.PutAttribute( kCDXProp_InterpretChemically, (INT8)m_interpretChemically );
	
	if (Known(has_showAtomQuery))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowQuery, (INT8)m_showAtomQuery);
	
	if (Known(has_showAtomStereo))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowStereo, (INT8)m_showAtomStereo );
	
	if (Known(has_showAtomEnhancedStereo))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowEnhancedStereo, (INT8)m_showAtomEnhancedStereo );
	
	if (Known(has_showAtomNumber))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowAtomNumber, (INT8)m_showAtomNumber );
	
	if (Known(has_showResidueID))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowResidueID, (INT8)m_showResidueID );
	
	if (Known(has_showBondQuery))
		sink_arg.PutAttribute( kCDXProp_Bond_ShowQuery, (INT8)m_showBondQuery );
	
	if (Known(has_showBondRxn))
		sink_arg.PutAttribute( kCDXProp_Bond_ShowRxn, (INT8)m_showBondRxn );
	
	if (Known(has_showBondStereo))
		sink_arg.PutAttribute( kCDXProp_Bond_ShowStereo, (INT8)m_showBondStereo );
	
	if (Known(has_showTerminalCarbons))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowTerminalCarbonLabels, (INT8)m_showTerminalCarbons );
	
	if (Known(has_showNonTerminalCarbons))
		sink_arg.PutAttribute( kCDXProp_Atom_ShowNonTerminalCarbonLabels, (INT8)m_showNonTerminalCarbons );
	
	if (Known(has_hideImplicitHydrogens))
		sink_arg.PutAttribute( kCDXProp_Atom_HideImplicitHydrogens, (INT8)m_hideImplicitHydrogens );
	
	if (m_labelStyle.family != UINT16(-1))
		sink_arg.PutAttribute( kCDXProp_LabelStyle, m_labelStyle );
	
	if (m_captionStyle.family != UINT16(-1))
		sink_arg.PutAttribute( kCDXProp_CaptionStyle, m_captionStyle );
	
	if (Known(has_hashSpacing))
		sink_arg.PutAttribute( kCDXProp_HashSpacing, m_hashSpacing );
	
	if (Known(has_marginWidth))
		sink_arg.PutAttribute( kCDXProp_MarginWidth, m_marginWidth );
	
	if (Known(has_lineWidth))
		sink_arg.PutAttribute( kCDXProp_LineWidth, m_lineWidth );
	
	if (Known(has_boldWidth))
		sink_arg.PutAttribute( kCDXProp_BoldWidth, m_boldWidth );
	
	if (Known(has_bondLength))
		sink_arg.PutAttribute( kCDXProp_BondLength, m_bondLength );
	
	if (Known(has_bondSpacing))
		sink_arg.PutAttribute( kCDXProp_BondSpacing, m_bondSpacing );
	
	if (Known(has_bondSpacingAbs))
		sink_arg.PutAttribute( kCDXProp_BondSpacingAbs, m_bondSpacingAbs );
	
	if (Known(has_chainAngle))
		sink_arg.PutAttribute( kCDXProp_ChainAngle, m_chainAngle );

	if (Known(has_labelJustification))
		sink_arg.PutAttribute( kCDXProp_LabelJustification, (UINT8) m_labelJustification );

	if (Known(has_captionJustification))
		sink_arg.PutAttribute( kCDXProp_CaptionJustification, (UINT8) m_captionJustification );

	if (Known(has_aminoAcidTermini))
		sink_arg.PutAttribute( kCDXProp_AminoAcidTermini, (UINT8) m_aminoAcidTermini );

	if (Known(has_showSequenceTermini))
		sink_arg.PutAttribute( kCDXProp_ShowSequenceTermini, (INT8)m_showSequenceTermini );
	
	if (Known(has_showSequenceBonds))
		sink_arg.PutAttribute( kCDXProp_ShowSequenceBonds, (INT8)m_showSequenceBonds );

	if (Known(has_showSequenceUnlinkedBranches))
		sink_arg.PutAttribute( kCDXProp_ShowSequenceUnlinkedBranches, (INT8)m_showSequenceUnlinkedBranches );
	
	if (Known(has_residueWrapCount))
		sink_arg.PutAttribute( kCDXProp_ResidueWrapCount, (INT8)m_residueWrapCount );
	
	if (Known(has_residueBlockCount))
		sink_arg.PutAttribute( kCDXProp_ResidueBlockCount, (INT8)m_residueBlockCount );
	
	if (Known(has_labelLineHeight))
	{
		// Put this out in 1440'ths of an inch, for compatibility with older ChemDraw versions
		if (m_labelLineHeight == kCDXLineHeight_Automatic)
			sink_arg.PutAttribute( kCDXProp_LabelLineHeight, INT16(kCDXLineHeight_Automatic) );
		else if (m_labelLineHeight != kCDXLineHeight_Variable)
			sink_arg.PutAttribute( kCDXProp_LabelLineHeight, INT16(m_labelLineHeight * 20.0 / 65536.0 + 0.5) );
	}

	if (Known(has_captionLineHeight))
	{
		// Put this out in 1440'ths of an inch, for compatibility with older ChemDraw versions
		if (m_captionLineHeight == kCDXLineHeight_Variable)
			sink_arg.PutAttribute( kCDXProp_CaptionLineHeight, INT16(kCDXLineHeight_Variable) );
		else if (m_captionLineHeight != kCDXLineHeight_Automatic)
			sink_arg.PutAttribute( kCDXProp_CaptionLineHeight, INT16(m_captionLineHeight * 20.0 / 65536.0 + 0.5) );
	}

	if (m_lineType != kCDXLineType_Solid)
		sink_arg.PutAttribute( kCDXProp_Line_Type, (INT16) m_lineType );

	if (m_fillType != kCDXFillType_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Fill_Type, (INT16) m_fillType );

}

void CDXGraphicObject::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	if (m_supersededBy != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_SupersededBy, m_supersededBy );

	if (KnownPosition())
		CDXMLPutAttribute(sink_arg, kCDXProp_2DPosition, Position() );

	if (KnownPosition3D())
		CDXMLPutAttribute(sink_arg, kCDXProp_3DPosition, Position3D() );

	if (KnownBoundingBox())
		CDXMLPutAttribute(sink_arg, kCDXProp_BoundingBox, BoundingBox() );

	if (KnownZOrder())
		CDXMLPutAttribute(sink_arg, kCDXProp_ZOrder, ZOrder() );

	if (KnownFGColor())
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundColor, FGColor() );

	if (KnownBGColor())
		CDXMLPutAttribute(sink_arg, kCDXProp_BackgroundColor, BGColor() );

	if (Known(has_fgAlpha))
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundAlpha, CDXValueToString(m_fgAlpha / 65535.0, 4));

	if (Known(has_bgAlpha))
		CDXMLPutAttribute(sink_arg, kCDXProp_BackgroundAlpha, CDXValueToString(m_bgAlpha / 65535.0, 4));

	if (HasHighlightColor())
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_HighlightColor, GetHighlightColor());
	}

    if (showPerspective)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_StructurePerspective, GetShowPerspective());
    }
    
	if (!m_visible)
		CDXMLPutAttribute(sink_arg, kCDXProp_Visible, m_visible);

	if (m_ignoreErrors)
		CDXMLPutAttribute(sink_arg, kCDXProp_IgnoreWarnings);

	if (m_errorDescription != NULL && !m_errorDescription->empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_ChemicalWarning, StringToUnicode("", m_errorDescription->str()));

	if (Known(has_fractionalWidths))
		CDXMLPutAttribute(sink_arg, kCDXProp_FractionalWidths, m_fractionalWidths );
	
	if (Known(has_interpretChemically))
		CDXMLPutAttribute(sink_arg, kCDXProp_InterpretChemically, m_interpretChemically );
	
	if (Known(has_showTerminalCarbons))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowTerminalCarbonLabels, m_showTerminalCarbons );
	
	if (Known(has_showNonTerminalCarbons))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowNonTerminalCarbonLabels, m_showNonTerminalCarbons );
	
	if (Known(has_hideImplicitHydrogens))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_HideImplicitHydrogens, m_hideImplicitHydrogens );
	
	if (Known(has_showAtomQuery))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowQuery, m_showAtomQuery );
	
	if (Known(has_showAtomStereo))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowStereo, m_showAtomStereo );
	
	if (Known(has_showAtomEnhancedStereo))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowEnhancedStereo, m_showAtomEnhancedStereo );
	
	if (Known(has_showAtomNumber))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowAtomNumber, m_showAtomNumber );
	
	if (Known(has_showResidueID))
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowResidueID, m_showResidueID );
	
	if (Known(has_showBondQuery))
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowQuery, m_showBondQuery );
	
	if (Known(has_showBondRxn))
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowRxn, m_showBondRxn );
	
	if (Known(has_showBondStereo))
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowStereo, m_showBondStereo );

	if (Known(has_labelStyle))
	{
		if (m_labelStyle.family != UINT16(-1))
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleFont, m_labelStyle.family );
		if (m_labelStyle.size != 0)
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleSize, m_labelStyle.size / 20.0 );
		if (m_labelStyle.face != 0)
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleFace, m_labelStyle.face );
		if (m_labelStyle.color != UINT16(-1) && m_labelStyle.color != 3)
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleColor, m_labelStyle.color );
	}

	if (Known(has_captionStyle))
	{
		if (m_captionStyle.family != UINT16(-1))
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleFont, m_captionStyle.family );
		if (m_captionStyle.size != 0)
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleSize, m_captionStyle.size / 20.0 );
		if (m_captionStyle.face != 0)
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleFace, m_captionStyle.face );
		if (m_captionStyle.color != UINT16(-1) && m_captionStyle.color != 3)
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleColor, m_captionStyle.color );
	}

	if (Known(has_hashSpacing))
		CDXMLPutAttribute(sink_arg, kCDXProp_HashSpacing, CDXCoordinateToString(m_hashSpacing) );
	
	if (Known(has_marginWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_MarginWidth, CDXCoordinateToString(m_marginWidth) );
	
	if (Known(has_lineWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_LineWidth, CDXCoordinateToString(m_lineWidth) );
	
	if (Known(has_boldWidth))
		CDXMLPutAttribute(sink_arg, kCDXProp_BoldWidth, CDXCoordinateToString(m_boldWidth) );
	
	if (Known(has_bondLength))
		CDXMLPutAttribute(sink_arg, kCDXProp_BondLength, CDXCoordinateToString(m_bondLength) );
	
	if (Known(has_bondSpacing))
		CDXMLPutAttribute(sink_arg, kCDXProp_BondSpacing, m_bondSpacing / 10.0 );
	
	if (Known(has_bondSpacingAbs))
		CDXMLPutAttribute(sink_arg, kCDXProp_BondSpacingAbs, CDXCoordinateToString(m_bondSpacingAbs) );
	
	if (Known(has_chainAngle))
		CDXMLPutAttribute(sink_arg, kCDXProp_ChainAngle, CdxAngleToDegrees(m_chainAngle) );

	if (Known(has_labelJustification))
		CDXMLPutAttribute(sink_arg, kCDXProp_LabelJustification, m_labelJustification );

	if (Known(has_captionJustification))
		CDXMLPutAttribute(sink_arg, kCDXProp_CaptionJustification, m_captionJustification );

	if (Known(has_aminoAcidTermini))
		CDXMLPutAttribute(sink_arg, kCDXProp_AminoAcidTermini, m_aminoAcidTermini );

	if (Known(has_showSequenceTermini))
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceTermini, m_showSequenceTermini );

	if (Known(has_showSequenceBonds))
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceBonds, m_showSequenceBonds );

	if (Known(has_showSequenceUnlinkedBranches))
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceUnlinkedBranches, m_showSequenceUnlinkedBranches );

	if (Known(has_residueWrapCount))
		CDXMLPutAttribute(sink_arg, kCDXProp_ResidueWrapCount, m_residueWrapCount );

	if (Known(has_residueBlockCount))
		CDXMLPutAttribute(sink_arg, kCDXProp_ResidueBlockCount, m_residueBlockCount );

	if (Known(has_labelLineHeight))
	{
		if (m_labelLineHeight == kCDXLineHeight_Automatic)
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelLineHeight, std::string(kCDXML_auto) );
		else if (m_labelLineHeight != kCDXLineHeight_Variable)
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelLineHeight, CDXCoordinateToString(m_labelLineHeight) );
	}

	if (Known(has_captionLineHeight))
	{
		if (m_captionLineHeight == kCDXLineHeight_Variable)
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionLineHeight, std::string(kCDXML_variable) );
		else if (m_captionLineHeight != kCDXLineHeight_Automatic)
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionLineHeight, CDXCoordinateToString(m_captionLineHeight) );
	}

	if (m_lineType != kCDXLineType_Solid)
		CDXMLPutAttribute(sink_arg, kCDXProp_Line_Type, m_lineType );

	if (m_fillType != kCDXFillType_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Fill_Type, m_fillType );

	if (m_chemicalProperties != NULL && !m_chemicalProperties->empty())
		CDXMLPutAttribute(sink_arg, *m_chemicalProperties );
}

// CDXGraphicObject::XMLWriteContent
//
// Override this to write content between the start-tag and end-tag.
void CDXGraphicObject::XMLWriteContent(XMLDataSink &sink_arg) const
{
	if (m_propertiesRepresented != NULL && !m_propertiesRepresented->empty())
		for (std::vector<CDXPropRep>::const_iterator i = m_propertiesRepresented->begin(); i != m_propertiesRepresented->end();  ++i)
			sink_arg.os << *i;
}

bool CDXGraphicObject::XMLNeedToWriteContent() const
{
	return m_propertiesRepresented != NULL && !m_propertiesRepresented->empty();
}

bool CDXGraphicObject::HasHighlightColor() const
{ 
	return (m_flags2 & hasHighlightColor) != 0; 
}

void CDXGraphicObject::SetHighlightColor(CDXColorIndex colorIndex)
{ 
	m_flags2 |= hasHighlightColor; 
	highlightColor = colorIndex;
}

CDXColorIndex CDXGraphicObject::GetHighlightColor() const
{ 
	return highlightColor; 
}

void CDXGraphicObject::SetShowPerspective(bool inShowPerspective)
{
    showPerspective = inShowPerspective;
}

bool CDXGraphicObject::GetShowPerspective() const
{
    return showPerspective;
}
