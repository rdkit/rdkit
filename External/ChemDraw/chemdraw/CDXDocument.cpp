// CommonCS/LibCommon/Src/CDXDocument.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2007, CambridgeSoft Corp., All Rights Reserved

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
#include "CDXDocumentPropertyCollection.h"
#include "CDXMLNames.h"
#include "CDXUnicode.h"
#if defined (_DEBUG) || defined (CFW_)
#include "CDXUtils.h"	// for CdxIsValid()
#endif // defined (_DEBUG) || defined (CFW_)
#include "XMLPutEnum.h"
#include <assert.h>
#include <stdlib.h>	// for atoi and friends
#include <string.h>	// for strncmp and friends
#include <sstream>
#include <fstream>
#include <stdexcept>

extern XMLPutEnum<CDXTextJustification> sXMLTextJustificationText;
extern XMLPutEnum<CDXAminoAcidTermini> sXMLAminoAcidTerminiText;

XMLPutEnum<CDXRxnAutonumberStyle>::Value sXMLRxnAutonumberStyleValues[] = {
    {kCDXRxnAutonumberStyle_Alpha,		"Alpha"},
    {kCDXRxnAutonumberStyle_Arabic,	    "Arabic"},
    {kCDXRxnAutonumberStyle_Roman,	    "Roman"}
};

XMLPutEnum<CDXRxnAutonumberStyle> sXMLRxnAutonumberStyleText(sXMLRxnAutonumberStyleValues, sizeof sXMLRxnAutonumberStyleValues, kCDXRxnAutonumberStyle_Roman);
void XMLPut(XMLDataSink &sink, CDXRxnAutonumberStyle  v) {	sink.os << sXMLRxnAutonumberStyleText.asTokens(v); ; }


// ***********************
// ** class CDXDocument **
// ***********************
//
// Specialization of CDXObject for CDXDocument objects

CDXDocument::CDXDocument( CDXObjectID id_arg )
	:	CDXObject(kCDXObj_Document, id_arg)
	,	m_windowIsZoomed		(false)
	,	m_fractionalWidths		(false)
	,	m_interpretChemically	(false)
	,	m_showAtomQuery			(false)
	,	m_showAtomStereo		(false)
	,	m_showAtomEnhancedStereo(false)
	,	m_showAtomNumber		(false)
	,	m_showResidueID			(false)
	,	m_showBondQuery			(false)
	,	m_showBondRxn			(false)
	,	m_showBondStereo		(false)
	,	m_showTerminalCarbons	(false)
	,	m_showNonTerminalCarbons(false)
	,	m_hideImplicitHydrogens	(false)
	,	m_magnification			(0)
	,	m_hashSpacing			(0)
	,	m_marginWidth			(0)
	,	m_lineWidth				(0)
	,	m_boldWidth				(0)
	,	m_bondLength			(0)
	,	m_bondSpacing			(0)
	,	m_bondSpacingAbs		(0)
	,	m_bondSpacingType		(kCDXBondSpacingRelative)
	,	m_chainAngle			(0)
	,	m_labelLineHeight		(0)
	,	m_captionLineHeight		(0)
	,	m_labelJustification	(kCDXTextJustification_Auto)
	,	m_captionJustification	(kCDXTextJustification_Auto)
	,	m_aminoAcidTermini		(kCDXAminoAcidTerminiHOH)
	,	m_showSequenceTermini	(true)
	,	m_showSequenceBonds		(true)
	,   m_showSequenceUnlinkedBranches    (false)
	,	m_residueWrapCount		(40)
	,	m_residueBlockCount		(40)
	,	m_fgColor				(0)
	,	m_bgColor				(0)
	,	m_fgAlpha				(0)
	,	m_bgAlpha				(0)
    ,   m_macPrintRecord        {0}
	,	m_flags					(0)
	,	m_flags2				(0)
    ,   m_lineType(kCDXLineType_Solid)
    ,   m_fillType(kCDXFillType_Unspecified)
    ,   m_rxnAutonumberStyle(kCDXRxnAutonumberStyle_Roman)
    ,   m_rxnAutonumberStart    (1)
    ,   m_rxnAutonumberConditions(false)
    ,   m_rxnAutonumberFormat   ("(#)")
{
}

CDXDocument::~CDXDocument()
{
}

CDXObject*	CDXDocument::Clone() const
{
	return new CDXDocument (*this);
}

CDXDocument::CDXDocument (const CDXDocument& src) :
    CDXObject              (src),
    m_creationProgram      (src.m_creationProgram),
    m_name                 (src.m_name),
    m_comment              (src.m_comment),
    m_cartridgeData        (src.m_cartridgeData),
    m_boundingBox          (src.m_boundingBox),
    m_windowPosition       (src.m_windowPosition),
    m_windowSize           (src.m_windowSize),
    m_windowIsZoomed       (src.m_windowIsZoomed),
    m_fractionalWidths     (src.m_fractionalWidths),
    m_interpretChemically  (src.m_interpretChemically),
    m_showAtomQuery        (src.m_showAtomQuery),
    m_showAtomStereo       (src.m_showAtomStereo),
    m_showAtomEnhancedStereo(src.m_showAtomEnhancedStereo),
    m_showAtomNumber       (src.m_showAtomNumber),
    m_showResidueID        (src.m_showResidueID),
    m_showBondQuery        (src.m_showBondQuery),
    m_showBondRxn          (src.m_showBondRxn),
    m_showBondStereo       (src.m_showBondStereo),
    m_showTerminalCarbons  (src.m_showTerminalCarbons),
    m_showNonTerminalCarbons(src.m_showNonTerminalCarbons),
    m_hideImplicitHydrogens(src.m_hideImplicitHydrogens),
    m_magnification        (src.m_magnification),
    m_labelStyle           (src.m_labelStyle),
    m_captionStyle         (src.m_captionStyle),
    m_hashSpacing          (src.m_hashSpacing),
    m_marginWidth          (src.m_marginWidth),
    m_lineWidth            (src.m_lineWidth),
    m_boldWidth            (src.m_boldWidth),
    m_bondLength           (src.m_bondLength),
    m_bondSpacing          (src.m_bondSpacing),
    m_bondSpacingAbs       (src.m_bondSpacingAbs),
    m_bondSpacingType      (src.m_bondSpacingType),
    m_chainAngle           (src.m_chainAngle),
    m_labelLineHeight      (src.m_labelLineHeight),
    m_captionLineHeight    (src.m_captionLineHeight),
    m_labelJustification   (src.m_labelJustification),
    m_aminoAcidTermini     (src.m_aminoAcidTermini),
    m_showSequenceTermini  (src.m_showSequenceTermini),
    m_showSequenceBonds    (src.m_showSequenceBonds),
    m_showSequenceUnlinkedBranches(src.m_showSequenceUnlinkedBranches),
    m_captionJustification (src.m_captionJustification),
    m_residueWrapCount     (src.m_residueWrapCount),
    m_residueBlockCount    (src.m_residueBlockCount),
    m_printMargins         (src.m_printMargins),
    m_fixInPlaceExtent     (src.m_fixInPlaceExtent),
    m_fixInPlaceGap        (src.m_fixInPlaceGap),
    m_colorTable           (src.m_colorTable),
    m_fgColor              (src.m_fgColor),
    m_bgColor              (src.m_bgColor),
    m_fgAlpha              (src.m_fgAlpha),
    m_bgAlpha              (src.m_bgAlpha),
    m_fontTable            (src.m_fontTable),
    m_flags                (src.m_flags),
    m_flags2               (src.m_flags2),
    m_chemicalProperties   (src.m_chemicalProperties),
    m_fillType             (src.m_fillType),
    m_lineType             (src.m_lineType),
    m_rxnAutonumberStyle   (src.m_rxnAutonumberStyle),
    m_rxnAutonumberConditions(src.m_rxnAutonumberConditions),
    m_rxnAutonumberStart   (src.m_rxnAutonumberStart),
    m_rxnAutonumberFormat  (src.m_rxnAutonumberFormat),
    m_monomerRenderingStyle(src.m_monomerRenderingStyle)
{
#ifdef _DEBUG
	assert (CDXIsValid (this));
#endif
	memcpy(m_macPrintRecord,  src.m_macPrintRecord, sizeof m_macPrintRecord);
}

std::string CDXDocument::XMLObjectName() const
{
	return kCDXML_CDXML;
}

void CDXDocument::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_ForegroundColor:
		m_fgColor = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_fgColor;
		break;
	case kCDXProp_BackgroundColor:
		m_bgColor = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_bgColor;
		break;
	case kCDXProp_ForegroundAlpha:
		m_fgAlpha = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_fgAlpha;
		break;
	case kCDXProp_BackgroundAlpha:
		m_bgAlpha = UINT16(src_arg.GetUINT(size_arg));
		m_flags2 |= has_bgAlpha;
		break;
	case kCDXProp_ColorTable:
		m_colorTable.Read(src_arg, size_arg);
		break;
	case kCDXProp_FontTable:
		m_fontTable.Read(src_arg, size_arg);
		break;
	case kCDXProp_MacPrintInfo:
		if (size_arg != sizeof m_macPrintRecord)
			throw std::runtime_error("MacPrintInfo wrong size");
		src_arg.GetBytes(m_macPrintRecord, size_arg);
		m_flags |= has_printRecord;
		break;
	case kCDXProp_CreationProgram:
		m_creationProgram.Read(src_arg, size_arg);
		break;

	case kCDXProp_Name:
		m_name.Read(src_arg, size_arg);
		break;
	case kCDXProp_Comment:
		m_comment.Read(src_arg, size_arg);
		break;
	case kCDXProp_CartridgeData:
		m_cartridgeData.assign(src_arg.GetString(size_arg));
		break;
	case kCDXProp_BoundingBox:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate top    = src_arg.GetINT32();
		CDXCoordinate left   = src_arg.GetINT32();
		CDXCoordinate bottom = src_arg.GetINT32();
		CDXCoordinate right  = src_arg.GetINT32();
		BoundingBox(CDXRectangle(top, left, bottom, right));
		m_flags |= has_boundingBox;
		break;
		}
	case kCDXProp_Window_IsZoomed:
		m_windowIsZoomed = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags |= has_windowIsZoomed;
		break;
	case kCDXProp_Window_Position:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_windowPosition.y = src_arg.GetINT32();
		m_windowPosition.x = src_arg.GetINT32();
		m_flags |= has_windowPosition;
		break;
		}
	case kCDXProp_Window_Size:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_windowSize.y = src_arg.GetINT32();
		m_windowSize.x = src_arg.GetINT32();
		m_flags |= has_windowSize;
		break;
		}
	case kCDXProp_PrintMargins:
		{
		if (size_arg != 16)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		m_printMargins.top    = src_arg.GetINT32();
		m_printMargins.left   = src_arg.GetINT32();
		m_printMargins.bottom = src_arg.GetINT32();
		m_printMargins.right  = src_arg.GetINT32();
		m_flags |= has_printMargins;
		break;
		}
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
		if ((size_arg == 2) && (m_labelLineHeight != kCDXLineHeight_Automatic) && (m_labelLineHeight != kCDXLineHeight_Variable))
			m_labelLineHeight = (m_labelLineHeight * 0x10000) / 20; // pre-6.0 files have this in 1440'ths of an inch
		m_flags |= has_labelLineHeight;
		break;
	case kCDXProp_CaptionLineHeight:
		m_captionLineHeight = src_arg.GetUINT(size_arg);
		if ((size_arg == 2) && (m_captionLineHeight != kCDXLineHeight_Automatic) && (m_captionLineHeight != kCDXLineHeight_Variable))
			m_captionLineHeight = (m_captionLineHeight * 0x10000) / 20; // pre-6.0 files have this in 1440'ths of an inch
		m_flags |= has_captionLineHeight;
		break;
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
		m_flags2 |= has_showAtomEnhancedStereo;
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
		m_flags2 |= has_showTerminalCarbons;
		break;
	case kCDXProp_Atom_ShowNonTerminalCarbonLabels:
		m_showNonTerminalCarbons = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_showNonTerminalCarbons;
		break;
	case kCDXProp_Atom_HideImplicitHydrogens:
		m_hideImplicitHydrogens = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
		m_flags2 |= has_hideImplicitHydrogens;
		break;
	case kCDXProp_Magnification:
		m_magnification = src_arg.GetUINT(size_arg);
		m_flags |= has_magnification;
		break;
	case kCDXProp_FixInplaceExtent:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_fixInPlaceExtent = CDXPoint2D(x, y);
		m_flags |= has_fixInPlaceExtent;
		break;
		}
	case kCDXProp_FixInplaceGap:
		{
		if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_fixInPlaceGap = CDXPoint2D(x, y);
		m_flags |= has_fixInPlaceGap;
		break;
		}
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
			m_chemicalProperties.push_back(CDXChemProp(attribTag_arg, str));
		}
		break;
    case kCDXProp_RxnAutonumber_Conditions:
        m_rxnAutonumberConditions = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
        m_flags2 |= has_rxnAutonumberConditions;
        break;
    case kCDXProp_RxnAutonumber_Start:
        m_rxnAutonumberStart = src_arg.GetUINT(size_arg);
        m_flags2 |= has_rxnAutonumberStart;
        break;
    case kCDXProp_RxnAutonumber_Style:
        m_rxnAutonumberStyle = (CDXRxnAutonumberStyle)src_arg.GetINT(size_arg);
        m_flags2 |= has_rxnAutonumberStyle;
        break;
    case kCDXProp_RxnAutonumber_Format:
        {
            CDXString str;
            str.Read(src_arg, size_arg);
            m_rxnAutonumberFormat = str;
            m_flags2 |= has_rxnAutonumberFormat;
        }
        break;
    case kCDXProp_MonomerRenderingStyle:
        m_monomerRenderingStyle.Read(src_arg, size_arg);
        break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXDocument::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_CreationProgram:
		m_creationProgram = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_Name:
		m_name = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_Comment:
		m_comment = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_CartridgeData:
		m_cartridgeData = HexToBinary(attribValue_arg);
		break;
	case kCDXProp_BoundingBox:
		BoundingBox(StringToCDXRectangle(attribValue_arg));
		m_flags |= has_boundingBox;
		break;
	case kCDXProp_Window_IsZoomed:
		m_windowIsZoomed = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags |= has_windowIsZoomed;
		break;
	case kCDXProp_Window_Position:
		{
		std::istringstream is(attribValue_arg);
		is >> m_windowPosition.x;
		is >> m_windowPosition.y;
		m_windowPosition.x *= 0x10000;
		m_windowPosition.y *= 0x10000;
		m_flags |= has_windowPosition;
		break;
		}
	case kCDXProp_Window_Size:
		{
		std::istringstream is(attribValue_arg);
		is >> m_windowSize.x;
		is >> m_windowSize.y;
		m_windowSize.x *= 0x10000;
		m_windowSize.y *= 0x10000;
		m_flags |= has_windowSize;
		break;
		}
	case kCDXProp_PrintMargins:
		m_printMargins = StringToCDXRectangle(attribValue_arg);
		m_flags |= has_printMargins;
		break;
	case kCDXProp_ChainAngle:
		m_chainAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
		m_flags |= has_chainAngle;
		break;
	case kCDXProp_BondSpacing:
		m_bondSpacing = (int)(10. * atof(attribValue_arg.c_str()));
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
		m_flags2 |= has_showAtomEnhancedStereo;
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
		m_flags2 |= has_showTerminalCarbons;
		break;
	case kCDXProp_Atom_ShowNonTerminalCarbonLabels:
		m_showNonTerminalCarbons = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_showNonTerminalCarbons;
		break;
	case kCDXProp_Atom_HideImplicitHydrogens:
		m_hideImplicitHydrogens = (attribValue_arg == "yes") || (attribValue_arg == "true");
		m_flags2 |= has_hideImplicitHydrogens;
		break;
	case kCDXProp_Magnification:
		m_magnification = atoi(attribValue_arg.c_str());
		m_flags |= has_magnification;
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
	case kCDXProp_MacPrintInfo:
	{
		std::string s = HexToBinary(attribValue_arg);
		if (s.size() != sizeof m_macPrintRecord)
			throw std::runtime_error("print record must be exactly 120 bytes");
		std::copy(s.begin(), s.end(), m_macPrintRecord);
		m_flags |= has_printRecord;
		break;
	}
	case kCDXProp_ForegroundColor:
	{
		m_fgColor = atoi(attribValue_arg.c_str());
		m_flags2 |= has_fgColor;
		break;
	}
	case kCDXProp_BackgroundColor:
	{
		m_bgColor = atoi(attribValue_arg.c_str());
		m_flags2 |= has_bgColor;
		break;
	}
	case kCDXProp_ForegroundAlpha:
	{
		m_fgAlpha = atoi(attribValue_arg.c_str()) * 65535;
		m_flags2 |= has_fgAlpha;
		break;
	}
	case kCDXProp_BackgroundAlpha:
	{
		m_bgAlpha = atoi(attribValue_arg.c_str()) * 65535;
		m_flags2 |= has_bgAlpha;
		break;
	}
	case kCDXProp_FixInplaceExtent:
	{
		m_fixInPlaceExtent = StringToCDXPoint2D(attribValue_arg);
		m_flags |= has_fixInPlaceExtent;
		break;
	}
	case kCDXProp_FixInplaceGap:
	{
		m_fixInPlaceGap = StringToCDXPoint2D(attribValue_arg);
		m_flags |= has_fixInPlaceGap;
		break;
	}
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
			m_chemicalProperties.push_back(CDXChemProp(attribTag_arg, str));
			break;
		}
    case kCDXProp_RxnAutonumber_Conditions:
        {
            m_rxnAutonumberConditions = attribValue_arg == "yes" || attribValue_arg == "true";
            m_flags2 |= has_rxnAutonumberConditions;
            break;
        }
    case kCDXProp_RxnAutonumber_Start:
        m_rxnAutonumberStart = atoi(attribValue_arg.c_str());
        m_flags2 |= has_rxnAutonumberStart;
        break;
    case kCDXProp_RxnAutonumber_Style:
        m_rxnAutonumberStyle = sXMLRxnAutonumberStyleText.lookup(attribValue_arg);
        m_flags2 |= has_rxnAutonumberStyle;
        break;
    case kCDXProp_RxnAutonumber_Format:
        {
            CDXString str;
            str = CDXString(UnicodeToString(attribValue_arg));
            m_rxnAutonumberFormat = str;
            m_flags2 |= has_rxnAutonumberFormat;
            break;
        }
    case kCDXProp_MonomerRenderingStyle:
        m_monomerRenderingStyle = CDXString(UnicodeToString(attribValue_arg));;
        break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXDocument::WriteAttributesTo(CDXDataSink &sink_arg) const
{
    // Allow our objects to prepare any text using our font table before we write everything
    CdxWalker i(*this);
    const CDXObject *cdxObject = NULL;
    while (i.Next(cdxObject))
    {
        assert(cdxObject);
        
        // Have to cast away the constness because CdxWalker only supports const objects
        const_cast<CDXObject*>(cdxObject)->PrepareTextForWrite(m_fontTable);
    }

	CDXObject::WriteAttributesTo(sink_arg);

	if (m_creationProgram.length() != 0)
	{
		sink_arg.PutAttributeCDXString( kCDXProp_CreationProgram, m_creationProgram );
	}

	if (m_name.length() != 0)
	{
		sink_arg.PutAttributeCDXString( kCDXProp_Name, m_name );
	}

	if (m_comment.length() != 0)
	{
		sink_arg.PutAttributeCDXString( kCDXProp_Comment, m_comment );
	}

	if (!m_cartridgeData.empty())
	{
		sink_arg.PutTag(kCDXProp_CartridgeData);
		size_t cartridgeDataSize = m_cartridgeData.size();
		if (cartridgeDataSize <= 0xFFFE)
		{
			sink_arg.Put( UINT16(cartridgeDataSize) );
		}
		else
		{
			sink_arg.Put( UINT16(0xFFFF) );
			sink_arg.Put( UINT32(cartridgeDataSize) );
		}
		sink_arg.Put((INT8 *)m_cartridgeData.data(), cartridgeDataSize);
	}

	if (Known(has_boundingBox))
	{
		sink_arg.PutAttribute( kCDXProp_BoundingBox, BoundingBox() );
	}

	if (Known(has_windowPosition))
	{
		sink_arg.PutAttribute( kCDXProp_Window_Position, CDXPoint2D(m_windowPosition.x, m_windowPosition.y) );
	}

	if (Known(has_windowSize))
	{
		sink_arg.PutAttribute( kCDXProp_Window_Size, CDXPoint2D(m_windowSize.x, m_windowSize.y) );
	}

	if (Known(has_windowIsZoomed))
	{
		sink_arg.PutAttribute( kCDXProp_Window_IsZoomed );
	}
	
	if (Known(has_fractionalWidths))
	{
		sink_arg.PutAttribute( kCDXProp_FractionalWidths, (INT8)m_fractionalWidths );
	}
	
	if (Known(has_interpretChemically))
	{
		sink_arg.PutAttribute( kCDXProp_InterpretChemically, (INT8)m_interpretChemically );
	}
	
	if (Known(has_showAtomQuery))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowQuery, (INT8)m_showAtomQuery);
	}
	
	if (Known(has_showAtomStereo))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowStereo, (INT8)m_showAtomStereo );
	}
	
	if (Known(has_showAtomEnhancedStereo))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowEnhancedStereo, (INT8)m_showAtomEnhancedStereo );
	}
	
	if (Known(has_showAtomNumber))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowAtomNumber, (INT8)m_showAtomNumber );
	}
	
	if (Known(has_showResidueID))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowResidueID, (INT8)m_showResidueID );
	}
	
	if (Known(has_showBondQuery))
	{
		sink_arg.PutAttribute( kCDXProp_Bond_ShowQuery, (INT8)m_showBondQuery );
	}
	
	if (Known(has_showBondRxn))
	{
		sink_arg.PutAttribute( kCDXProp_Bond_ShowRxn, (INT8)m_showBondRxn );
	}
	
	if (Known(has_showBondStereo))
	{
		sink_arg.PutAttribute( kCDXProp_Bond_ShowStereo, (INT8)m_showBondStereo );
	}
	
	if (Known(has_showTerminalCarbons))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowTerminalCarbonLabels, (INT8)m_showTerminalCarbons );
	}
	
	if (Known(has_showNonTerminalCarbons))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_ShowNonTerminalCarbonLabels, (INT8)m_showNonTerminalCarbons );
	}
	
	if (Known(has_hideImplicitHydrogens))
	{
		sink_arg.PutAttribute( kCDXProp_Atom_HideImplicitHydrogens, (INT8)m_hideImplicitHydrogens );
	}
	
	if (Known(has_magnification))
	{
		sink_arg.PutAttribute( kCDXProp_Magnification, m_magnification );
	}
	
	if (m_labelStyle.family != UINT16(-1))
	{
		sink_arg.PutAttribute( kCDXProp_LabelStyle, m_labelStyle );
	}
	
	if (m_captionStyle.family != UINT16(-1))
	{
		sink_arg.PutAttribute( kCDXProp_CaptionStyle, m_captionStyle );
	}
	
	if (Known(has_hashSpacing))
	{
		sink_arg.PutAttribute( kCDXProp_HashSpacing, m_hashSpacing );
	}
	
	if (Known(has_marginWidth))
	{
		sink_arg.PutAttribute( kCDXProp_MarginWidth, m_marginWidth );
	}
	
	if (Known(has_lineWidth))
	{
		sink_arg.PutAttribute( kCDXProp_LineWidth, m_lineWidth );
	}
	
	if (Known(has_boldWidth))
	{
		sink_arg.PutAttribute( kCDXProp_BoldWidth, m_boldWidth );
	}
	
	if (Known(has_bondLength))
	{
		sink_arg.PutAttribute( kCDXProp_BondLength, m_bondLength );
	}
	
	if (Known(has_bondSpacing))
	{
		sink_arg.PutAttribute( kCDXProp_BondSpacing, m_bondSpacing );
	}
	
	if (Known(has_bondSpacingAbs))
	{
		sink_arg.PutAttribute( kCDXProp_BondSpacingAbs, m_bondSpacingAbs );
	}
	
	if (Known(has_chainAngle))
	{
		sink_arg.PutAttribute( kCDXProp_ChainAngle, m_chainAngle );
	}
	
	if (Known(has_labelJustification))
	{
		sink_arg.PutAttribute( kCDXProp_LabelJustification, (UINT8) m_labelJustification );
	}

	if (Known(has_captionJustification))
	{
		sink_arg.PutAttribute( kCDXProp_CaptionJustification, (UINT8) m_captionJustification );
	}

	if (Known(has_aminoAcidTermini))
	{
		sink_arg.PutAttribute( kCDXProp_AminoAcidTermini, (UINT8) m_aminoAcidTermini );
	}

	if (Known(has_showSequenceTermini))
	{
		sink_arg.PutAttribute( kCDXProp_ShowSequenceTermini, (INT8)m_showSequenceTermini );
	}
	
	if (Known(has_showSequenceBonds))
	{
		sink_arg.PutAttribute( kCDXProp_ShowSequenceBonds, (INT8)m_showSequenceBonds );
	}

	if (Known(has_showSequenceUnlinkedBranches))
	{
		sink_arg.PutAttribute( kCDXProp_ShowSequenceUnlinkedBranches, (INT8)m_showSequenceUnlinkedBranches );
	}
	
	if (Known(has_residueWrapCount))
	{
		sink_arg.PutAttribute( kCDXProp_ResidueWrapCount, (INT8)m_residueWrapCount );
	}
	
	if (Known(has_residueBlockCount))
	{
		sink_arg.PutAttribute( kCDXProp_ResidueBlockCount, (INT8)m_residueBlockCount );
	}
			
	if (Known(has_labelLineHeight))
	{
		// Put this out in 1440'ths of an inch, for compatibility with older ChemDraw versions
		if (m_labelLineHeight == kCDXLineHeight_Automatic)
		{
			sink_arg.PutAttribute( kCDXProp_LabelLineHeight, INT16(kCDXLineHeight_Automatic) );
		}
		else if (m_labelLineHeight != kCDXLineHeight_Variable)
		{
			sink_arg.PutAttribute( kCDXProp_LabelLineHeight, INT16(m_labelLineHeight * 20.0 / 65536.0 + 0.5) );
		}
	}

	if (Known(has_captionLineHeight))
	{
		// Put this out in 1440'ths of an inch, for compatibility with older ChemDraw versions
		if (m_captionLineHeight == kCDXLineHeight_Variable)
		{
			sink_arg.PutAttribute( kCDXProp_CaptionLineHeight, INT16(kCDXLineHeight_Variable) );
		}
		else if (m_captionLineHeight != kCDXLineHeight_Automatic)
		{
			sink_arg.PutAttribute( kCDXProp_CaptionLineHeight, INT16(m_captionLineHeight * 20.0 / 65536.0 + 0.5) );
		}
	}

	if (Known(has_printMargins))
	{
		sink_arg.PutAttribute( kCDXProp_PrintMargins, m_printMargins );
	}

	if (Known(has_fgColor))
	{
		sink_arg.PutAttribute( kCDXProp_ForegroundColor, m_fgColor );
	}

	if (Known(has_bgColor))
	{
		sink_arg.PutAttribute( kCDXProp_BackgroundColor, m_bgColor );
	}

	if (Known(has_fgAlpha))
	{
		sink_arg.PutAttribute(kCDXProp_ForegroundAlpha, m_fgAlpha);
	}

	if (Known(has_bgAlpha))
	{
		sink_arg.PutAttribute(kCDXProp_BackgroundAlpha, m_bgAlpha);
	}

	if (m_colorTable.NumColors() != 0)
	{
		m_colorTable.Write(sink_arg);
	}

	if (m_fontTable.NumFonts() != 0)
	{
		m_fontTable.Write(sink_arg);
	}

	if (Known(has_printRecord))
	{
		sink_arg.PutTag(kCDXProp_MacPrintInfo);
		sink_arg.Put(UINT16(sizeof m_macPrintRecord));
		sink_arg.Put((INT8 *)m_macPrintRecord, sizeof m_macPrintRecord);
	}

	if (Known(has_fixInPlaceExtent))
	{
		sink_arg.PutAttribute(kCDXProp_FixInplaceExtent, m_fixInPlaceExtent);
	}

	if (Known(has_fixInPlaceGap))
	{
		sink_arg.PutAttribute(kCDXProp_FixInplaceGap, m_fixInPlaceGap);
	}

	if (!m_chemicalProperties.empty())
	{
		sink_arg.PutAttribute( m_chemicalProperties );
	}
    
    if (Known(has_rxnAutonumberStart))
	{
        sink_arg.PutAttribute(kCDXProp_RxnAutonumber_Start, m_rxnAutonumberStart);
	}
    
    if (Known(has_rxnAutonumberConditions))
	{
        sink_arg.PutAttribute(kCDXProp_RxnAutonumber_Conditions, (INT8)m_rxnAutonumberConditions);
	}
    
    if (Known(has_rxnAutonumberStyle))
	{
        sink_arg.PutAttribute(kCDXProp_RxnAutonumber_Style, (UINT8)m_rxnAutonumberStyle);
	}
    
    if (Known(has_rxnAutonumberFormat))
	{
        sink_arg.PutAttributeCDXString(kCDXProp_RxnAutonumber_Format, m_rxnAutonumberFormat);
	}

    if (m_monomerRenderingStyle.length() != 0)
    {
        sink_arg.PutAttributeCDXString(kCDXProp_MonomerRenderingStyle, m_monomerRenderingStyle);
    }
}

void CDXDocument::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	if (m_creationProgram.length() != 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_CreationProgram, StringToUnicode("", m_creationProgram.str()) );
	}

	if (m_name.length() != 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Name, StringToUnicode("", m_name.str()) );
	}

	if (m_comment.length() != 0)
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Comment, StringToUnicode("", m_comment.str()) );
	}

	if (!m_cartridgeData.empty())
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_CartridgeData, BinaryToHex(m_cartridgeData) );
	}

	if (Known(has_boundingBox))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BoundingBox, BoundingBox() );
	}

	if (Known(has_windowPosition))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Window_Position, Point32((m_windowPosition.x * 0x10000), (m_windowPosition.y * 0x10000)) );
	}

	if (Known(has_windowSize))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Window_Size, Point32((m_windowSize.x * 0x10000), (m_windowSize.y * 0x10000)) );
	}

	if (Known(has_windowIsZoomed))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Window_IsZoomed, m_windowIsZoomed );
	}
	
	if (Known(has_fractionalWidths))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_FractionalWidths, m_fractionalWidths );
	}
	
	if (Known(has_interpretChemically))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_InterpretChemically, m_interpretChemically );
	}
	
	if (Known(has_showAtomQuery))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowQuery, m_showAtomQuery );
	}
	
	if (Known(has_showAtomStereo))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowStereo, m_showAtomStereo );
	}
	
	if (Known(has_showAtomEnhancedStereo))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowEnhancedStereo, m_showAtomEnhancedStereo );
	}
	
	if (Known(has_showAtomNumber))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowAtomNumber, m_showAtomNumber );
	}
	
	if (Known(has_showResidueID))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowResidueID, m_showResidueID );
	}
	
	if (Known(has_showBondQuery))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowQuery, m_showBondQuery );
	}
	
	if (Known(has_showBondRxn))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowRxn, m_showBondRxn );
	}
	
	if (Known(has_showBondStereo))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_ShowStereo, m_showBondStereo );
	}

	if (Known(has_showTerminalCarbons))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowTerminalCarbonLabels, m_showTerminalCarbons );
	}
	
	if (Known(has_showNonTerminalCarbons))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowNonTerminalCarbonLabels, m_showNonTerminalCarbons );
	}
	
	if (Known(has_hideImplicitHydrogens))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Atom_HideImplicitHydrogens, m_hideImplicitHydrogens );
	}
	
	if (Known(has_magnification))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Magnification, m_magnification );
	}
	
	if (Known(has_labelStyle))
	{
		if (m_labelStyle.family != UINT16(-1))
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleFont, m_labelStyle.family );
		}
		if (m_labelStyle.size != 0)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleSize, m_labelStyle.size / 20.0 );
		}
		if (m_labelStyle.face != 0)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleFace, m_labelStyle.face );
		}
		if ((m_labelStyle.color != UINT16(-1)) && (m_labelStyle.color != 3))
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelStyleColor, m_labelStyle.color );
		}
	}

	if (Known(has_captionStyle))
	{
		if (m_captionStyle.family != UINT16(-1))
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleFont, m_captionStyle.family );
		}
		if (m_captionStyle.size != 0)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleSize, m_captionStyle.size / 20.0 );
		}
		if (m_captionStyle.face != 0)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleFace, m_captionStyle.face );
		}
		if ((m_captionStyle.color != UINT16(-1)) && (m_captionStyle.color != 3))
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionStyleColor, m_captionStyle.color );
		}
	}

	if (Known(has_hashSpacing))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_HashSpacing, CDXCoordinateToString(m_hashSpacing) );
	}
	
	if (Known(has_marginWidth))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_MarginWidth, CDXCoordinateToString(m_marginWidth) );
	}
	
	if (Known(has_lineWidth))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_LineWidth, CDXCoordinateToString(m_lineWidth) );
	}
	
	if (Known(has_boldWidth))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BoldWidth, CDXCoordinateToString(m_boldWidth) );
	}
	
	if (Known(has_bondLength))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BondLength, CDXCoordinateToString(m_bondLength) );
	}
	
	if (Known(has_bondSpacing))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BondSpacing, m_bondSpacing / 10.0 );
	}
	
	if (Known(has_bondSpacingAbs))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BondSpacingAbs, CDXCoordinateToString(m_bondSpacingAbs) );
	}
	
	if (Known(has_chainAngle))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ChainAngle, CdxAngleToDegrees(m_chainAngle) );
	}
	
	if (Known(has_labelJustification))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_LabelJustification, m_labelJustification );
	}

	if (Known(has_captionJustification))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_CaptionJustification, m_captionJustification );
	}

	if (Known(has_aminoAcidTermini))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_AminoAcidTermini, m_aminoAcidTermini );
	}

	if (Known(has_showSequenceTermini))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceTermini, m_showSequenceTermini );
	}

	if (Known(has_showSequenceBonds))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceBonds, m_showSequenceBonds );
	}

	if (Known(has_showSequenceUnlinkedBranches))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ShowSequenceUnlinkedBranches, m_showSequenceUnlinkedBranches );
	}

	if (Known(has_residueWrapCount))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ResidueWrapCount, m_residueWrapCount );
	}

	if (Known(has_residueBlockCount))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ResidueBlockCount, m_residueBlockCount );
	}

	if (Known(has_labelLineHeight))
	{
		if (m_labelLineHeight == kCDXLineHeight_Automatic)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelLineHeight, std::string(kCDXML_auto) );
		}
		else if (m_labelLineHeight != kCDXLineHeight_Variable)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_LabelLineHeight, CDXCoordinateToString(m_labelLineHeight) );
		}
	}

	if (Known(has_captionLineHeight))
	{
		if (m_captionLineHeight == kCDXLineHeight_Variable)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionLineHeight, std::string(kCDXML_variable) );
		}
		else if (m_captionLineHeight != kCDXLineHeight_Automatic)
		{
			CDXMLPutAttribute(sink_arg, kCDXProp_CaptionLineHeight, CDXCoordinateToString(m_captionLineHeight) );
		}
	}

	if (Known(has_printMargins))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_PrintMargins, m_printMargins );
	}

	if (Known(has_printRecord))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_MacPrintInfo, BinaryToHex(std::string(m_macPrintRecord, sizeof m_macPrintRecord)) );
	}

	if (Known(has_fixInPlaceExtent))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_FixInplaceExtent, m_fixInPlaceExtent );
	}

	if (Known(has_fixInPlaceGap))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_FixInplaceGap, m_fixInPlaceGap );
	}

	if (!m_chemicalProperties.empty())
	{
		CDXMLPutAttribute(sink_arg, m_chemicalProperties );
	}

	if (Known(has_fgColor))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundColor, m_fgColor );
	}

	if (Known(has_bgColor))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BackgroundColor, m_bgColor );
	}

	if (Known(has_fgAlpha))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_ForegroundAlpha, CDXValueToString(m_fgAlpha / 65535.0, 4));
	}

	if (Known(has_bgAlpha))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_BackgroundAlpha, CDXValueToString(m_bgAlpha / 65535.0, 4));
	}

    if (Known(has_rxnAutonumberStart))
	{
        CDXMLPutAttribute(sink_arg, kCDXProp_RxnAutonumber_Start, m_rxnAutonumberStart);
	}
    
    if (Known(has_rxnAutonumberConditions))
	{
        CDXMLPutAttribute(sink_arg, kCDXProp_RxnAutonumber_Conditions, m_rxnAutonumberConditions);
	}
    
    if (Known(has_rxnAutonumberStyle))
	{
        CDXMLPutAttribute(sink_arg, kCDXProp_RxnAutonumber_Style, m_rxnAutonumberStyle);
	}
    
    if (Known(has_rxnAutonumberFormat))
	{
        CDXMLPutAttribute(sink_arg, kCDXProp_RxnAutonumber_Format, StringToUnicode("", m_rxnAutonumberFormat.str()) );
	}

    if (m_monomerRenderingStyle.length() != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_MonomerRenderingStyle, StringToUnicode("", m_monomerRenderingStyle.str()) );
    }
}

void CDXDocument::XMLWriteContent(XMLDataSink &sink_arg) const
{
	if (m_colorTable.NumColors() != 0)
	{
		m_colorTable.XMLWrite(sink_arg);
	}
	if (m_fontTable.NumFonts() != 0)
	{
		CDXPushFontTable(&m_fontTable);
		m_fontTable.XMLWrite(sink_arg);
	}
}

void CDXDocument::XMLWriteContentDone(XMLDataSink &sink_arg) const
{
	if (m_fontTable.NumFonts() != 0)
	{
		CDXPopFontTable(&m_fontTable);
	}
}

bool CDXDocument::XMLNeedToWriteContent() const
{
	return (m_colorTable.NumColors() != 0) || (m_fontTable.NumFonts() != 0);
}

void CDXDocument::FinishReading()
{
    // Allow our objects to finalise any text using our font table
    CdxWalker i(*this);
    const CDXObject *cdxObject = NULL;
    while (i.Next(cdxObject))
    {
        assert(cdxObject);
        
        // Have to cast away the constness because CdxWalker only supports const objects
        const_cast<CDXObject*>(cdxObject)->FinalizeTextAfterRead(m_fontTable);
    }
}

const CDXDocumentPropertyCollection* CDXDocument::GetDocumentProperties() const
{
    std::pair<CDXObjectsByTag::const_iterator, CDXObjectsByTag::const_iterator> range;

    const CDXObjectsByTag* contents = this->GetContents();
    if (!contents)
    {
        return NULL;
    }

    range = contents->equal_range(kCDXObj_DocumentProperties);

    if (range.first == range.second)
    {
        return NULL;
    }

    return static_cast<const CDXDocumentPropertyCollection*>(range.first->second);
}

CDXDocumentPropertyCollection* CDXDocument::GetDocumentProperties()
{
    std::pair<CDXObjectsByTag::const_iterator, CDXObjectsByTag::const_iterator> range;

    const CDXObjectsByTag* contents = this->GetContents();
    if (!contents)
    {
        this->AddChild(new CDXDocumentPropertyCollection());
        contents = this->GetContents();
    }

    range = contents->equal_range(kCDXObj_DocumentProperties);

    if (range.first == range.second)
    {
        this->AddChild(new CDXDocumentPropertyCollection());
        range = contents->equal_range(kCDXObj_DocumentProperties);
    }

    return static_cast<CDXDocumentPropertyCollection*>(range.first->second);
}

/*
 * Normalizes the text encoding for all objects in the document
 */
void CDXDocument::NormalizeTextEncoding()
{
    FinishReading();
}

// **********************************
// ***** CDXReadDocFromStorage  *****
// **********************************
// Construct a CDX object tree from the given input stream.
//
CDXDocument* CDXReadDocFromStorage (CDXDataSource &src, bool doThrow /* = true */)
{
	try
	{
		char hdrBuff[kCDX_HeaderLength];
		src.GetBytes(hdrBuff, kCDX_HeaderLength);

		// Verify header is a CDX header
		if (strncmp (hdrBuff, kCDX_HeaderString, kCDX_HeaderStringLen))
        {
        	throw std::runtime_error ("Not a CDX file");
        }

		if (hdrBuff[8] == '\1' && hdrBuff[9] == '\2' && hdrBuff[10] == '\3' && hdrBuff[11] == '\4')
        {
			src.SetBigEndian();
        }
        
		// Make a factory; just use the standard classes in CDXStdObjects.cpp
		CDXObjectFactory factory;
		CDXAddStandardFactories(factory);

		// Read one object and everything it contains
		CDXObject	*pTop = factory.ReadOneObject(src);
        CDXDocument	*pDoc = dynamic_cast<CDXDocument*>(pTop);
        if (!pDoc)
        {
            delete pTop;
        }
        
		return pDoc;
	} // try
	catch (const std::exception &)
	{
		if (doThrow)
			throw;
		return NULL;
	}
}
//-----------------------------------------
CDXDocument* CDXReadDocFromStorage (CDXDataSource &src, string& errMsg)
{
	errMsg.clear();
	try { CDXDocument* result = CDXReadDocFromStorage (src, true);  return result; }
	catch (std::bad_exception& ex)	{ errMsg = ex.what(); }
	catch (...)						{ errMsg = "Error reading CDX"; }
	return NULL;
}

// ********************************************************************
// **** CFW51_CDX   ChemDraw 5 couldn't "WriteTo" from the top doc ****
// ********************************************************************
//
void CFW51_CDX (const CDXObject *pObj, CDXDataSink &sink)
{
	pObj->WriteAttributesTo(sink);		// Attributes handled by the subclass
	if (pObj->m_contents != NULL)
	{
		pObj->m_contents->WriteTo(sink);		// The contained objects
	}
	sink.PutTag (kCDXProp_EndObject);
}

// *********************************
// ***** CDXWriteDocToStorage  *****
// *********************************
// Serialize a CDX object tree to the given output stream.
//
void CDXWriteDocToStorage (const CDXDocument* pDoc, CDXDataSink &sink)
{
	if (!pDoc)
	{
		return;	// empty file
	}
	CDXWriteHeader (sink);
	// ChemDraw 5.1 does not expect a leading document object.  We omit it, starting the
	// stream off with (after the header) the doc attributes.
	const bool fTargetIsCD5 = false;	// ChemDraw 7 has trouble too -heh 7/24/01; CD 8 is ok -heh 4/13/02
	if (fTargetIsCD5)
	{
		CFW51_CDX (pDoc, sink);
	}
	else	// CD 5 won't accept leading doc object; 6 does
	{
		pDoc->WriteTo (sink);
	}
	sink.PutTag(kCDXProp_EndObject);	// Won't hurt 8.0, and 7.0 really wants it.
}

// *******************************************************************
// **** ReadCdxDiskFile   Read a CDX object tree from a disk file ****
// *******************************************************************
//
CDXDocument*  ReadCdxDiskFile (const char *filename, bool doThrow)
{	// If doThrow is False, exceptions will be quashed, and a Null pointer returned.
	std::ifstream	is (filename, std::ios::binary);
	if (is.peek() == EOF)
	{
		return NULL;
	}
	CDXistream	source (is);

	return CDXReadDocFromStorage (source, doThrow);
}

// *********************************************************************
// **** WriteCdxDiskFile   Write a disk file from a CDX object tree ****
// *********************************************************************
//
 bool WriteCdxDiskFile (const CDXDocument *pDoc, const char *filename, std::string* pErrorMsg)
{
	std::ofstream	os (filename, std::ios::out | std::ios::trunc | std::ios::binary);
	if (!pDoc)
	{
		if (pErrorMsg)	*pErrorMsg = std::string ("Unable to access file ") + filename;
		{
		    return false;
		}
	}
	CDXostream	sink (os);
	CDXWriteDocToStorage (pDoc, sink);
	return true;
}

// *********************************************************************
// **** WriteCdxDiskFile   Write a disk file from a CDX object tree ****
// *********************************************************************
//
 bool WriteCdxmlDiskFile (const CDXDocument *pDoc, const char *filename, std::string* pErrorMsg)
{
	std::ofstream	dumpFile (filename, ios::binary);
	if (!pDoc)
	{
		if (pErrorMsg)	*pErrorMsg = std::string ("Unable to access file ") + filename;
		return false;
	}
	std::ostringstream os;
	os << kCDXML_HeaderString;
	XMLDataSink	ds (os);
	pDoc->XMLWrite (ds);
	dumpFile << os.str();
	return true;
}
