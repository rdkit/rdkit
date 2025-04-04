// CommonCS/LibCommon/Src/CDXNode.cpp
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
#include "CDXUtils.h"
#include "cs_auto_buffer.h"
#include <sstream>
#include <stdexcept>
using namespace std;

namespace
{
    /**
     *  Finds the attachment node in the specified groups or fragments
     *
     *  @param groupsOrFragments The groups or fragments in which to search for the attachment node
     *  @param attachmentIndex The index of the attachment
     *  @return The attachment node, or nullptr if not found
     */
    const CDXNode* FindAttachmentNodeInGroupsOrFragments(CDXObjectsRange groupsOrFragments, size_t attachmentIndex)
    {
        for (const auto groupOrFragment : groupsOrFragments)
        {
            auto groupOrFragmentObject = dynamic_cast<CDXGroupOrFragment*>(groupOrFragment.second);
            if (groupOrFragmentObject == nullptr)
            {
                continue;
            }

            const auto attachmentNode = groupOrFragmentObject->GetAttachmentNode(attachmentIndex);
            if (attachmentNode == nullptr)
            {
                continue;
            }

            return attachmentNode;
        }

        return nullptr;
    }
}

XMLPutEnum<CDXNodeType>::Value sXMLNodeTypeValues[] = {
	{kCDXNodeType_Unspecified,					"Unspecified"},
	{kCDXNodeType_Element,						"Element"},
	{kCDXNodeType_ElementList,					"ElementList"},
	{kCDXNodeType_ElementListNickname,			"ElementListNickname"},
	{kCDXNodeType_Nickname,						"Nickname"},
	{kCDXNodeType_Fragment,						"Fragment"},
	{kCDXNodeType_Formula,						"Formula"},
	{kCDXNodeType_GenericNickname,				"GenericNickname"},
	{kCDXNodeType_AnonymousAlternativeGroup,	"AnonymousAlternativeGroup"},
	{kCDXNodeType_NamedAlternativeGroup,		"NamedAlternativeGroup"},
	{kCDXNodeType_MultiAttachment,				"MultiAttachment"},
	{kCDXNodeType_VariableAttachment,			"VariableAttachment"},
	{kCDXNodeType_ExternalConnectionPoint,		"ExternalConnectionPoint"},
	{kCDXNodeType_LinkNode,						"LinkNode"},
    {kCDXNodeType_Monomer,                      "Monomer"}
};

XMLPutEnum<CDXNodeType> sXMLNodeType(sXMLNodeTypeValues, sizeof sXMLNodeTypeValues, kCDXNodeType_Unspecified);
void XMLPut(XMLDataSink &sink, CDXNodeType  v) {	sink.os << sXMLNodeType.lookup(v); }

XMLPutEnum<CDXRadical>::Value sXMLRadicalValues[] = {
	{kCDXRadical_None,		"None"},
	{kCDXRadical_Singlet,	"Singlet"},
	{kCDXRadical_Doublet,	"Doublet"},
	{kCDXRadical_Triplet,	"Triplet"}
};

XMLPutEnum<CDXRadical> sXMLRadical(sXMLRadicalValues, sizeof sXMLRadicalValues, kCDXRadical_None);
void XMLPut(XMLDataSink &sink, CDXRadical  v) {	sink.os << sXMLRadical.lookup(v); }

XMLPutEnum<CDXRingBondCount>::Value sXMLRingBondCountValues[] = {
	{kCDXRingBondCount_Unspecified,		"Unspecified"},
	{kCDXRingBondCount_NoRingBonds,		"NoRingBonds"},
	{kCDXRingBondCount_AsDrawn,			"AsDrawn"},
	{kCDXRingBondCount_SimpleRing,		"SimpleRing"},
	{kCDXRingBondCount_Fusion,			"Fusion"},
	{kCDXRingBondCount_SpiroOrHigher,	"SpiroOrHigher"}
};

XMLPutEnum<CDXRingBondCount> sXMLRingBondCount(sXMLRingBondCountValues, sizeof sXMLRingBondCountValues, kCDXRingBondCount_Unspecified);
void XMLPut(XMLDataSink &sink, CDXRingBondCount  v) {	sink.os << sXMLRingBondCount.lookup(v); }

XMLPutEnum<CDXUnsaturation>::Value sXMLUnsaturationValues[] = {
	{kCDXUnsaturation_Unspecified,		"Unspecified"},
	{kCDXUnsaturation_MustBeAbsent,		"MustBeAbsent"},
	{kCDXUnsaturation_MustBePresent,	"MustBePresent"}
};

XMLPutEnum<CDXUnsaturation> sXMLUnsaturation(sXMLUnsaturationValues, sizeof sXMLUnsaturationValues, kCDXUnsaturation_Unspecified);
void XMLPut(XMLDataSink &sink, CDXUnsaturation  v) {	sink.os << sXMLUnsaturation.lookup(v); }

XMLPutEnum<CDXReactionStereo>::Value sXMLReactionStereoValues[] = {
	{kCDXReactionStereo_Unspecified,	"Unspecified"},
	{kCDXReactionStereo_Inversion,		"Inversion"},
	{kCDXReactionStereo_Retention,		"Retention"}
};

XMLPutEnum<CDXReactionStereo> sXMLReactionStereo(sXMLReactionStereoValues, sizeof sXMLReactionStereoValues, kCDXReactionStereo_Unspecified);
void XMLPut(XMLDataSink &sink, CDXReactionStereo  v) {	sink.os << sXMLReactionStereo.lookup(v); }

XMLPutEnum<CDXTranslation>::Value sXMLTranslationValues[] = {
	{kCDXTranslation_Equal,		"Equal"},
	{kCDXTranslation_Broad,		"Broad"},
	{kCDXTranslation_Narrow,	"Narrow"},
	{kCDXTranslation_Any,		"Any"}
};

XMLPutEnum<CDXTranslation> sXMLTranslation(sXMLTranslationValues, sizeof sXMLTranslationValues, kCDXTranslation_Equal);
void XMLPut(XMLDataSink &sink, CDXTranslation  v) {	sink.os << sXMLTranslation.lookup(v); }

XMLPutEnum<CDXAbundance>::Value sXMLAbundanceValues[] = {
	{kCDXAbundance_Unspecified,	"Unspecified"},
	{kCDXAbundance_Any,			"Any"},
	{kCDXAbundance_Natural,		"Natural"},
	{kCDXAbundance_Enriched,	"Enriched"},
	{kCDXAbundance_Deficient,	"Deficient"},
	{kCDXAbundance_Nonnatural,	"Nonnatural"}
};

XMLPutEnum<CDXAbundance> sXMLAbundance(sXMLAbundanceValues, sizeof sXMLAbundanceValues, kCDXAbundance_Unspecified);
void XMLPut(XMLDataSink &sink, CDXAbundance  v) {	sink.os << sXMLAbundance.lookup(v); }

XMLPutEnum<CDXExternalConnectionType>::Value sXMLExternalConnectionTypeValues[] = {
	{kCDXExternalConnection_Unspecified,	"Unspecified"},
	{kCDXExternalConnection_Diamond,		"Diamond"},
	{kCDXExternalConnection_Star,			"Star"},
	{kCDXExternalConnection_PolymerBead,	"PolymerBead"},
	{kCDXExternalConnection_Wavy,			"Wavy"},
	{kCDXExternalConnection_Residue,		"Residue"},
	{kCDXExternalConnection_Peptide,		"Peptide"},
	{kCDXExternalConnection_DNA,			"DNA"},
	{kCDXExternalConnection_RNA,			"RNA"},
	{kCDXExternalConnection_Terminus,		"Terminus"},
	{kCDXExternalConnection_Sulfide,		"Sulfide"},
	{kCDXExternalConnection_Nucleotide,		"Nucleotide"},
    {kCDXExternalConnection_UnlinkedBranch, "UnlinkedBranch"}
};

XMLPutEnum<CDXExternalConnectionType> sXMLExternalConnectionType(sXMLExternalConnectionTypeValues, sizeof sXMLExternalConnectionTypeValues, kCDXExternalConnection_Unspecified);
void XMLPut(XMLDataSink &sink, CDXExternalConnectionType  v) {	sink.os << sXMLExternalConnectionType.lookup(v); }

XMLPutEnum<CDXAtomGeometry>::Value sXMLAtomGeometryValues[] = {
	{kCDXAtomGeometry_Unknown,				"Unknown"},
	{kCDXAtomGeometry_1Ligand,				"1"},
	{kCDXAtomGeometry_Linear,				"Linear"},
	{kCDXAtomGeometry_Bent,					"Bent"},
	{kCDXAtomGeometry_TrigonalPlanar,		"TrigonalPlanar"},
	{kCDXAtomGeometry_TrigonalPyramidal,	"TrigonalPyramidal"},
	{kCDXAtomGeometry_SquarePlanar,			"SquarePlanar"},
	{kCDXAtomGeometry_Tetrahedral,			"Tetrahedral"},
	{kCDXAtomGeometry_TrigonalBipyramidal,	"TrigonalBipyramidal"},
	{kCDXAtomGeometry_SquarePyramidal,		"SquarePyramidal"},
	{kCDXAtomGeometry_5Ligand,				"5"},
	{kCDXAtomGeometry_Octahedral,			"Octahedral"},
	{kCDXAtomGeometry_6Ligand,				"6"},
	{kCDXAtomGeometry_7Ligand,				"7"},
	{kCDXAtomGeometry_8Ligand,				"8"},
	{kCDXAtomGeometry_9Ligand,				"9"},
	{kCDXAtomGeometry_10Ligand,				"10"}
};

XMLPutEnum<CDXAtomGeometry> sXMLAtomGeometry(sXMLAtomGeometryValues, sizeof sXMLAtomGeometryValues, kCDXAtomGeometry_Unknown);
void XMLPut(XMLDataSink &sink, CDXAtomGeometry  v) {	sink.os << sXMLAtomGeometry.lookup(v); }

XMLPutEnum<CDXAtomCIPType>::Value sXMLAtomCIPTypeValues[] = {
    {kCDXCIPAtom_Undetermined, "U"},
    {kCDXCIPAtom_None,         "N"},
    {kCDXCIPAtom_R,            "R"},
    {kCDXCIPAtom_S,            "S"},
    {kCDXCIPAtom_r,            "r"},
    {kCDXCIPAtom_s,            "s"},
    {kCDXCIPAtom_Unspecified,  "u"},
    {kCDXCIPAtom_M,            "M"},
    {kCDXCIPAtom_P,            "P"},
    {kCDXCIPAtom_AtropisomerUnspecified, "a"},
    {kCDXCIPAtom_Allene_M, "-"},
    {kCDXCIPAtom_Allene_P, "+"},
    {kCDXCIPAtom_AlleneUnspecified, "?"}
};

XMLPutEnum<CDXAtomCIPType> sXMLAtomCIPType(sXMLAtomCIPTypeValues, sizeof sXMLAtomCIPTypeValues, kCDXCIPAtom_Undetermined);
void XMLPut(XMLDataSink &sink, CDXAtomCIPType  v) {	sink.os << sXMLAtomCIPType.lookup(v); }

XMLPutEnum<CDXLabelDisplay>::Value sXMLLabelDisplayValues[] = {
	{kCDXLabelDisplay_Auto,			"Auto"},
	{kCDXLabelDisplay_Left,			"Left"},
	{kCDXLabelDisplay_Center,		"Center"},
	{kCDXLabelDisplay_Right,		"Right"},
	{kCDXLabelDisplay_Above,		"Above"},
	{kCDXLabelDisplay_Below,		"Below"},
	{kCDXLabelDisplay_BestInitial,	"Best"}
};

XMLPutEnum<CDXLabelDisplay> sXMLLabelDisplay(sXMLLabelDisplayValues, sizeof sXMLLabelDisplayValues, kCDXLabelDisplay_Auto);
void XMLPut(XMLDataSink &sink, CDXLabelDisplay  v);
void XMLPut(XMLDataSink &sink, CDXLabelDisplay  v) {	sink.os << sXMLLabelDisplay.lookup(v); }

XMLPutEnum<CDXEnhancedStereoType>::Value sXMLEnhancedStereoValues[] = {
	{kCDXEnhancedStereo_Unspecified,"Unspecified"},
	{kCDXEnhancedStereo_None,		"None"},
	{kCDXEnhancedStereo_Absolute,	"Absolute"},
	{kCDXEnhancedStereo_Or,			"Or"},
	{kCDXEnhancedStereo_And,		"And"}
};

XMLPutEnum<CDXEnhancedStereoType> sXMLEnhancedStereo(sXMLEnhancedStereoValues, sizeof sXMLEnhancedStereoValues, kCDXEnhancedStereo_Unspecified);
void XMLPut(XMLDataSink &sink, CDXEnhancedStereoType  v) {	sink.os << sXMLEnhancedStereo.lookup(v); }

// *******************
// ** class CDXNode **
// *******************
//
// Specialization of CDXObject for CDXNode objects

CDXNode::CDXNode(CDXObjectID id_arg)
	: CDXGraphicObject		( kCDXObj_Node, id_arg )
	, m_charge				( 0 )
	, m_nodeType			( kCDXNodeType_Element )
	, m_elementNum			( UINT16(6) ) // Carbon
	, m_numHydrogens		( kNumHydrogenUnspecified )
	, m_isotope				( kCDXIsotope_Natural )
	, m_radical				( kCDXRadical_None )
	, m_substituentCode		( kCDXProp_EndObject )
	, m_substituentCount	( 0 )
	, m_ringBondCount		( kCDXRingBondCount_Unspecified )
	, m_unsaturation		( kCDXUnsaturation_Unspecified )
	, m_rxnStereo			( kCDXReactionStereo_Unspecified )
	, m_translation			( kCDXTranslation_Equal )
	, m_implicitHydrogens	( false )
	, m_abnormalValence		( false )
	, m_restrictRxnChange	( false )
    , m_isotopicAbundance	( kCDXAbundance_Unspecified )
    , m_externalConnectionType	( kCDXExternalConnection_Unspecified )
    , m_externalConnectionNum	( 0 )
	, m_geometry			( kCDXAtomGeometry_Unknown )
	, m_CIP					( kCDXCIPAtom_Undetermined )
	, m_altGroupID			( 0 )
	, m_labelDisplay		( kCDXLabelDisplay_Auto )
	, m_hStereo				( 0 )
	, m_bondOrdering		( NULL )
	, m_attachments			( NULL )
    , m_hydrogenBondAttachmentAtoms     ( {} )
    , m_hydrogenBonds       ( {} )
	, m_elementList			( NULL )
	, m_genericList			( NULL )
	, m_negativeList		( false )
	, m_linkCountLow		( 1 )
	, m_linkCountHigh		( 1 )
	, m_enhancedStereoType	( kCDXEnhancedStereo_Unspecified )
	, m_enhancedStereoGroupNum ( 0 )
	, m_needsClean			( false )
    , showAtomID            (false)
    , atomID                (0)
{
}

CDXNode::CDXNode(const CDXNode& src)
    : CDXGraphicObject(src)
    , m_charge(src.m_charge)
    , m_nodeType(src.m_nodeType)
    , m_elementNum(src.m_elementNum)
    , m_numHydrogens(src.m_numHydrogens)
    , m_isotope(src.m_isotope)
    , m_radical(src.m_radical)
    , m_substituentCode(src.m_substituentCode)
    , m_substituentCount(src.m_substituentCount)
    , m_ringBondCount(src.m_ringBondCount)
    , m_unsaturation(src.m_unsaturation)
    , m_rxnStereo(src.m_rxnStereo)
    , m_translation(src.m_translation)
    , m_implicitHydrogens(src.m_implicitHydrogens)
    , m_abnormalValence(src.m_abnormalValence)
    , m_restrictRxnChange(src.m_restrictRxnChange)
    , m_isotopicAbundance(src.m_isotopicAbundance)
    , m_externalConnectionType(src.m_externalConnectionType)
    , m_externalConnectionNum(src.m_externalConnectionNum)
    , m_geometry(src.m_geometry)
    , m_CIP(src.m_CIP)
    , m_altGroupID(src.m_altGroupID)
    , m_labelDisplay(src.m_labelDisplay)
    , m_hStereo(src.m_hStereo)
    , m_genericNickname(src.m_genericNickname)
    , m_atNum(src.m_atNum)
    , m_residueID(src.m_residueID)
    , m_bondOrdering((src.m_bondOrdering == NULL) ? NULL : new vector<CDXObjectID>(*src.m_bondOrdering))
    , m_attachments((src.m_attachments == NULL) ? NULL : new vector<CDXObjectID>(*src.m_attachments))
    , m_hydrogenBondAttachmentAtoms(src.m_hydrogenBondAttachmentAtoms)
    , m_hydrogenBonds(src.m_hydrogenBonds)
    ,	m_elementList		( (src.m_elementList == NULL) ? NULL : new vector<INT16>(*src.m_elementList) )
	,	m_genericList		( (src.m_genericList == NULL) ? NULL : new vector<string>(*src.m_genericList) )
	,	m_negativeList		(src.m_negativeList)
	,	m_linkCountLow		(src.m_linkCountLow)
	,	m_linkCountHigh		(src.m_linkCountHigh)
	,	m_enhancedStereoType(src.m_enhancedStereoType)
	,	m_enhancedStereoGroupNum(src.m_enhancedStereoGroupNum)
	,   m_needsClean		(src.m_needsClean)
    ,   showAtomID          (src.showAtomID)
    ,   atomID              (src.atomID)
{
}

CDXNode::~CDXNode()
{
	DeleteAndNull(m_bondOrdering);
	DeleteAndNull(m_attachments);
	DeleteAndNull(m_elementList);
	DeleteAndNull(m_genericList);
}

CDXObject*	CDXNode::Clone() const
{
	return new CDXNode (*this);
}

void CDXNode::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Node_Type:
		m_nodeType = (CDXNodeType) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Node_Element:
		m_elementNum = (UINT16) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Node_NeedsClean:
		this->m_needsClean = true;
		break;
	case kCDXProp_Atom_NumHydrogens:
		m_numHydrogens = (UINT16) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Atom_Charge:
		m_charge = src_arg.GetINT(size_arg);
		// If the charge is stored as an INT8, then it's an integer; convert it to a FIX24
		if (size_arg == 1)
			m_charge <<= 24;
		break;
	case kCDXProp_Atom_Radical:
		m_radical = (CDXRadical) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Atom_Isotope:
		m_isotope = (INT16) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_RestrictSubstituentsUpTo:
	case kCDXProp_Atom_RestrictSubstituentsExactly:
	case kCDXProp_Atom_RestrictFreeSites:
		m_substituentCode = attribTag_arg;
		m_substituentCount = (INT8) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_RestrictRingBondCount:
		m_ringBondCount = (CDXRingBondCount) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_RestrictUnsaturatedBonds:
		m_unsaturation = (CDXUnsaturation) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_RestrictRxnStereo:
		m_rxnStereo = static_cast<CDXReactionStereo>(src_arg.GetINT(size_arg));
		break;
	case kCDXProp_Atom_RestrictImplicitHydrogens:
		m_implicitHydrogens = true;
		break;
	case kCDXProp_Atom_AbnormalValence:
		m_abnormalValence = true;
		break;
	case kCDXProp_Atom_RestrictRxnChange:
		m_restrictRxnChange = true;
		break;
	case kCDXProp_Atom_Translation:
		m_translation = static_cast<CDXTranslation>(src_arg.GetINT(size_arg));
		break;
	case kCDXProp_Atom_IsotopicAbundance:
		m_isotopicAbundance = static_cast<CDXAbundance>(src_arg.GetINT(size_arg));
		break;
	case kCDXProp_Atom_ExternalConnectionType:
		m_externalConnectionType = static_cast<CDXExternalConnectionType>(src_arg.GetINT(size_arg));
		break;
    case kCDXProp_Atom_ExternalConnectionNum:
        m_externalConnectionNum = src_arg.GetINT(size_arg);
        break;
	case kCDXProp_Atom_Geometry:
		m_geometry = (CDXAtomGeometry) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_CIPStereochemistry:
		m_CIP = (CDXAtomCIPType) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_AltGroupID:
		m_altGroupID = src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_BondOrdering:
	{
		DeleteAndNull(m_bondOrdering);
		m_bondOrdering = new vector<CDXObjectID>;
		CDXReadItems(*m_bondOrdering, size_arg/sizeof(UINT32), src_arg);
		break;
	}
    case kCDXProp_Node_HydrogenBondAttachmentAtoms:
    {
        CDXReadItems(m_hydrogenBondAttachmentAtoms, size_arg/sizeof(UINT32), src_arg);
        break;
    }
    case kCDXProp_Node_HydrogenBonds:
    {
        CDXReadItems(m_hydrogenBonds, size_arg/sizeof(UINT32), src_arg);
        break;
    }
	case kCDXProp_Node_LabelDisplay:
		m_labelDisplay = (CDXLabelDisplay) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_GenericNickname:
		m_genericNickname.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_genericNickname[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_genericNickname [0];
			int	numSkipBytes = (int)m_genericNickname.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_genericNickname.erase(0, numSkipBytes);
		}
		break;
	case kCDXProp_Atom_AtomNumber:
	{
		if (size_arg > 0)
		{
			auto_buffer<char> buf((unsigned int)size_arg); 
			src_arg.GetBytes(buf.ptr(), size_arg); 

			m_atNum.assign(buf.ptr(), size_arg);

			// The count of styles is normally zero, but it might not be.
			if (size_arg > 1  &&  m_atNum[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
			{
				INT16 numStyles = m_atNum[0];
				int	numSkipBytes = (int)m_atNum.length();
				if (numSkipBytes > 2 + numStyles * 10)
					numSkipBytes = 2 + numStyles * 10;
				m_atNum.erase(0, numSkipBytes);
			}
		}
		break;
	}
	case kCDXProp_Atom_ResidueID:
	{
		if (size_arg > 0)
		{
			auto_buffer<char> buf((unsigned int)size_arg); 
			src_arg.GetBytes(buf.ptr(), size_arg); 

			m_residueID.assign(buf.ptr(), size_arg);

			// The count of styles is normally zero, but it might not be.
			if (size_arg > 1  &&  m_residueID[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
			{
				INT16 numStyles = m_residueID[0];
				int	numSkipBytes = (int)m_residueID.length();
				if (numSkipBytes > 2 + numStyles * 10)
					numSkipBytes = 2 + numStyles * 10;
				m_residueID.erase(0, numSkipBytes);
			}
		}
		break;
	}
	case kCDXProp_Node_Attachments:
	{
		size_t nItems = size_arg/sizeof(UINT32);
		if (size_arg % 4 == 2)	// a 2-byte count is present
		{
			nItems = src_arg.GetINT16();
			if (nItems * 4 + 2 != size_arg)
				throw std::runtime_error("kCDXProp_Node_Attachments: size does not match count");
		}
		DeleteAndNull(m_attachments);
		m_attachments = new vector<CDXObjectID>;
		CDXReadItems(*m_attachments, nItems, src_arg);
		break;
	}
	case kCDXProp_Atom_ElementList:
	{
		int nItems = src_arg.GetINT16();
		if (nItems < 0)
		{
			m_negativeList = true;
			nItems = -nItems;
		}
		DeleteAndNull(m_elementList);
		m_elementList = new vector<INT16>;
		CDXReadItems(*m_elementList, nItems, src_arg);
		break;
	}
	case kCDXProp_Atom_GenericList:
	{
		int nItems = src_arg.GetINT16();
		if (nItems < 0)
		{
			m_negativeList = true;
			nItems = -nItems;
		}
		string tempStr(src_arg.GetString(size_arg - 2));
		DeleteAndNull(m_genericList);
		m_genericList = new vector<string>;
		for (int i = 0; i < nItems && tempStr.size() > 4; ++i)
		{
			INT16 thisNameLen = *(INT16 *)tempStr.data();
			INT16 numStyles = *(INT16 *)(tempStr.data() + 2);
			if (10 * numStyles + 4 > tempStr.size())
				break;
			tempStr = tempStr.substr(10 * numStyles + 4);
			m_genericList->push_back(tempStr.substr(0, thisNameLen - 2 - numStyles * 10));
			tempStr = tempStr.substr(m_genericList->back().size());
		}
		break;
	}
	case kCDXProp_Atom_HDot:
	case kCDXProp_Atom_HDash:
		if (size_arg == 0)
			m_hStereo = attribTag_arg;
		else if (src_arg.GetINT(size_arg))			// Written out as a boolean value - only set the member if it's non-zero
			m_hStereo = attribTag_arg;
		break;
	case kCDXProp_Atom_LinkCountLow:
		m_linkCountLow = (UINT16) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Atom_LinkCountHigh:
		m_linkCountHigh = (UINT16) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Atom_EnhancedStereoType:
		m_enhancedStereoType = (CDXEnhancedStereoType) src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Atom_EnhancedStereoGroupNum:
		m_enhancedStereoGroupNum = (INT8) src_arg.GetINT(size_arg);
		break;

    case kCDXProp_Atom_ShowAtomID:
        showAtomID = src_arg.GetINT(size_arg) != 0;
        break;

    case kCDXProp_Atom_AtomID:
        atomID = src_arg.GetUINT32();
        break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXNode::XMLStoreAttribute(CDXTag attribTag_arg, const string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Node_Type:
		m_nodeType = sXMLNodeType.lookup(attribValue_arg);
		break;
	case kCDXProp_Node_Element:
		m_elementNum = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Node_NeedsClean:
		m_needsClean = attribValue_arg == "yes" || attribValue_arg == "true";
		break;
	case kCDXProp_Atom_NumHydrogens:
		m_numHydrogens = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_Charge:
		m_charge = (INT32)(atof(attribValue_arg.c_str()) * 0x1000000);
		break;
	case kCDXProp_Atom_Radical:
		m_radical = sXMLRadical.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_Isotope:
		m_isotope = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_RestrictSubstituentsUpTo:
	case kCDXProp_Atom_RestrictSubstituentsExactly:
	case kCDXProp_Atom_RestrictFreeSites:
		m_substituentCode = attribTag_arg;
		m_substituentCount = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_RestrictRingBondCount:
		m_ringBondCount = sXMLRingBondCount.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_RestrictUnsaturatedBonds:
		m_unsaturation = sXMLUnsaturation.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_RestrictRxnStereo:
		m_rxnStereo = sXMLReactionStereo.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_RestrictImplicitHydrogens:
		m_implicitHydrogens = attribValue_arg == "yes" || attribValue_arg == "true";
		break;
	case kCDXProp_Atom_AbnormalValence:
		m_abnormalValence = attribValue_arg == "yes" || attribValue_arg == "true";
		break;
	case kCDXProp_Atom_RestrictRxnChange:
		m_restrictRxnChange = attribValue_arg == "yes" || attribValue_arg == "true";
		break;
	case kCDXProp_Atom_Translation:
		m_translation = sXMLTranslation.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_IsotopicAbundance:
		m_isotopicAbundance = sXMLAbundance.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_ExternalConnectionType:
		m_externalConnectionType = sXMLExternalConnectionType.lookup(attribValue_arg);
        break;
    case kCDXProp_Atom_ExternalConnectionNum:
        m_externalConnectionNum = atoi(attribValue_arg.c_str());
        break;
	case kCDXProp_Atom_Geometry:
		m_geometry = sXMLAtomGeometry.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_CIPStereochemistry:
		m_CIP = sXMLAtomCIPType.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_AltGroupID:
		m_altGroupID = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_BondOrdering:
	{
		istringstream is(attribValue_arg);
		DeleteAndNull(m_bondOrdering);
		m_bondOrdering = new vector<CDXObjectID>;
		copy(istream_iterator<CDXObjectID>(is),
			 istream_iterator<CDXObjectID>(),
			 back_inserter(*m_bondOrdering));
		break;
	}
    case kCDXProp_Node_HydrogenBondAttachmentAtoms:
    {
        istringstream is(attribValue_arg);
        copy(istream_iterator<CDXObjectID>(is),
            istream_iterator<CDXObjectID>(),
            back_inserter(m_hydrogenBondAttachmentAtoms));
        break;
    }
    case kCDXProp_Node_HydrogenBonds:
    {
        istringstream is(attribValue_arg);
        copy(istream_iterator<CDXObjectID>(is),
            istream_iterator<CDXObjectID>(),
            back_inserter(m_hydrogenBonds));
        break;
    }

	case kCDXProp_Node_LabelDisplay:
		m_labelDisplay = sXMLLabelDisplay.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_GenericNickname:
		m_genericNickname = attribValue_arg;
		break;
	case kCDXProp_Atom_AtomNumber:
		m_atNum = attribValue_arg;
		break;
	case kCDXProp_Atom_ResidueID:
		m_residueID = attribValue_arg;
		break;
	case kCDXProp_Node_Attachments:
	{
		istringstream is(attribValue_arg);
		DeleteAndNull(m_attachments);
		m_attachments = new vector<CDXObjectID>;
		copy(istream_iterator<CDXObjectID>(is),
			 istream_iterator<CDXObjectID>(),
			 back_inserter(*m_attachments));
		break;
	}
	case kCDXProp_Atom_ElementList:
	{
		int startPos = 0;
		if (attribValue_arg.substr(0, 3) == "NOT")
		{
			m_negativeList = true;
			startPos = 3;
		}
		istringstream is(attribValue_arg.substr(startPos, string::npos));
		DeleteAndNull(m_elementList);
		m_elementList = new vector<INT16>;
		copy(istream_iterator<INT16>(is),
			 istream_iterator<INT16>(),
			 back_inserter(*m_elementList));
		break;
	}
	case kCDXProp_Atom_GenericList:
	{
		int startPos = 0;
		if (attribValue_arg.substr(0, 3) == "NOT")
		{
			m_negativeList = true;
			startPos = 3;
		}
		istringstream is(attribValue_arg.substr(startPos, string::npos));
		DeleteAndNull(m_genericList);
		m_genericList = new vector<string>;
		copy(istream_iterator<string>(is),
			 istream_iterator<string>(),
			 back_inserter(*m_genericList));
		break;
	}
	case kCDXProp_Atom_HDot:
	case kCDXProp_Atom_HDash:
		if (attribValue_arg == "yes" || attribValue_arg == "true")
			m_hStereo = attribTag_arg;
		break;
	case kCDXProp_Atom_LinkCountLow:
		m_linkCountLow = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_LinkCountHigh:
		m_linkCountHigh = atoi(attribValue_arg.c_str());
		break;
	case kCDXProp_Atom_EnhancedStereoType:
		m_enhancedStereoType = sXMLEnhancedStereo.lookup(attribValue_arg);
		break;
	case kCDXProp_Atom_EnhancedStereoGroupNum:
		m_enhancedStereoGroupNum = atoi(attribValue_arg.c_str());
		break;

    case kCDXProp_Atom_ShowAtomID:
        showAtomID = (attribValue_arg == "yes") || (attribValue_arg == "true");
        break;

    case kCDXProp_Atom_AtomID:
        atomID = atoi(attribValue_arg.c_str());
        break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXNode::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if ( m_nodeType != kCDXNodeType_Element )
		sink_arg.PutAttribute( kCDXProp_Node_Type, UINT16(m_nodeType) );

	if ( (m_nodeType == kCDXNodeType_Element || m_nodeType == kCDXNodeType_LinkNode) && m_elementNum != 6 )
		sink_arg.PutAttribute( kCDXProp_Node_Element, UINT16(m_elementNum) );

	if ( m_numHydrogens != UINT16(-1) )
		sink_arg.PutAttribute( kCDXProp_Atom_NumHydrogens, UINT16(m_numHydrogens) );

	if ( m_charge != 0 )
	{
		if ((m_charge & 0x00FFFFFF) == 0)
			sink_arg.PutAttribute( kCDXProp_Atom_Charge, INT8(m_charge >> 24) );
		else
			sink_arg.PutAttribute( kCDXProp_Atom_Charge, UINT32(m_charge) );
	}

	if ( m_radical != 0 )
		sink_arg.PutAttribute( kCDXProp_Atom_Radical, UINT8(m_radical) );

	if ( m_isotope != 0 )
		sink_arg.PutAttribute( kCDXProp_Atom_Isotope, INT16(m_isotope) );

	if ( m_substituentCode != kCDXProp_EndObject )
		sink_arg.PutAttribute( m_substituentCode, INT8(m_substituentCount) );

	if (m_ringBondCount != kCDXRingBondCount_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Atom_RestrictRingBondCount, INT8(m_ringBondCount) );

	if (m_unsaturation != kCDXUnsaturation_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Atom_RestrictUnsaturatedBonds, (INT8) m_unsaturation);

	if (m_rxnStereo != kCDXReactionStereo_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Atom_RestrictRxnStereo, (INT8) m_rxnStereo);

	if (m_translation != kCDXTranslation_Equal)
		sink_arg.PutAttribute( kCDXProp_Atom_Translation, (INT8) m_translation);

	if (m_isotopicAbundance != kCDXAbundance_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Atom_IsotopicAbundance, (INT8) m_isotopicAbundance);

	if (m_implicitHydrogens)
		sink_arg.PutAttribute(kCDXProp_Atom_RestrictImplicitHydrogens);

	if (m_externalConnectionType)
    {
		sink_arg.PutAttribute(kCDXProp_Atom_ExternalConnectionType, (INT8) m_externalConnectionType);
    }
    
    if ((kCDXNodeType_ExternalConnectionPoint == m_nodeType)
        && m_externalConnectionNum > 0)
    {
        sink_arg.PutAttribute(kCDXProp_Atom_ExternalConnectionNum, (INT8) m_externalConnectionNum);
    }
    
	if (m_abnormalValence)
		sink_arg.PutAttribute(kCDXProp_Atom_AbnormalValence);

	if (m_restrictRxnChange)
		sink_arg.PutAttribute(kCDXProp_Atom_RestrictRxnChange);

	if (m_needsClean)
		sink_arg.PutAttribute(kCDXProp_Node_NeedsClean);

	if (m_geometry != kCDXAtomGeometry_Unknown)
		sink_arg.PutAttribute(kCDXProp_Atom_Geometry, INT8(m_geometry));

	if (m_CIP != kCDXCIPAtom_Undetermined)
		sink_arg.PutAttribute(kCDXProp_Atom_CIPStereochemistry, INT8(m_CIP));

	if (m_altGroupID != 0)
		sink_arg.PutAttribute(kCDXProp_Atom_AltGroupID, m_altGroupID);

	if (m_bondOrdering != NULL && !m_bondOrdering->empty())
	{
		sink_arg.PutTag( kCDXProp_Atom_BondOrdering );
		sink_arg.Put( UINT16( 4 * m_bondOrdering->size() ) );
		Put( sink_arg, *m_bondOrdering );
	}

    if (!m_hydrogenBondAttachmentAtoms.empty())
    {
        sink_arg.PutTag(kCDXProp_Node_HydrogenBondAttachmentAtoms);
        sink_arg.Put(UINT16(4 * m_hydrogenBondAttachmentAtoms.size()));
        Put(sink_arg, m_hydrogenBondAttachmentAtoms);
    }

    if (!m_hydrogenBonds.empty())
    {
        sink_arg.PutTag(kCDXProp_Node_HydrogenBonds);
        sink_arg.Put(UINT16(4 * m_hydrogenBonds.size()));
        Put(sink_arg, m_hydrogenBonds);
    }

	if (m_labelDisplay != kCDXLabelDisplay_Auto)
		sink_arg.PutAttribute(kCDXProp_Node_LabelDisplay, INT8(m_labelDisplay));

	if (!m_genericNickname.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_Atom_GenericNickname );
        CDXPut(sink_arg, m_genericNickname);
	}

	if (!m_atNum.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_Atom_AtomNumber );
        CDXPut(sink_arg, m_atNum);
	}

	if (!m_residueID.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_Atom_ResidueID );
        CDXPut(sink_arg, m_residueID);
	}

	if (m_attachments && !m_attachments->empty())	
	{
		sink_arg.PutTag( kCDXProp_Node_Attachments );
		sink_arg.Put( UINT16( 2 + 4 * m_attachments->size() ) );
		sink_arg.Put( UINT16( m_attachments->size() ) );
		Put( sink_arg, *m_attachments );
	}

	if (m_elementList != NULL && !m_elementList->empty())	
	{
		sink_arg.PutTag( kCDXProp_Atom_ElementList );
		sink_arg.Put( UINT16( 2 + 2 * m_elementList->size() ) );
#if 0	// The ff. line generates a compiler warning about negating an unsigned value
		sink_arg.Put( UINT16( m_negativeList ? -m_elementList->size() : m_elementList->size() ) );
#else
		const INT16	val = (INT16)m_elementList->size();
		sink_arg.Put( UINT16( m_negativeList ? -val : val));
#endif
		Put( sink_arg, *m_elementList );
	}

	if (m_genericList != NULL && !m_genericList->empty())	
	{
		sink_arg.PutTag( kCDXProp_Atom_GenericList );
		UINT16 argSize = 2;
		vector<string>::const_iterator s;
		for (s = m_genericList->begin(); s != m_genericList->end(); ++s)
			argSize += 4 + s->size();
		sink_arg.Put( argSize );
#if 0	// The ff. line generates a compiler warning about negating an unsigned value
		sink_arg.Put( UINT16( m_negativeList ? -m_genericList->size() : m_genericList->size() ) );
#else
		const INT16	val = (INT16)m_genericList->size();
		sink_arg.Put( UINT16( m_negativeList ? -val : val));
#endif
		for (s = m_genericList->begin(); s != m_genericList->end(); ++s)
		{
			sink_arg.Put( UINT16( s->size() + 2 ) );
			sink_arg.Put( UINT16 (0) );	// the styles
			sink_arg.Put( (INT8 *)s->data(), s->size() );
		}
	}

	if (m_hStereo != 0)
		sink_arg.PutAttribute(m_hStereo);

	if ( m_linkCountLow != 1 )
		sink_arg.PutAttribute( kCDXProp_Atom_LinkCountLow, UINT16(m_linkCountLow) );

	if ( m_linkCountHigh != 1 )
		sink_arg.PutAttribute( kCDXProp_Atom_LinkCountHigh, UINT16(m_linkCountHigh) );

	if ( m_enhancedStereoType != kCDXEnhancedStereo_Unspecified )
		sink_arg.PutAttribute( kCDXProp_Atom_EnhancedStereoType, INT8(m_enhancedStereoType) );

	if ( (m_enhancedStereoType == kCDXEnhancedStereo_Or || m_enhancedStereoType == kCDXEnhancedStereo_And) && m_enhancedStereoGroupNum > 0 )
		sink_arg.PutAttribute( kCDXProp_Atom_EnhancedStereoGroupNum, INT16(m_enhancedStereoGroupNum) );

    if (showAtomID)
    {
        sink_arg.PutAttribute(kCDXProp_Atom_ShowAtomID, (UINT8)showAtomID);
    }

    if (atomID > 0)
    {
        sink_arg.PutAttribute(kCDXProp_Atom_AtomID, atomID);
    }
}

string CDXNode::XMLObjectName() const
{
	return kCDXML_node;
}

void CDXNode::XMLWriteAttributes(XMLDataSink& sink_arg) const
{
    CDXGraphicObject::XMLWriteAttributes(sink_arg);

    if (m_nodeType != kCDXNodeType_Element)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_Type, m_nodeType);
    }

    if (((m_nodeType == kCDXNodeType_Element) || (m_nodeType == kCDXNodeType_LinkNode)) && (m_elementNum != 6))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_Element, m_elementNum);
    }

    if (m_numHydrogens != UINT16(-1))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_NumHydrogens, m_numHydrogens);
    }

    if (m_charge != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_Charge, m_charge / double(0x1000000));
    }

    if (m_radical != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_Radical, m_radical);
    }

    if (m_isotope != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_Isotope, m_isotope);
    }

    if (m_substituentCode != kCDXProp_EndObject)
    {
        CDXMLPutAttribute(sink_arg, m_substituentCode, m_substituentCount);
    }

    if (m_ringBondCount != kCDXRingBondCount_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_RestrictRingBondCount, m_ringBondCount);
    }

    if (m_unsaturation != kCDXUnsaturation_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_RestrictUnsaturatedBonds, m_unsaturation);
    }

    if (m_rxnStereo != kCDXReactionStereo_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_RestrictRxnStereo, m_rxnStereo);
    }

    if (m_translation != kCDXTranslation_Equal)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_Translation, m_translation);
    }

    if (m_isotopicAbundance != kCDXAbundance_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_IsotopicAbundance, m_isotopicAbundance);
    }

    if (m_externalConnectionType)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ExternalConnectionType, m_externalConnectionType);
    }

    if ((kCDXNodeType_ExternalConnectionPoint == m_nodeType)
        && m_externalConnectionNum > 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ExternalConnectionNum, m_externalConnectionNum);
    }

    if (m_implicitHydrogens)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_RestrictImplicitHydrogens);
    }

    if (m_abnormalValence)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_AbnormalValence);
    }

    if (m_restrictRxnChange)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_RestrictRxnChange);
    }

    if (m_needsClean)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_NeedsClean);
    }

    if (m_geometry != kCDXAtomGeometry_Unknown)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_Geometry, m_geometry);
    }

    if (m_CIP != kCDXCIPAtom_Undetermined)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_CIPStereochemistry, m_CIP);
    }

    if (m_altGroupID != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_AltGroupID, m_altGroupID);
    }

    if (m_bondOrdering != NULL && !m_bondOrdering->empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_BondOrdering, *m_bondOrdering);
    }

    if (m_labelDisplay != kCDXLabelDisplay_Auto)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_LabelDisplay, m_labelDisplay);
    }

    if (!m_genericNickname.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_GenericNickname, m_genericNickname);
    }

    if (!m_atNum.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_AtomNumber, m_atNum);
    }

    if (!m_residueID.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ResidueID, m_residueID);
    }

    if (m_attachments && !m_attachments->empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_Attachments, *m_attachments);
    }

    if (!m_hydrogenBondAttachmentAtoms.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_HydrogenBondAttachmentAtoms, m_hydrogenBondAttachmentAtoms);
    }

    if (!m_hydrogenBonds.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Node_HydrogenBonds, m_hydrogenBonds);
    }

	if (m_elementList != NULL && !m_elementList->empty())	
	{
		// Do this one the hard way, so we can slip in the NOT if we need to
		sink_arg.os << " " << CDXMLAttributeName(kCDXProp_Atom_ElementList) << "=\"";
		if (m_negativeList)
			sink_arg.os << "NOT ";
		XMLPut(sink_arg, *m_elementList);
		sink_arg.os << "\"" << GetTextEOL();
	}

	if (m_genericList != NULL && !m_genericList->empty())	
	{
		// Do this one the hard way, so we can slip in the NOT if we need to
		sink_arg.os << " " << CDXMLAttributeName(kCDXProp_Atom_GenericList) << "=\"";
		if (m_negativeList)
			sink_arg.os << "NOT ";
		XMLPut(sink_arg, *m_genericList);
		sink_arg.os << "\"" << GetTextEOL();
	}

    if (m_hStereo != 0)
    {
        CDXMLPutAttribute(sink_arg, m_hStereo);
    }

	if ( m_linkCountLow != 1 )
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_LinkCountLow, m_linkCountLow);
    }

    if (m_linkCountHigh != 1)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_LinkCountHigh, m_linkCountHigh);
    }

    if (m_enhancedStereoType != kCDXEnhancedStereo_Unspecified)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_EnhancedStereoType, m_enhancedStereoType);
    }

    if (((m_enhancedStereoType == kCDXEnhancedStereo_Or) || (m_enhancedStereoType == kCDXEnhancedStereo_And)) && (m_enhancedStereoGroupNum > 0))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_EnhancedStereoGroupNum, m_enhancedStereoGroupNum);
    }

    if (showAtomID)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_ShowAtomID);
    }

    if (atomID > 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Atom_AtomID, atomID);
    }
}

/**
 *  Gets the attachment node for the specified bond
 *
 *  @param bondID The ID of the bond to get the attachment node for
 *  @return The attachment node, or nullptr if not found
 */
const CDXNode* CDXNode::GetAttachmentNodeForBond(CDXObjectID bondID) const
{
    // If this node contains one or more fragments, dig into them to find 
    // the node within it that the specified bond attaches to
    CDXObjectsRange fragments = ContainedObjects(kCDXObj_Fragment);
    CDXObjectsRange groups = ContainedObjects(kCDXObj_Group);

    // No contained fragments or groups
    if ((fragments.size() == 0) && (groups.size() == 0))
    {
        return this;
    }

    // If there are no bonds, just return
    if (m_alphaBonds.size() == 0)
    {
        return nullptr;
    }

    size_t attachmentIndex = 0;

    // If there's just one bond, we don't need a bond ordering
    if (m_alphaBonds.size() > 1)
    {
        if (m_bondOrdering == nullptr)
        {
            return nullptr;
        }

        // Find which element in the bond ordering contains this bond ID
        const auto bondIterator = std::find(m_bondOrdering->begin(), m_bondOrdering->end(), bondID);
        if (bondIterator == m_bondOrdering->end())
        {
            return nullptr;
        }

        attachmentIndex = bondIterator - m_bondOrdering->begin();
    }

    const CDXNode* attachmentNode = nullptr;

    // Find the attachment node in all groups
    if (groups.size() > 0)
    {
        attachmentNode = FindAttachmentNodeInGroupsOrFragments(groups, attachmentIndex);
    }

    // Find the attachment node in all fragments, if not found in all groups
    if ((fragments.size() > 0) && (attachmentNode == nullptr))
    {
        attachmentNode = FindAttachmentNodeInGroupsOrFragments(fragments, attachmentIndex);
    }

    return attachmentNode;
}
