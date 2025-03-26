// CommonCS/LibCommon/Src/CDXBond.cpp
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
#include "cs_assert.h"
#include <sstream>
#include <stdexcept>

XMLPutEnum<CDXBondOrder>::Value sXMLBondOrderValues[] = {
	{kCDXBondOrder_Single,		"1"},
	{kCDXBondOrder_Double,		"2"},
	{kCDXBondOrder_Triple,		"3"},
	{kCDXBondOrder_Quadruple,	"4"},
	{kCDXBondOrder_Quintuple,	"5"},
	{kCDXBondOrder_Sextuple,	"6"},
	{kCDXBondOrder_Half,		"0.5"},
	{kCDXBondOrder_OneHalf,		"1.5"},
	{kCDXBondOrder_TwoHalf,		"2.5"},
	{kCDXBondOrder_ThreeHalf,	"3.5"},
	{kCDXBondOrder_FourHalf,	"4.5"},
	{kCDXBondOrder_FiveHalf,	"5.5"},
	{kCDXBondOrder_Dative,		"dative"},
	{kCDXBondOrder_Ionic,		"ionic"},
	{kCDXBondOrder_Hydrogen,	"hydrogen"},
	{kCDXBondOrder_ThreeCenter,	"threecenter"}
};
XMLPutEnum<CDXBondOrder> sXMLBondOrder(sXMLBondOrderValues, sizeof sXMLBondOrderValues, kCDXBondOrder_Single);

void XMLPut(XMLDataSink &sink, CDXBondOrder  v)
{
	if (UINT16(v) == UINT16(-1))
		sink.os << "any";
	else
		sink.os << sXMLBondOrder.asTokens(v);
}

XMLPutEnum<CDXConnectivity>::Value sXMLConnectivityValues[] = {
	{kCDXConnectivity_Unspecified,		"Unspecified"},
	{kCDXConnectivity_Linear,			"Linear"},
	{kCDXConnectivity_Bridged,			"Bridged"},
	{kCDXConnectivity_Staggered,		"Staggered"},
	{kCDXConnectivity_Cyclic,			"Cyclic"}
};

XMLPutEnum<CDXConnectivity> sXMLConnectivity(sXMLConnectivityValues, sizeof sXMLConnectivityValues, kCDXConnectivity_Linear);
void XMLPut(XMLDataSink &sink, CDXConnectivity  v) {	sink.os << sXMLConnectivity.lookup(v); }

XMLPutEnum<CDXBondDisplay>::Value sXMLBondDisplayValues[] = {
	{kCDXBondDisplay_Solid,				"Solid"},
	{kCDXBondDisplay_Dash,				"Dash"},
	{kCDXBondDisplay_Hash,				"Hash"},
	{kCDXBondDisplay_WedgedHashBegin,	"WedgedHashBegin"},
	{kCDXBondDisplay_WedgedHashEnd,		"WedgedHashEnd"},
	{kCDXBondDisplay_Bold,				"Bold"},
	{kCDXBondDisplay_WedgeBegin,		"WedgeBegin"},
	{kCDXBondDisplay_WedgeEnd,			"WedgeEnd"},
	{kCDXBondDisplay_Wavy,				"Wavy"},
	{kCDXBondDisplay_HollowWedgeBegin,	"HollowWedgeBegin"},
	{kCDXBondDisplay_HollowWedgeEnd,	"HollowWedgeEnd"},
	{kCDXBondDisplay_WavyWedgeBegin,	"WavyWedgeBegin"},
	{kCDXBondDisplay_WavyWedgeEnd,		"WavyWedgeEnd"},
	{kCDXBondDisplay_Dot,				"Dot"},
	{kCDXBondDisplay_DashDot,			"DashDot"},
    {kCDXBondDisplay_DottedHydrogen,    "DottedHydrogen"}
};

XMLPutEnum<CDXBondDisplay> sXMLBondDisplay(sXMLBondDisplayValues, sizeof sXMLBondDisplayValues, kCDXBondDisplay_Solid);
void XMLPut(XMLDataSink &sink, CDXBondDisplay  v) {	sink.os << sXMLBondDisplay.lookup(v); }

XMLPutEnum<CDXBondDoublePosition>::Value sXMLBondDoublePositionValues[] = {
	{kCDXBondDoublePosition_UserCenter,	"Center"},
	{kCDXBondDoublePosition_UserRight,	"Right"},
	{kCDXBondDoublePosition_UserLeft,	"Left"},
	{kCDXBondDoublePosition_AutoCenter,	"Center"},
	{kCDXBondDoublePosition_AutoRight,	"Right"},
	{kCDXBondDoublePosition_AutoLeft,	"Left"}
};

XMLPutEnum<CDXBondDoublePosition> sXMLBondDoublePosition(sXMLBondDoublePositionValues, sizeof sXMLBondDoublePositionValues, kCDXBondDoublePosition_AutoCenter);
void XMLPut(XMLDataSink &sink, CDXBondDoublePosition  v) {	sink.os << sXMLBondDoublePosition.lookup(v); }

XMLPutEnum<CDXBondTopology>::Value sXMLBondTopologyValues[] = {
	{kCDXBondTopology_Unspecified,	"Unspecified"},
	{kCDXBondTopology_Ring,			"Ring"},
	{kCDXBondTopology_Chain,		"Chain"},
	{kCDXBondTopology_RingOrChain,	"RingOrChain"}
};

XMLPutEnum<CDXBondTopology> sXMLBondTopology (sXMLBondTopologyValues, sizeof sXMLBondTopologyValues, kCDXBondTopology_Unspecified);
void XMLPut(XMLDataSink &sink, CDXBondTopology  v) {	sink.os << sXMLBondTopology.lookup(v); }

XMLPutEnum<CDXBondReactionParticipation>::Value sXMLBondReactionParticipationValues[] = {
	{kCDXBondReactionParticipation_Unspecified,			"Unspecified"},
	{kCDXBondReactionParticipation_ReactionCenter,		"ReactionCenter"},
	{kCDXBondReactionParticipation_MakeOrBreak,			"MakeOrBreak"},
	{kCDXBondReactionParticipation_ChangeType,			"ChangeType"},
	{kCDXBondReactionParticipation_MakeAndChange,		"MakeAndChange"},
	{kCDXBondReactionParticipation_NotReactionCenter,	"NotReactionCenter"},
	{kCDXBondReactionParticipation_NoChange,			"NoChange"},
	{kCDXBondReactionParticipation_Unmapped,			"Unmapped"}
};

XMLPutEnum<CDXBondReactionParticipation> sXMLBondReactionParticipation(sXMLBondReactionParticipationValues, sizeof sXMLBondReactionParticipationValues, kCDXBondReactionParticipation_Unspecified);
void XMLPut(XMLDataSink &sink, CDXBondReactionParticipation  v) {	sink.os << sXMLBondReactionParticipation.lookup(v); }

XMLPutEnum<CDXBondCIPType>::Value sXMLBondCIPTypeValues[] = {
	{kCDXCIPBond_Undetermined,	"U"},
	{kCDXCIPBond_None,			"N"},
	{kCDXCIPBond_E,				"E"},
	{kCDXCIPBond_Z,				"Z"}
};

XMLPutEnum<CDXBondCIPType> sXMLBondCIPType(sXMLBondCIPTypeValues, sizeof sXMLBondCIPTypeValues, kCDXCIPBond_Undetermined);
void XMLPut(XMLDataSink &sink, CDXBondCIPType  v) {	sink.os << sXMLBondCIPType.lookup(v); }

// ***********************
// ** class CDXBond **
// ***********************
//
// Specialization of CDXObject for CDXBond objects

CDXBond::CDXBond(CDXObjectID id_arg)
	: CDXGraphicObject	(kCDXObj_Bond, id_arg)
	, m_beginNodeID		(0)
	, m_endNodeID		(0)
	, m_bondOrder		(kCDXBondOrder_Single)
	, m_display			(kCDXBondDisplay_Solid)
	, m_display2		(kCDXBondDisplay_Solid)
	, m_doublePosition	(kCDXBondDoublePosition_AutoCenter)
	, m_topology		(kCDXBondTopology_Unspecified)
	, m_rxnParticipation(kCDXBondReactionParticipation_Unspecified)
	, m_CIP				(kCDXCIPBond_Undetermined)
	, m_beginNode		(NULL)
	, m_endNode			(NULL)
	, m_bondOrdering	(NULL)
    , m_crossingBonds	(NULL)
    , m_beginAttach		(-1)
    , m_endAttach		(-1)
    , m_connectivity	(kCDXConnectivity_Unspecified)
    , m_beginExternalNum(0)
    , m_endExternalNum	(0)
{
}

CDXBond::CDXBond (const CDXBond &src)
	:	CDXGraphicObject	(src)
	,	m_beginNodeID		(src.m_beginNodeID)
	,	m_endNodeID			(src.m_endNodeID)
	,	m_bondOrder			(src.m_bondOrder)
	,	m_display			(src.m_display)
	,	m_display2			(src.m_display2)
	,	m_doublePosition	(src.m_doublePosition)
	,	m_topology			(src.m_topology)
	,	m_rxnParticipation	(src.m_rxnParticipation)
	,	m_CIP				(src.m_CIP)
	,	m_beginNode			(NULL)
	,	m_endNode			(NULL)
	,	m_bondOrdering		((src.m_bondOrdering == NULL) ? NULL : new vector<CDXObjectID>(*src.m_bondOrdering))
	,	m_crossingBonds		((src.m_crossingBonds == NULL) ? NULL : new vector<CDXObjectID>(*src.m_crossingBonds))
	,	m_beginAttach		(src.m_beginAttach)
	,	m_endAttach			(src.m_endAttach)
	,   m_connectivity		(src.m_connectivity)
    ,   m_beginExternalNum  (src.m_beginExternalNum)
    ,   m_endExternalNum    (src.m_endExternalNum)
    ,   isRouted            (src.isRouted)
{
}

CDXBond::~CDXBond()
{
	DeleteAndNull(m_bondOrdering);
	DeleteAndNull(m_crossingBonds);
}

CDXObject*	CDXBond::Clone() const
{
	return new CDXBond (*this);
}

std::string CDXBond::XMLObjectName() const
{
	return kCDXML_bond;
}

void CDXBond::FinishReading()
{
	if (m_beginNodeID == 0 || m_endNodeID == 0)
		throw std::runtime_error("Bonds must have a Begin and an End");
//	12/01 (7.0.1): it can happen that nodes are equal, just ignore .. jdd
//	if (m_beginNodeID == m_endNodeID)
//		throw std::runtime_error("Bonds must have a distinct Begin and End");
}

void CDXBond::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bond_Begin:
		m_beginNodeID = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_End:
		m_endNodeID = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_Order:
		m_bondOrder = (CDXBondOrder) src_arg.GetUINT(size_arg);
		if (size_arg == 1 && m_bondOrder == 0xFF)
			m_bondOrder = (CDXBondOrder) -1;
		else if (size_arg == 2 && m_bondOrder == 0xFFFF)
			m_bondOrder = (CDXBondOrder) -1;
		assert(sizeof(CDXBondOrder) <= 4);
		break;
	case kCDXProp_Bond_Connectivity:
		m_connectivity = (CDXConnectivity) src_arg.GetUINT(size_arg);
        break;
    case kCDXProp_Bond_Connectivity_Routed:
        SetIsRouted((size_arg != 0) && src_arg.GetINT(size_arg) != 0);
        break;
    case kCDXProp_Bond_BeginExternalNum:
        m_beginExternalNum = src_arg.GetINT(size_arg);
        break;
    case kCDXProp_Bond_EndExternalNum:
        m_endExternalNum = src_arg.GetINT(size_arg);
        break;
	case kCDXProp_Bond_Display:
		m_display = (CDXBondDisplay) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_Display2:
		m_display2 = (CDXBondDisplay) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_DoublePosition:
		m_doublePosition = (CDXBondDoublePosition) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_RestrictTopology:
		m_topology = (CDXBondTopology) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_RestrictRxnParticipation:
		m_rxnParticipation = (CDXBondReactionParticipation) src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_CIPStereochemistry:
		m_CIP = (CDXBondCIPType)src_arg.GetINT(size_arg);
		break;
	case kCDXProp_Bond_BeginAttach:
		m_beginAttach = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_EndAttach:
		m_endAttach = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Bond_BondOrdering:
	{
		DeleteAndNull(m_bondOrdering);
		m_bondOrdering = new vector<CDXObjectID>;
		CDXReadItems(*m_bondOrdering, size_arg/sizeof(UINT32), src_arg);
		break;
	}
	case kCDXProp_Bond_CrossingBonds:
	{
		DeleteAndNull(m_crossingBonds);
		m_crossingBonds = new vector<CDXObjectID>;
		CDXReadItems(*m_crossingBonds, size_arg/sizeof(UINT32), src_arg);
		break;
	}
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}


static void ReadSomeIDs(const std::string &value_arg, std::vector<CDXObjectID> &dest)
{
	// This verbose function could be written a lot more simply as follows.
	// The current longer implementation is intended simply as a (signficant) performance optimization
	//	std::istringstream is(value_arg);
	//	m_crossingBonds.clear();
	//	std::copy(std::istream_iterator<CDXObjectID>(is),
	//			  std::istream_iterator<CDXObjectID>(),
	//			  std::back_inserter(m_crossingBonds));

	if (value_arg.empty())
		return;

	dest.clear();

	int numIDs = 0;
	bool prevWasDigit = false;
	string::const_iterator c;
	for (c = value_arg.begin(); c != value_arg.end(); ++c)
	{
		if (isdigit((unsigned char)*c))
		{
			if (!prevWasDigit)
			{
				++numIDs;
				prevWasDigit = true;
			}
		}
		else if (prevWasDigit)
		{
			prevWasDigit = false;
		}
	}

	dest.resize(numIDs);

	prevWasDigit = false;
	int idNum = -1;
	for (c = value_arg.begin(); c != value_arg.end(); ++c)
	{
		if (isdigit((unsigned char)*c))
		{
			if (!prevWasDigit)
			{
				dest[++idNum] = (*c - '0');
				prevWasDigit = true;
			}
			else
			{
				dest[idNum] = dest[idNum] * 10 + (*c - '0');
			}
		}
		else if (prevWasDigit)
		{
			prevWasDigit = false;
		}
	}
}


void CDXBond::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Bond_Begin:
		m_beginNodeID = atoi(value_arg.c_str());
		break;
	case kCDXProp_Bond_End:
		m_endNodeID = atoi(value_arg.c_str());
		break;
	case kCDXProp_Bond_Order:
		if (value_arg == "any")
			m_bondOrder = kCDXBondOrder_Any;
		else
			m_bondOrder = sXMLBondOrder.lookupTokens(value_arg);
		break;
	case kCDXProp_Bond_Connectivity:
		m_connectivity = sXMLConnectivity.lookup(value_arg);
		break;
    case kCDXProp_Bond_Connectivity_Routed:
        SetIsRouted(value_arg == "yes");
        break;
    case kCDXProp_Bond_BeginExternalNum:
        m_beginExternalNum = atoi(value_arg.c_str());
        break;
    case kCDXProp_Bond_EndExternalNum:
        m_endExternalNum = atoi(value_arg.c_str());
        break;
	case kCDXProp_Bond_Display:
		m_display = sXMLBondDisplay.lookup(value_arg);
		break;
	case kCDXProp_Bond_Display2:
		m_display2 = sXMLBondDisplay.lookup(value_arg);
		break;
	case kCDXProp_Bond_DoublePosition:
		m_doublePosition = sXMLBondDoublePosition.lookup(value_arg);
		break;
	case kCDXProp_Bond_RestrictTopology:
		m_topology = sXMLBondTopology.lookup(value_arg);
		break;
	case kCDXProp_Bond_RestrictRxnParticipation:
		m_rxnParticipation = sXMLBondReactionParticipation.lookup(value_arg);
		break;
	case kCDXProp_Bond_CIPStereochemistry:
		m_CIP = sXMLBondCIPType.lookup(value_arg);
		break;
	case kCDXProp_Bond_BeginAttach:
		m_beginAttach = atoi(value_arg.c_str());
		break;
	case kCDXProp_Bond_EndAttach:
		m_endAttach = atoi(value_arg.c_str());
		break;
	case kCDXProp_Bond_BondOrdering:
	{
		DeleteAndNull(m_bondOrdering);
		m_bondOrdering = new vector<CDXObjectID>;
		ReadSomeIDs(value_arg, *m_bondOrdering);
		break;
	}
	case kCDXProp_Bond_CrossingBonds:
	{
		DeleteAndNull(m_crossingBonds);
		m_crossingBonds = new vector<CDXObjectID>;
		ReadSomeIDs(value_arg, *m_crossingBonds);
		break;
	}
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, value_arg);
		break;
	}
}


void CDXBond::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_Bond_Begin, m_beginNodeID );
	sink_arg.PutAttribute( kCDXProp_Bond_End,   m_endNodeID );

	if (m_bondOrder != kCDXBondOrder_Single)
		sink_arg.PutAttribute( kCDXProp_Bond_Order, INT16(m_bondOrder) );
	if (m_connectivity != kCDXConnectivity_Unspecified)
    {
        if (isRouted)
        {
            sink_arg.PutAttribute( kCDXProp_Bond_Connectivity_Routed, UINT8(isRouted));
        }
        
        sink_arg.PutAttribute( kCDXProp_Bond_Connectivity, UINT16(m_connectivity) );
    }
    
    if (m_beginExternalNum > 0)
    {
        sink_arg.PutAttribute(kCDXProp_Bond_BeginExternalNum, INT8(m_beginExternalNum));
    }
    
    if (m_endExternalNum > 0)
    {
        sink_arg.PutAttribute(kCDXProp_Bond_EndExternalNum, INT8(m_endExternalNum));
    }
    
	if (m_display != kCDXBondDisplay_Solid)
		sink_arg.PutAttribute( kCDXProp_Bond_Display, UINT16(m_display) );
	if (m_display2 != kCDXBondDisplay_Solid)
		sink_arg.PutAttribute( kCDXProp_Bond_Display2, UINT16(m_display2) );
	if (m_doublePosition != kCDXBondDoublePosition_AutoCenter)
		sink_arg.PutAttribute( kCDXProp_Bond_DoublePosition, UINT16(m_doublePosition) );
	if (m_topology != kCDXBondTopology_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Bond_RestrictTopology, UINT8(m_topology) );			// verify size -heh 5/9/99
	if (m_rxnParticipation != kCDXBondReactionParticipation_Unspecified)
		sink_arg.PutAttribute( kCDXProp_Bond_RestrictRxnParticipation, UINT8 (m_rxnParticipation) );	// verify size - heh 5/9/99
	if (m_CIP != kCDXCIPBond_Undetermined)
		sink_arg.PutAttribute( kCDXProp_Bond_CIPStereochemistry, INT8 (m_CIP) );	// verify size - heh 5/9/99
	if (m_beginAttach >= 0)
		sink_arg.PutAttribute( kCDXProp_Bond_BeginAttach, INT8 (m_beginAttach) );
	if (m_endAttach >= 0)
		sink_arg.PutAttribute( kCDXProp_Bond_EndAttach, INT8 (m_endAttach) );
	if (m_bondOrdering != NULL && !m_bondOrdering->empty())
	{
		sink_arg.PutTag( kCDXProp_Bond_BondOrdering );
		sink_arg.Put( UINT16( 4 * m_bondOrdering->size() ) );
		Put( sink_arg, *m_bondOrdering );
	}
	if (m_crossingBonds != NULL && !m_crossingBonds->empty())
	{
		sink_arg.PutTag( kCDXProp_Bond_CrossingBonds );
		sink_arg.Put( UINT16( 4 * m_crossingBonds->size() ) );
		Put( sink_arg, *m_crossingBonds );
	}
}

void CDXBond::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Begin, m_beginNodeID );
	CDXMLPutAttribute(sink_arg, kCDXProp_Bond_End,   m_endNodeID );

	if (m_bondOrder != kCDXBondOrder_Single)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Order, m_bondOrder );
	if (m_connectivity != kCDXConnectivity_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Connectivity, m_connectivity );
    
    if (IsRouted())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Connectivity_Routed, isRouted );
    }
    
    if (m_beginExternalNum > 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bond_BeginExternalNum, m_beginExternalNum);
    }
    
    if (m_endExternalNum > 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Bond_EndExternalNum, m_endExternalNum);
    }
    
    if (m_display != kCDXBondDisplay_Solid)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Display, m_display );
	if (m_display2 != kCDXBondDisplay_Solid)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_Display2, m_display2 );
	if (m_doublePosition == kCDXBondDoublePosition_UserCenter ||
		m_doublePosition == kCDXBondDoublePosition_UserRight ||
		m_doublePosition == kCDXBondDoublePosition_UserLeft)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_DoublePosition, m_doublePosition );
	if (m_topology != kCDXBondTopology_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_RestrictTopology, m_topology );
	if (m_rxnParticipation != kCDXBondReactionParticipation_Unspecified)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_RestrictRxnParticipation, m_rxnParticipation );
	if (m_CIP != kCDXCIPBond_Undetermined)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_CIPStereochemistry, m_CIP );
	if (m_bondOrdering != NULL && !m_bondOrdering->empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_BondOrdering, *m_bondOrdering);
	if (m_crossingBonds != NULL && !m_crossingBonds->empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_CrossingBonds, *m_crossingBonds);
	if (m_beginAttach >= 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_BeginAttach, m_beginAttach );
	if (m_endAttach >= 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Bond_EndAttach, m_endAttach );
}

void CDXBond::Connects(CDXNode* node1, CDXNode* node2)
{
    if (!node1 || !node2)
    {
        // Incomplete information so we defer until the CDXPage is finished
        return;
    }

    m_beginNode = node1;
    m_endNode = node2;

    m_beginNodeID = node1->GetObjectID();
    m_endNodeID = node2->GetObjectID();

    node1->m_alphaBonds.push_back(this);
    node2->m_alphaBonds.push_back(this);

    node1->m_alphaAtoms.push_back(node2);
    node2->m_alphaAtoms.push_back(node1);
}

/**
 * Does this bond require that both atoms reside share the same parent?
 */
bool CDXBond::MustBothAtomsResideInSameFragment() const
{
    return (m_display2 != kCDXBondDisplay_DottedHydrogen);
}

CDXNode* CDXBond::OtherNode(CDXNode* thisNode)
{
	if (m_beginNode == thisNode)
	{
		return m_endNode;
	}

	if (m_endNode == thisNode)
	{
		return m_beginNode;
	}

	return NULL;
}
