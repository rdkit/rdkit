// CommonCS/LibCommon/Src/CDXStdObjects.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 1986-2007, CambridgeSoft Corp., All Rights Reserved

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
#include "CDXDocumentProperty.h"
#include "CDXDocumentPropertyCollection.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <ostream>
#include <stdexcept>
#include "cs_lock.h"
#include "CDXColoredMolecularArea.h"

static cs::CriticalSection s_lock;


// *******************************
// ** class CDXObjectDefinition **
// *******************************
//
// Specialization of CDXObject for CDXObjectDefinition objects

CDXObjectDefinition::CDXObjectDefinition(CDXObjectID id_arg)
 : CDXObject(kCDXObj_ObjectDefinition, id_arg)
{
}

CDXObjectDefinition::CDXObjectDefinition (const CDXObjectDefinition &src)
	:	CDXObject(src)
{
}

CDXObjectDefinition::~CDXObjectDefinition()
{
}

CDXObject*	CDXObjectDefinition::Clone() const
{
	return new CDXObjectDefinition (*this);
}

std::string CDXObjectDefinition::XMLObjectName() const
{
	return kCDXML_objectdefinition;
}

// *******************************
// *** CDXCreateStandardObject ***
// *******************************

static CDXObject *CDXCreateStandardObject(CDXTag tag_arg, CDXObjectID id_arg)
{
	CDXObject *obj;
	// Create an instance of the appropriate subclass for this tag
	switch(tag_arg)
	{
	case kCDXObj_Document:				obj = new CDXDocument(id_arg); 				break;
	case kCDXObj_Page:					obj = new CDXPage(id_arg); 					break;
	case kCDXObj_Group:					obj = new CDXGroup(id_arg); 				break;
	case kCDXObj_Fragment:				obj = new CDXFragment(id_arg);				break;
	case kCDXObj_Node:					obj = new CDXNode(id_arg);					break;
	case kCDXObj_Bond:					obj = new CDXBond(id_arg);					break;
	case kCDXObj_Text:					obj = new CDXText(id_arg);					break;
	case kCDXObj_Graphic:				obj = new CDXGraphic(id_arg);				break;
	case kCDXObj_Curve:					obj = new CDXCurve(id_arg);					break;
	case kCDXObj_EmbeddedObject:		obj = new CDXEmbeddedObject(id_arg);		break;
	case kCDXObj_NamedAlternativeGroup:	obj = new CDXNamedAlternativeGroup(id_arg);	break;
	case kCDXObj_TemplateGrid:			obj = new CDXTemplateGrid(id_arg);			break;
	case kCDXObj_ReactionStep:			obj = new CDXReactionStep(id_arg);			break;
	case kCDXObj_ReactionScheme:		obj = new CDXReactionScheme(id_arg);		break;
	case kCDXObj_ObjectDefinition:		obj = new CDXObjectDefinition(id_arg);		break;
	case kCDXObj_Spectrum:				obj = new CDXSpectrum(id_arg);				break;
	case kCDXObj_ObjectTag:				obj = new CDXObjectTag(id_arg);				break;
	case kCDXObj_Sequence:				obj = new CDXSequence(id_arg);				break;
	case kCDXObj_CrossReference:		obj = new CDXCrossReference(id_arg);		break;
	case kCDXObj_Splitter:				obj = new CDXSplitter(id_arg);				break;
	case kCDXObj_Table:					obj = new CDXTable(id_arg);					break;
	case kCDXObj_BracketedGroup:		obj = new CDXBracketedGroup(id_arg);		break;
	case kCDXObj_BracketAttachment:		obj = new CDXBracketAttachment(id_arg);		break;
	case kCDXObj_CrossingBond:			obj = new CDXCrossingBond(id_arg);			break;
	case kCDXObj_Border:				obj = new CDXBorder(id_arg);				break;
	case kCDXObj_Geometry:				obj = new CDXGeometry(id_arg);				break;
	case kCDXObj_Constraint:			obj = new CDXConstraint(id_arg);			break;
	case kCDXObj_TLCPlate:				obj = new CDXTLCPlate(id_arg);				break;
	case kCDXObj_GEPPlate:				obj = new CDXGEPPlate(id_arg);				break;
	case kCDXObj_TLCLane:				obj = new CDXTLCLane(id_arg);				break;
	case kCDXObj_GEPLane:				obj = new CDXGEPLane(id_arg);				break;
	case kCDXObj_TLCSpot:				obj = new CDXTLCSpot(id_arg);				break;
	case kCDXObj_GEPBand:				obj = new CDXGEPBand(id_arg);				break;
	case kCDXObj_Marker:				obj = new CDXBandMarker(id_arg);			break;
	case kCDXObj_ChemicalProperty:		obj = new CDXChemicalProperty(id_arg);		break;
	case kCDXObj_Arrow:					obj = new CDXArrow(id_arg);					break;
	case kCDXObj_StoichiometryGrid:		obj = new CDXStoichiometryGrid(id_arg);		break;
	case kCDXObj_SGComponent:			obj = new CDXSGComponent(id_arg);			break;
	case kCDXObj_SGDatum:				obj = new CDXSGDatum(id_arg);				break;
	case kCDXObj_BioShape:				obj = new CDXBioShape(id_arg);				break;
	case kCDXObj_PlasmidMap:			obj = new CDXPlasmidMap(id_arg);			break;
	case kCDXObj_PlasmidMarker:			obj = new CDXPlasmidMarker(id_arg);			break;
	case kCDXObj_PlasmidRegion:			obj = new CDXPlasmidRegion(id_arg);			break;
    case kCDXObj_RLogic:				obj = new CDXRLogic(id_arg);				break;
    case kCDXObj_RLogicItem:			obj = new CDXRLogicItem(id_arg);			break;
    case kCDXObj_Annotation:			obj = new CDXAnnotation(id_arg);			break;
    case kCDXObj_DocumentProperties:	obj = new CDXDocumentPropertyCollection(id_arg);	break;
    case kCDXObj_Property:			    obj = new CDXDocumentProperty(id_arg);			    break;
    case kCDXObj_ColoredMolecularArea:  obj = new CDXColoredMolecularArea(id_arg);  break;
	default: throw std::logic_error("Unrecognized tag in CDXCreateStandardObject");
	}
	return obj;
}

void CDXAddStandardFactories(CDXObjectFactory &factory_arg)
{
	PROTECT_GLOBAL_AND_STATIC_DATA(s_lock);
	// We use the same factory function for all the subclasses
	factory_arg.AddFactory(kCDXObj_Document,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Page,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Group,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Fragment,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Node,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Bond,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Text,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Graphic,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Curve,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_EmbeddedObject,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_NamedAlternativeGroup,	CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_TemplateGrid,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_ReactionStep,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_ReactionScheme,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_ObjectDefinition,		CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Spectrum,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_ObjectTag,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Sequence,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_CrossReference,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Splitter,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Table,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_BracketedGroup,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_BracketAttachment,		CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_CrossingBond,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Border,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Geometry,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Constraint,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_TLCPlate,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_GEPPlate,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_TLCLane,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_GEPLane,					CDXCreateStandardObject);	
	factory_arg.AddFactory(kCDXObj_TLCSpot,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_GEPBand,					CDXCreateStandardObject);	
	factory_arg.AddFactory(kCDXObj_Marker,					CDXCreateStandardObject);		
	factory_arg.AddFactory(kCDXObj_ChemicalProperty,		CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_Arrow,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_StoichiometryGrid,		CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_SGComponent,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_SGDatum,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_BioShape,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_PlasmidMap,				CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_PlasmidMarker,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_PlasmidRegion,			CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_RLogic,					CDXCreateStandardObject);
	factory_arg.AddFactory(kCDXObj_RLogicItem,				CDXCreateStandardObject);
    factory_arg.AddFactory(kCDXObj_Annotation,              CDXCreateStandardObject);
    factory_arg.AddFactory(kCDXObj_DocumentProperties,      CDXCreateStandardObject);
    factory_arg.AddFactory(kCDXObj_Property,                CDXCreateStandardObject);
    factory_arg.AddFactory(kCDXObj_ColoredMolecularArea,    CDXCreateStandardObject);
}

std::ostream & operator<<(std::ostream &os, const CDXPropRep &p)
{
	return os << "<" << kCDXML_representElement
		<< " " << kCDXML_attribute << "=\"" << CDXMLAttributeName(p.propTag) << "\""
		<< " " << kCDXML_object << "=\"" << p.objectID << "\""
		<< "/>";
}

void ConnectBonds(CDXObject *obj)
{
	CDXObjectsRange bonds = obj->ContainedObjects(kCDXObj_Bond);

	// Surprisingly, the longer chunk of code is quite a bit faster than the commented-out section.
	// The benefits to calling reserve() are tremendous compared to a bunch of push_back() calls.
	//for (CDXObjectsByTag::const_iterator i = bonds.begin();  i != bonds.end();  ++i)
	//{
	//	CDXBond *b = dynamic_cast<CDXBond *>(GetObject(i));
	//	if (b == 0) throw std::logic_error("Object with bond tag is not a bond");
	//	CDXNode *node1 = dynamic_cast<CDXNode *>(obj->FindByID(b->m_beginNodeID));
	//	CDXNode *node2 = dynamic_cast<CDXNode *>(obj->FindByID(b->m_endNodeID));
	//	if (node1 != NULL && node2 != NULL)
	//		// throw std::runtime_error("IDs at ends of bond do not represent atoms");
	//		b->Connects(node1, node2);
	//}

	if (bonds.begin() != bonds.end())
	{
		// Count the number of bonds per atom
		map<CDXObjectID, int> bondsPerAtom;
		CDXObjectsByTag::const_iterator i;
		for (i = bonds.begin();  i != bonds.end();  ++i)
		{
			CDXBond *b = dynamic_cast<CDXBond *>(GetObject(i));
			if (b == 0) throw std::logic_error("Object with bond tag is not a bond");
			++bondsPerAtom[b->m_beginNodeID];
			++bondsPerAtom[b->m_endNodeID];
		}

		// Reserve storage space for the required number of bonds
		for (map<CDXObjectID, int>::const_iterator bpa = bondsPerAtom.begin(); bpa != bondsPerAtom.end(); ++bpa)
		{
			CDXNode *node = dynamic_cast<CDXNode *>(obj->FindByID(bpa->first));
			if (node != NULL)
			{
				node->m_alphaAtoms.reserve(bpa->second);
				node->m_alphaBonds.reserve(bpa->second);
			}
		}

		// Now store the information on the connecting bonds
		for (i = bonds.begin();  i != bonds.end();  ++i)
		{
			CDXBond *b = dynamic_cast<CDXBond *>(GetObject(i));
			if (b == 0) throw std::logic_error("Object with bond tag is not a bond");
			CDXNode *node1 = dynamic_cast<CDXNode *>(obj->FindByID(b->m_beginNodeID));
			CDXNode *node2 = dynamic_cast<CDXNode *>(obj->FindByID(b->m_endNodeID));
			b->Connects(node1, node2);
		}
	}
}

void CDXPage::FinishReading()
{
    CDXObjectsRange fragments = ContainedObjects(kCDXObj_Fragment);
    for (CDXObjectsByTag::const_iterator i = fragments.begin(); i != fragments.end(); ++i)
    {
        const auto fragment = GetObject(i);

        // Look for bonds that have yet to be completely resolved and resolve them now at the page level
        CDXObjectsRange bonds = fragment->ContainedObjects(kCDXObj_Bond);
        for (CDXObjectsByTag::const_iterator b = bonds.begin(); b != bonds.end(); ++b)
        {
            CDXBond* bond = dynamic_cast<CDXBond*>(b->second);

            if (bond)
            {
                const bool fullyResolved = bond->m_beginNode && bond->m_endNode;
                if (fullyResolved)
                {
                    continue;
                }

                if (!bond->m_beginNode)
                {
                    CDXNode* node = dynamic_cast<CDXNode*>(FindByIDRecursive(bond->m_beginNodeID));
                    if (node != bond->m_beginNode)
                    {
                        bond->m_beginNode = node;
                    }
                }

                if (!bond->m_endNode)
                {
                    CDXNode* node = dynamic_cast<CDXNode*>(FindByIDRecursive(bond->m_endNodeID));
                    if (node != bond->m_endNode)
                    {
                        bond->m_endNode = node;
                    }
                }

                bond->Connects(bond->m_beginNode, bond->m_endNode);
            }
        }
    }
}