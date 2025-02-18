/*
File:		CDXUtils.cpp
Purpose:	Utility functions for manipulating CDX.
Copyright:	(c) 2002-2009, CambridgeSoft Corp., All Rights Reserved

--------------------------------------------------------------------------------------------+
HEH  04/07/09	BreakDiscontinuousSetsIntoFragments(): Handle case where parent is a group.	|
HEH  02/27/07	Restore CDXBuildObjectMap().												|
JSB  01/16/07	Remove CDXObjectAddress, CDXFindHighestUsedID, CDXCompileIdsInUse,			|
				CDXBuildObjectMap, CDXTranslate, CDXNumPagesInDocument, CDXGetPage,			|
				CDXFindFirstEnclosingObjectOfType, CDXMergeDocuemnt, CDDXAddPage,			|
				CDXAddFragment, CDXAddObject, CDXAddText CDXFormObjectsAddress,				|
				CDXFindObjectByAddress, CDXObjectsAddressString, AnyArrows,					|
				GetNodesTextNode, SetNodesText.												|
HEH  07/22/04	CSBR-46228: CDXIsValid(): Don't flag dupl ID's as error if not atom or bond.|
HEH  07/14/04	Added CDXFindFirstEnclosingObjectOfType() (two overloads).					|
HEH  04/12/04	CDXReadOrComputeBoundingBox() gets a fIgnoreObjectTags argument.			|
HEH  03/14/04	Move CdxFindMedianBondLength() from CdxToCcMolecule.cpp to CdxUtils.cpp,	|
				and consolidate with preexisting CdxGetMedianBondLength().					|
HEH  03/09/04	Added CDXGetMedianBondLength().												|
HEH  03/08/03	Replace SizeOfCDXObjectsRange() by CDXObjectsRange::size().					|
HEH  03/06/03	1. CSBR-33644: CDXGetTextExtent(): Handle case of font size of "-1".		|
				2. ReadOrCompute_Recurse(): Must pursue text children of nodes w/positions.	|
HEH  02/28/03	CDXFindObjectByAddress(): Search entire tree instead of giving up.			|
HEH  02/27/03	Added CDXObjectAddressString().												|
HEH  01/22/03	CDXIsValid(): Permit bond's nodes to be null.								|
HEH  12/10/02	Moved GetNodesTextNode() from CdxToCcMolecule.cpp to CdxUtils.cpp.			|
				Added SetNodesText().														|
HEH  10/25/02	Moved CDXBuildObjectMap() from CdxMol.cpp to CdxUtils.cpp.  Added			|
				"fIgnoreIdZero" arg. to it and to CDXCompileIdsInUse().						|
HEH  12/23/01	Added CDXCompileIdsInUse(), CDXReadOrComputeBoundingBox(), CDXTranslate(),	|
				CDXNumPagesInDocument(), CDXGetPage(), CDXMergeDocument(), CDXAddPage(),	|
				CDXAddObject(), CDXAddText(), CDXGetTextExtent(), CDXFormObjectsAddress(),	|
				CDXFindObjectByAddress().													|
HEH  12/23/01	Created by moving CDXIsValid() here from CDXDocument.cpp.					|
--------------------------------------------------------------------------------------------+
*/

#include "cs_assert.h"
#include "cs_charUtils.h"
#include "cs_math.h"

#include "CDXUtils.h"
#include "CDMap.h"
#include "UTF8String.h"

#include <iostream>	// temp for clog

namespace
{
    /**
     * Check whether the attached data property is present in the fragment
     *
     * @param fragment The fragment
     *
     * @return true if the attached-date property is present, otherwise false
     */
    bool ValidateFragmentHasAttachedDataProperty(CDXFragment* fragment)
    {
        CDXPage* cdxPage = dynamic_cast<CDXPage*>(fragment->GetParent());
        ASSERT (cdxPage != NULL);
        
        CDXObjectsRange pageObjects = cdxPage->ContainedObjects(kCDXObj_ChemicalProperty);
        for (auto object : pageObjects)
        {
            CDXChemicalProperty *property = dynamic_cast<CDXChemicalProperty*>(object.second);
            if (property && (std::find(property->m_propertyTypes.begin(), property->m_propertyTypes.end(), kCDXChemicalPropertyTypeAttachedData) !=  property->m_propertyTypes.end()))
            {
                CDXObjectsRange titles = fragment->ContainedObjects(kCDXObj_Text);
                for (auto title : titles)
                {
                    // Check the display ID is associated with the property
                    CDXText* text = dynamic_cast<CDXText*>(title.second);
                    if (text && (text->GetObjectID() == property->m_displayID))
                    {
                        return true;
                    }
                }
            }
        }
        
        return false;
    }
}

// *************************************************************************************
// ******* CDXIsValid   Quick, non-exhaustive integrity check of an object tree ********
// *************************************************************************************
//
bool CDXIsValid (const CDXObject* pObj, std::set<const CDXObject*>* pAlreadyChecked)	// typically pObj is a CDXDocument*
{
	std::set<const CDXObject*>	alreadyChecked_local;
	if (!pObj)
		goto Corrupt;
	if (!pAlreadyChecked)
	{
		pAlreadyChecked = &alreadyChecked_local;

#ifdef _WINDOWS
		ASSERT (_CrtCheckMemory());
#endif

	}
	if (pAlreadyChecked->find (pObj) != pAlreadyChecked->end())
		return true;
	if (!pObj->IsValid())
		return false;

	try
	{
#ifdef _AFX
		if (!AfxIsValidAddress (pObj, sizeof (CDXObject), true))	// T=check write access
			goto Corrupt;
#endif
		const CDXObjectID	id = pObj->GetObjectID();
		if ((int)id < 0)
			goto Corrupt;

		// Recurse on children
		CDXObjectsRange children = pObj->ContainedObjects();
		set<CDXObjectID>	siblings;
		for (CDXObjectsByTag::const_iterator itChildren = children.begin();  itChildren != children.end();  itChildren++)
		{
			const CDXTag		tag = itChildren->first;
			const CDXObject*	pChild = itChildren->second;
			const CDXObjectID	childId = pChild->GetObjectID();
#ifdef _AFX
		if (!AfxIsValidAddress (pChild, sizeof (CDXObject), true))	// T=check write access
			goto Corrupt;
#endif
			if (pChild->GetParent() != pObj)
				goto Corrupt;
			if (childId  &&  siblings.find (childId) != siblings.end())
			{
				if (tag == kCDXObj_Node  ||  tag == kCDXObj_Bond)
					goto Corrupt;
#if _HEH_
				ASSERT (false);
#endif
			}
			siblings.insert (childId);
			if (!CDXIsValid (pChild, pAlreadyChecked))
				goto Corrupt;
			// Check bond's nodes are defined
			const CDXBond*	pBond = dynamic_cast<const CDXBond*>(pChild);
			if (pBond)
			{
				if (pBond->m_beginNode  &&  !CDXIsValid (pBond->m_beginNode, pAlreadyChecked))
					goto Corrupt;
				if (pBond->m_endNode    &&  !CDXIsValid (pBond->m_endNode, pAlreadyChecked))
					goto Corrupt;
			}
		}
		pAlreadyChecked->insert (pObj);
	}
	catch (...)
	{ goto Corrupt; }

	return true;

Corrupt:
	return false;
}

/*
+===========================================================================================+
| CDXReadOrComputeBoundingBox	Returns the bounding box occupied by the given object,		|
|								including any contained objects.  If an object has a		|
|								defined bounding box, that can be used directly.			|
|								Otherwise, the box is manually computed.					|
|								Is like CDXObject::BoundsIncludingChildren(), but attempts	|
|								to calculate bounding boxes when absent.					|
|								Penned by HEH, 12/19/01.									|
+-------------------------------------------------------------------------------------------+
|[R-]pObj				The top-level object, for which the bounding box is being computed.	|
|[R-]fIgnoreCachedInfo	If False, an object's stored bounding box information is used when	|
|						present and manual calculation is avoided.  If True, such cached	|
|						information is ignored, and the calculation is exhaustive.			|
|[R-]fIgnoreObjectTags	If True, object tags (such as atom numbers and stereo indicators)	|
|						nested under a node will not be considered.							|
+===========================================================================================+
*/
typedef pair<CDXRectangle,bool>	RectPair;
static void  ReadOrCompute_Recurse (RectPair& rp, const CDXObject* pObj, bool fIgnoreCachedInfo, bool fIgnoreObjectTags)
{
	if (!pObj)
		return;
	RectPair	rp2 = make_pair (CDXRectangle (0,0,0,0), false);
	// General graphic object
	const CDXGraphicObject* pGo = dynamic_cast<const CDXGraphicObject*>(pObj);
	if (pGo)
	{
		if (!fIgnoreCachedInfo  &&  pGo->Known (CDXGraphicObject::has_boundingBox))
		{
			rp2.second = true;
			rp2.first = pGo->BoundingBox();
		}
		else if (pGo->Known (CDXGraphicObject::has_2dPosition))
		{
			rp2.second = true;
			CDXRectangle&	rect2 = rp2.first;
			CDXPoint2D		pos = pGo->Position();
			rect2.left = rect2.right  = pos.x;
			rect2.top  = rect2.bottom = pos.y;
			// If a string, figure its dimension.
			const CDXText*		pText = dynamic_cast<const CDXText*>(pGo);
			if (pText)
			{
				const CDXPoint2D	extent = CDXGetTextExtent (pText->GetText());
				rect2.right  += extent.x;
				rect2.bottom += extent.y;
			}
		}
	}
#if 0
	// A page's bounding box refers (?) to its potential size, not its contained
	// objects.
	// Page
	else if (!fIgnoreCachedInfo)
	{
		const CDXPage*	pPage = dynamic_cast<const CDXPage*>(pObj);
		if (pPage)
		{
			if (pPage->Known (CDXPage::has_boundingBox))
			{
				rp2.second = true;
				rp2.first = pPage->BoundingBox();
			}
		}
	}
#endif

	// Incorporate new rectangle into what we've already got.
	if (rp2.second)
	{
		if (!rp.second)
			rp.first = rp2.first;
		else
			rp.first |= rp2.first;
		rp.second = true;
//		if (!fIgnoreCachedInfo)
//			return;	// No need to explore children
	}
	
	// Explore children -- but exclude contracted labels, etc
	CDXObjectsRange	range = pObj->ContainedObjects();
	for (CDXObjectsByTag::const_iterator it = range.begin();  it != range.end();  it++)
	{
		const CDXObject*		pChild = it->second;
		// Ignore tags
		if (fIgnoreObjectTags  &&  pChild->GetTag() == kCDXObj_ObjectTag)
			continue;
		// Exclude contracted labels, etc
		if (dynamic_cast<const CDXNode *>(pObj) != NULL && pChild->GetTag() != kCDXObj_ObjectTag)
			continue;
		ReadOrCompute_Recurse (rp, pChild, fIgnoreCachedInfo, fIgnoreObjectTags);
	}
} // ReadOrCompute_Recurse()
//------------------
CDXRectangle CDXReadOrComputeBoundingBox (const CDXObject* pObj, bool fIgnoreCachedInfo /* = false */, bool fIgnoreObjectTags /* = false */)
{
	RectPair	rp = make_pair (CDXRectangle (0,0,0,0), false);
	ReadOrCompute_Recurse (rp, pObj, fIgnoreCachedInfo, fIgnoreObjectTags);
	return rp.first;	// will be (0,0,0,0) if no bounding boxes exist.
}


/*
+===========================================================================================+
| CDXBuildObjectMap		Scans a CDX object tree, building a map relating every object ID	|
|						to its CDXObject.  If object ID's duplicate, this map is not very	|
|						useful.  Penned by HEH, 8/21/99.									|
|						Cf. CDXCompileIdsInUse().											|
+-------------------------------------------------------------------------------------------+
| top					The object at the top of the tree to be scanned; is typically the	|
|						CDXDocument* itself.												|
| objMap				Returned as the answer.												|
| pDuplicateIdsPresent	If provided, is returned False iff no ID is used more than once,	|
|						ignoring the special ID value zero.									|
| fIgnoreIdZero			If True, "pDulicateIdsPresent" is not affected by duplicate ID's	|
|						of value zero.														|
+===========================================================================================+
*/
static void CDXBuildObjectMap_Recurse	( const CDXObject*						pObj
										, map<CDXObjectID,const CDXObject*>&	objMap
										, bool*									dup
										, bool									fIgnoreIdZero
										)
{
	CDXObjectsRange	children = pObj->ContainedObjects();
	for (CDXObjectsByTag::const_iterator  itChildren = children.begin();  itChildren != children.end();  itChildren++)
	{
		const CDXObject*	pChildObj = (itChildren->second);
		CDXObjectID			id = pChildObj->GetObjectID();
		if (objMap.find (id) != objMap.end())
		{
			if (id  ||  !fIgnoreIdZero)
				*dup = true;
		}
		else
			objMap [id] = pChildObj;
		CDXBuildObjectMap_Recurse (pChildObj, objMap, dup, fIgnoreIdZero);
	}
}
// - - - - - - - - - - - - - - - - - - -
void CDXBuildObjectMap	( const CDXObject&						top
						, map<CDXObjectID,const CDXObject*>&	objMap
						, bool *								pDuplicateIdsPresent
						, bool									fIgnoreIdZero	// true
						)
{
	bool	dups_buff;
	if (!pDuplicateIdsPresent)
		pDuplicateIdsPresent = &dups_buff;
	*pDuplicateIdsPresent = false;
	objMap.clear();
	objMap [top.GetObjectID()] = &top;
	CDXBuildObjectMap_Recurse (&top, objMap, pDuplicateIdsPresent, fIgnoreIdZero);
}



/*
+===========================================================================================+
| CDXGetTextExtent		Returns the approximate dimensions of the string, in				|
|						CDXCoordinate's.  Penned by HEH, 12/20/01.							|
+-------------------------------------------------------------------------------------------+
|[R-]str				The string of interest.												|
|<FRV>					The X and Y extent of the string, in CDXCoordinates, bundled into	|
|						a CDXPoint2D.														|
+===========================================================================================+
*/
CDXPoint2D CDXGetTextExtent (const CDXString& str)
{
	const double	kCharWidth = 0.6;	// width is 0.6 times height, on average.  Crude!
	CDXPoint2D	extent (0,0);
	int	lastStart = 0,
		lastSize  = 10 * kNumUnitsPerPoint;	// Whatever the default font size is (style units, i.e. 1/20 point)
	for (int i = 0;  i < str.nstyles();  i++)
	{
		const CDXStyle&	style = str.style (i);
//		int	family = style.family;
//		CDXFontTableEntry	font = pDoc->GetFontTable()...;
		const int	len = style.startChar - lastStart;
		if (len)
		{
			extent.x += cs::Round (len * lastSize * kCharWidth);
			extent.y = max (extent.y, (CDXCoordinate)lastSize);
		}
		lastStart = style.startChar;
		if (style.size  &&  style.size != (UINT16)-1)
			lastSize = style.size;
	}

	// Get the last segment.
	const ptrdiff_t	len = str.length() - lastStart;
	if (len)
	{
		extent.x += cs::Round (len * lastSize * kCharWidth);
		extent.y = max (extent.y, (CDXCoordinate)lastSize);
	}

	// Convert from 1/20 points to CDXCoordinate's
	extent.x *= 65535;
	extent.y *= 65535;
	return extent;
} // CDXGetTextExtent()



/*
+===========================================================================================+
| CdxFindMedianBondLength	Return the median bond length of those present at & under the	|
|							given CDXObject.  Penned by HEH, 12/18/01.  Mean added 3/8/04.	|
+===========================================================================================+
*/
double CdxFindMedianBondLength (const CDXObject	&top,
								Dimensions2or3	dimension /* = kIn2D */,
								double			*pMean /* = NULL */)
{
	ASSERT (CDXIsValid (&top));
	// Skip over groups using a tree-walking helper class.
	CdxGetChildrenIgnoringGroups	noGroups (&top, kCDXObj_Bond);

	vector<double>	bdLens;
	double			sum = 0;
	const CDXObject	*pObj;
	while (noGroups.NextChild (pObj))
	{
		const CDXBond*	pBond = FAST_dynamic_cast<const CDXBond*>(pObj);
		if (pBond == NULL)
			continue;
		const CDXNode	*pBeginNode = pBond->m_beginNode,
						*pEndNode   = pBond->m_endNode;
		bool	fIn3D = dimension == kIn3D  &&  pBeginNode->KnownPosition3D()  &&  pEndNode->KnownPosition3D();
		double	bdLen_sq;
		if (fIn3D)
		{
			const CDXPoint3D	pos_1 = pBeginNode->Position3D(),
								pos_2 = pEndNode  ->Position3D();
			bdLen_sq = double (pos_1.x - pos_2.x) * (pos_1.x - pos_2.x) +
					   double (pos_1.x - pos_2.y) * (pos_1.x - pos_2.y) +
					   double (pos_1.x - pos_2.z) * (pos_1.x - pos_2.z);
		}
		else if (pBeginNode->KnownPosition()  &&  pEndNode->KnownPosition())
		{
			const CDXPoint2D	pos_1 = pBeginNode->Position(),
								pos_2 = pEndNode  ->Position();
			bdLen_sq = double (pos_1.x - pos_2.x) * (pos_1.x - pos_2.x) +
					   double (pos_1.y - pos_2.y) * (pos_1.y - pos_2.y);
		}
		else
			continue;
		bdLens.push_back (bdLen_sq);
		if (pMean)
			sum += sqrt (bdLen_sq);
	}
	if (bdLens.empty())
		return 0.;
	if (pMean)
		*pMean = sum / bdLens.size();
	nth_element (bdLens.begin(), bdLens.begin() + (bdLens.size() / 2), bdLens.end());
	const double	median = sqrt (bdLens[bdLens.size() / 2]);
	return median;
}

/*
+===================================================================================================+
| AnyChargeObjects			Penned by JSB.	Is there a charge object anywhere in the object tree?	|
+===================================================================================================+
*/
bool AnyChargeObjects(const CDXObject *cdxObj)
{
	ASSERT(cdxObj != NULL);
	const CDXGraphic *cdxGraphic = dynamic_cast<const CDXGraphic *>(cdxObj);
	const CDXText *cdxText = dynamic_cast<const CDXText *>(cdxObj);
	if (cdxGraphic != NULL && cdxGraphic->m_graphicType == kCDXGraphicType_Symbol)
	{
		if (cdxGraphic->m_symbolType == kCDXSymbolType_LonePair			|| cdxGraphic->m_symbolType == kCDXSymbolType_Electron
			|| cdxGraphic->m_symbolType == kCDXSymbolType_RadicalCation	|| cdxGraphic->m_symbolType == kCDXSymbolType_RadicalAnion
			|| cdxGraphic->m_symbolType == kCDXSymbolType_CirclePlus	|| cdxGraphic->m_symbolType == kCDXSymbolType_CircleMinus
			|| cdxGraphic->m_symbolType == kCDXSymbolType_Plus			|| cdxGraphic->m_symbolType == kCDXSymbolType_Minus
            || cdxGraphic->m_symbolType == kCDXSymbolType_LonePairBar)
		{
			return true;
		}
	}
	else if (cdxText != NULL)
	{
        const string &str = cdxText->GetText().str();
        if ((str.length() == 1) && cs::IsChargeSign(str[0])
            && !(cdxText->Known(cdxText->has_interpretChemically) && !cdxText->m_interpretChemically))
        {
            return true;
        }
	}

	if (dynamic_cast<const CDXNode *>(cdxObj) == NULL)
	{
		// For each object in this one (but skip atoms because we don't want to burrow into contracted labels)...
		CDXObjectsRange containedObjects = cdxObj->ContainedObjects();
		for (CDXObjectsByTag::const_iterator objIter = containedObjects.begin();  objIter != containedObjects.end();  ++objIter)
			if (AnyChargeObjects(objIter->second))
				return true;
	}

	return false;
}

/*
+===================================================================================================+
| AnyAbsRacRelMarkers		Penned by JSB.															|
+===================================================================================================+
*/
bool AnyAbsRacRelMarkers(const CDXObject *cdxObj, bool &hasAbsMarker, bool &hasRacMarker, bool &hasRelMarker)
{
	hasAbsMarker = false;
	hasRacMarker = false;
	hasRelMarker = false;
	CdxGetChildrenIgnoringGroups	ArrowsWoGroups (cdxObj, kCDXObj_Graphic);
	const CDXObject *pObj = NULL;
	while (ArrowsWoGroups.NextChild (pObj))
	{
		const CDXGraphic		*grobj = dynamic_cast<const CDXGraphic*>(pObj);

		if (grobj != NULL && grobj->m_graphicType == kCDXGraphicType_Symbol)
		{
			if (grobj->m_symbolType == kCDXSymbolType_Absolute)
				hasAbsMarker = true;
			else if (grobj->m_symbolType == kCDXSymbolType_Racemic)
				hasRacMarker = true;
			else if (grobj->m_symbolType == kCDXSymbolType_Relative)
				hasRelMarker = true;
		}
	} // while arrows
	if ((hasRacMarker && hasRelMarker) || (hasRacMarker && hasAbsMarker) || (hasAbsMarker && hasRelMarker))
	{
		// Can't have more than one marker, so let's ignore them all
		hasAbsMarker = false;
		hasRacMarker = false;
		hasRelMarker = false;
	}

	return hasAbsMarker || hasRacMarker || hasRelMarker;
}

/*
+===========================================================================================+
|																							|
| GetIDOfExternalBond	For a bond connected to an external connection point within a		|
|						fragment, return the ID of the corresponding external bond			|
|						connecting to the fragment (else return the ID of this bond)		|
|																							|
+===========================================================================================+
*/
CDXObjectID GetIDOfExternalBond(const CDXBond *pObjBond, const CDXBond **pObjExternalBond, const CDXNode *containingNodeOverride)
{
	const CDXObjectID id = pObjBond->GetObjectID();
	if (pObjExternalBond != NULL)
		*pObjExternalBond = pObjBond;
	
	if (pObjBond->m_beginNode == NULL || pObjBond->m_endNode == NULL)
		return id;

	if (pObjBond->m_beginNode->m_nodeType != kCDXNodeType_ExternalConnectionPoint
		&& pObjBond->m_endNode->m_nodeType != kCDXNodeType_ExternalConnectionPoint)
	{
		return id;
	}

	const CDXFragment *frag = dynamic_cast<const CDXFragment *>(pObjBond->GetParent());
	if (frag != NULL)
	{
		const CDXNode *node = containingNodeOverride;
		if (node == NULL)
			node = dynamic_cast<const CDXNode *>(pObjBond->GetParent()->GetParent());
		if (node != NULL && !node->m_alphaBonds.empty())
		{
			if (frag->m_connectionOrdering.size() == node->m_alphaBonds.size())
			{
				for (int i = 0; i < frag->m_connectionOrdering.size(); ++i)
				{
					// The fragment's connection ordering is the order of the attachment points, not the bonds connected to them
					const CDXNode *attachmentNode = dynamic_cast<const CDXNode *>(frag->FindByID(frag->m_connectionOrdering[i]));
					if (attachmentNode != NULL && !attachmentNode->m_alphaBonds.empty())
					{
						ASSERT(attachmentNode->m_alphaBonds.size() == 1);	// can it ever be more than 1?

						if (pObjBond->GetObjectID() == attachmentNode->m_alphaBonds[0]->GetObjectID())
						{
							// OK, we've found a matching internal bond, and we know it's in the i'th position.
							// Let's find the i'th corresponding external bond.

							if (node->m_bondOrdering != NULL && node->m_bondOrdering->size() > i)
							{
								for (vector<CDXBond *>::const_iterator alphaBond = node->m_alphaBonds.begin(); alphaBond != node->m_alphaBonds.end(); ++alphaBond)
									if ((*alphaBond)->GetObjectID() == (*node->m_bondOrdering)[i])
										return GetIDOfExternalBond(*alphaBond, pObjExternalBond);
							}

							return GetIDOfExternalBond(node->m_alphaBonds[i], pObjExternalBond);
						}
					}
				}

				// If we expected to be able to figure things out from the connection ordering, it would be really strange to fail.
				// If we do fail, we might as well return the original id.
				return id;
			}

			CDXObjectsRange	cont_bonds = frag->ContainedObjects (kCDXObj_Bond);
			int attachNum = 0;
			for (CDXObjectsByTag::const_iterator  itBond = cont_bonds.begin();  itBond != cont_bonds.end();  itBond++)
			{
				const CDXBond *b = dynamic_cast<CDXBond*>(itBond->second);
				if (b == NULL)
					continue;
				if (b == pObjBond)
				{
					if (attachNum >= node->m_alphaBonds.size())
					{
						return id;
					}
					else 
					{
						if (node->m_bondOrdering != NULL && attachNum < node->m_bondOrdering->size())
						{
							CDXObjectID bondID = (*node->m_bondOrdering)[attachNum];
							for (vector<CDXBond *>::const_iterator alphaBond = node->m_alphaBonds.begin(); alphaBond != node->m_alphaBonds.end(); ++alphaBond)
								if ((*alphaBond)->GetObjectID() == bondID)
									return GetIDOfExternalBond(*alphaBond, pObjExternalBond);
						}
						
						return GetIDOfExternalBond(node->m_alphaBonds[attachNum], pObjExternalBond);
					}
				}
				if (b->m_beginNode->m_nodeType == kCDXNodeType_ExternalConnectionPoint || b->m_endNode->m_nodeType == kCDXNodeType_ExternalConnectionPoint)
					++attachNum;
			}

			if (pObjExternalBond != NULL)
				*pObjExternalBond = node->m_alphaBonds[attachNum];
			return node->m_alphaBonds[attachNum]->GetObjectID();
		}
	}
	return id;
}

/*
+===========================================================================================+
|																							|
| GetInteriorConnectionForNode	For a node representing a contracted fragment, find the		|
|								interior node to which an exterior bond is logically		|
|								connected.													|
|								This function is recursive because CDX allows unlimited		|
|								levels of nested fragments.															|
|																							|
+===========================================================================================+
*/
CDXNode *GetInteriorConnectionForNode(map<const CDXNode *, vector<CDXNode *> > &mapAtToInteriorConnections, CDXNode *externalNode, const CDXBond *pObjBond)
{
	if (mapAtToInteriorConnections.find(externalNode) != mapAtToInteriorConnections.end())
	{
		const CDXObjectID id = pObjBond->GetObjectID();
		std::vector<CDXNode *> &interiorConnections = mapAtToInteriorConnections[externalNode];
		ASSERT(interiorConnections.size() == externalNode->m_alphaBonds.size());
		size_t maxConnection = min(interiorConnections.size(), externalNode->m_alphaBonds.size());
		for (size_t i = 0; i < maxConnection; ++i)
		{
			if (externalNode->m_bondOrdering != NULL && !externalNode->m_bondOrdering->empty())
			{
				if (i < externalNode->m_bondOrdering->size() && (*externalNode->m_bondOrdering)[i] == id)
					return GetInteriorConnectionForNode(mapAtToInteriorConnections, interiorConnections[i], externalNode->m_alphaBonds[i]);
			}
			else if (externalNode->m_alphaBonds[i] == pObjBond)
			{
				ASSERT(!interiorConnections[i]->m_alphaBonds.empty());
				if (interiorConnections[i]->m_alphaBonds.empty())
					return externalNode;
				if (interiorConnections[i]->m_alphaBonds.size() > i)
					return GetInteriorConnectionForNode(mapAtToInteriorConnections, interiorConnections[i], interiorConnections[i]->m_alphaBonds[i]);
				else
					return GetInteriorConnectionForNode(mapAtToInteriorConnections, interiorConnections[i], interiorConnections[i]->m_alphaBonds[0]);
			}
		}
	}

	return externalNode;
}

CDXGroup* BreakDiscontinuousSetsIntoFragments(
	CDXFragment* topLevelFragment,
	CDXObjectID& nextID)
{
	// Recursive algorithm which maps atom and bonds of discontinuous sets
	// to unique fragment objects.
	CDXObject* pPageOrGroupContainer = topLevelFragment->GetParent();
	ASSERT (pPageOrGroupContainer->GetTag() == kCDXObj_Page  ||
			pPageOrGroupContainer->GetTag() == kCDXObj_Group);
	const bool	fContainerIsGroup = pPageOrGroupContainer->GetTag() == kCDXObj_Group;

	/*
	+---------------------------------------------------------------------------------------+
	|	Part 1.	Identify the molecules, i.e. discrete connected atoms/bonds, and for each	|
	|			create a CDXFragment, store that in fragmentList, and record (in			|
	|			objectToFragmentMap) all the CDX objects that are associated with that		|
	|			fragment.																	|
	+---------------------------------------------------------------------------------------+
	*/
	vector<CDXFragment*> fragmentList;
	CDMap<CDXObject*, CDXFragment*> objectToFragmentMap;

	CDXObjectsRange nodes = topLevelFragment->ContainedObjects(kCDXObj_Node);
	for (CDXObjectsByTag::const_iterator nodeIter = nodes.begin();
		nodeIter != nodes.end();
		nodeIter++)
	{
		CDXNode* nextNode = (CDXNode*) nodeIter->second;

		if (!objectToFragmentMap.Contains(nextNode))
		{
			CDXFragment* newFragment = new CDXFragment(++nextID);
			fragmentList.push_back(newFragment);

			// Identify all the atoms and bonds, etc. that belong to this fragment.
			BreakDiscontinuousSetsIntoFragmentsHelper(
				nextNode,
				newFragment,
				objectToFragmentMap);
		}
	}

	// Abandon the process if there are less than 2 continuous sets
    // and there is no attached-data property is present.
	if ((fragmentList.size() < 2) &&
        (!ValidateFragmentHasAttachedDataProperty(topLevelFragment)))
	{
		if (fragmentList.size() == 1)
		{
			--nextID;
			delete fragmentList[0];
		}
		return NULL;
	}

	/*
	+---------------------------------------------------------------------------------------+
	|	Part 2.	Split the original fragment into several.									|
	+---------------------------------------------------------------------------------------+
	*/
	// Move children to a temporary group so they won't get cloned.
	CDXGroup tempGroup(0);
	topLevelFragment->TransferChildrenTo(&tempGroup);

	// Does not delete the object.
	// The original top level fragment must be removed before the new top level group is added
	// because they have the same CDXObjectID and the page or group container cannot contain
	// two objects with the same CDXObjectID at the same time.
	pPageOrGroupContainer->RemoveChild(topLevelFragment, false);

	// Create a new top level group which will contain the fragment
	// objects and anything else that was contained by the original
	// top level fragment object.
	CDXGroup* topLevelGroup;
	if (fContainerIsGroup)
		topLevelGroup = (CDXGroup*)pPageOrGroupContainer;
	else	// it's a page
	{
		topLevelGroup = new CDXGroup(*topLevelFragment);
		pPageOrGroupContainer->AddChild(topLevelGroup);
	}

	DeleteAndNull(topLevelFragment); 

	// Make a separate list that will not be affected by removal during
	// iteration.
	vector<CDXObject*> reParentList;
	CDXObjectsRange objects = tempGroup.ContainedObjects();
	for (CDXObjectsByTag::const_iterator objectIter = objects.begin();
		objectIter != objects.end();
		++objectIter)
	{
		reParentList.push_back(objectIter->second);
	}

	// Reparent objects to the new top level group.
	for (int i = 0; i < reParentList.size(); i++)
	{
		CDXObject* nextObject = reParentList[i];

		CDXObject* newParent = topLevelGroup;
		if (objectToFragmentMap.Contains(nextObject))
			newParent = objectToFragmentMap[nextObject];

		tempGroup.TransferChildTo(
			nextObject,
			newParent);
	}

	// Transfer the fragment objects to the container.
	for (int i = 0; i < fragmentList.size(); i++)
		topLevelGroup->AddChild(fragmentList[i]);

    return topLevelGroup;
}

void BreakDiscontinuousSetsIntoFragmentsHelper
(
	CDXNode* nextNode,
	CDXFragment* assignedFragment,
	CDMap<CDXObject*, CDXFragment*>& objectToFragmentMap)
{
	if (!objectToFragmentMap.Contains(nextNode))
	{
		objectToFragmentMap.insert(nextNode, assignedFragment);

		for (int i = 0; i < nextNode->m_alphaBonds.size(); i++)
		{
			CDXBond* nextBond = nextNode->m_alphaBonds[i];            
			if (!objectToFragmentMap.Contains(nextBond))
			{
				objectToFragmentMap.insert(nextBond, assignedFragment);

                const bool followThisBond = nextBond->MustBothAtomsResideInSameFragment();
                if (followThisBond)
                {
                    // Recurse if this bond requires it
                    CDXNode* otherNode = nextBond->OtherNode(nextNode);                
                    BreakDiscontinuousSetsIntoFragmentsHelper(otherNode, assignedFragment, objectToFragmentMap);                                                                                        
                }
			}
		}
	}
}

void FixNestingOfBrackets(CDXPage *cdxPage)
{
	// Molfiles typically don't have appropriate bracket hierarchies
	CdxGetChildrenIgnoringGroups	brGroups(cdxPage, kCDXObj_BracketedGroup);
	vector<CDXBracketedGroup *> bracketedGroups;
	const CDXObject *pObj;
	while (brGroups.NextChild(pObj))
	{
		const CDXBracketedGroup *pBracketedGroup = dynamic_cast<const CDXBracketedGroup*>(pObj);
		bracketedGroups.push_back(const_cast<CDXBracketedGroup *>(pBracketedGroup));
	}

	if (bracketedGroups.size() <= 1)
		return;

	for (vector<CDXBracketedGroup *>::iterator bracketedGroupInner = bracketedGroups.begin(); bracketedGroupInner != bracketedGroups.end(); ++bracketedGroupInner)
	{
		if (*bracketedGroupInner == NULL)
			continue;

		CDXBracketedGroup *bestOuter = NULL;
		for (vector<CDXBracketedGroup *>::iterator bracketedGroupOuter = bracketedGroups.begin(); bracketedGroupOuter != bracketedGroups.end(); ++bracketedGroupOuter)
		{
			if (*bracketedGroupOuter == NULL)
				continue;

			// Can't nest inside itself
			if (*bracketedGroupOuter == *bracketedGroupInner)
				continue;

			// Can't nest inside a smaller group
			if ((*bracketedGroupOuter)->m_bracketedObjects.size() < (*bracketedGroupInner)->m_bracketedObjects.size())
				continue;

			// Must nest inside the smallest container
			if (bestOuter != NULL && bestOuter->m_bracketedObjects.size() < (*bracketedGroupOuter)->m_bracketedObjects.size())
				continue;

			bool allInnerObjectsContainedInOuter = true;
			for (vector<CDXObjectID>::const_iterator id = (*bracketedGroupInner)->m_bracketedObjects.begin(); id != (*bracketedGroupInner)->m_bracketedObjects.end(); ++id)
				if (find((*bracketedGroupOuter)->m_bracketedObjects.begin(), (*bracketedGroupOuter)->m_bracketedObjects.end(), *id) == (*bracketedGroupOuter)->m_bracketedObjects.end())
				{
					allInnerObjectsContainedInOuter = false;
					break;
				}
			if (allInnerObjectsContainedInOuter)
				bestOuter = *bracketedGroupOuter;
		}

		if (bestOuter != NULL)
		{
			(*bracketedGroupInner)->GetParent()->TransferChildTo(*bracketedGroupInner, bestOuter);
			*bracketedGroupInner = NULL;
		}
	}
}

/*
+===========================================================================================+
|																							|
| MoveObjects_Recurse	Moves the specified CDXObject along with all the contained objects  |
|						Checks for objects of type CDXGraphicObject having a known Position |
|						or BoundingBox.														|
|																							|
+===========================================================================================+
*/
void MoveObjects_Recurse(CDXObject* objRoot, CDXCoordinate xShift, CDXCoordinate yShift)
{
	CDXObjectsRange range = objRoot->ContainedObjects();
	for(CDXObjectsByTag::const_iterator it = range.begin(); it != range.end(); it++)
	{
		CDXObject* obj = it->second;
		if(CDXGraphicObject* go = dynamic_cast<CDXGraphicObject*>(obj))
		{
			// move any CDXGraphicsObject having Position or BoundingBox
			if(go->KnownPosition())
			{
				CDXPoint2D pos = go->Position();
				pos.Offset(xShift, yShift);
				go->Position(pos);
			}
			else if(go->KnownBoundingBox())
			{
				CDXRectangle bb = go->BoundingBox();
				bb.Offset(xShift, yShift);
				go->BoundingBox(bb);
			}

			// special processing for certain objects that do not rely only on position or bounding box
			// and have their own custom position data
			if(CDXCurve* curve = dynamic_cast<CDXCurve*>(go))
			{
				// shift each curve point by shift values
				for(std::vector<CDXPoint2D>::iterator it = curve->m_curvePoints.begin(); it != curve->m_curvePoints.end(); it++)
					it->Offset(xShift, yShift);
			}
			else if(CDXGraphic* graphic = dynamic_cast<CDXGraphic*>(go))
			{
				graphic->m_3dHead.Offset(xShift, yShift, 0);
				graphic->m_3dTail.Offset(xShift, yShift, 0);
				graphic->m_center.Offset(xShift, yShift, 0);
				graphic->m_majorAxisEnd.Offset(xShift, yShift, 0);
				graphic->m_minorAxisEnd.Offset(xShift, yShift, 0);
			}
			else if(CDXArrow* arrow = dynamic_cast<CDXArrow*>(go))
			{
				arrow->m_3dHead.Offset(xShift, yShift, 0);
				arrow->m_3dTail.Offset(xShift, yShift, 0);
				arrow->m_center.Offset(xShift, yShift, 0);
				arrow->m_majorAxisEnd.Offset(xShift, yShift, 0);
				arrow->m_minorAxisEnd.Offset(xShift, yShift, 0);
			}
			else if(CDXNamedAlternativeGroup* nag = dynamic_cast<CDXNamedAlternativeGroup*>(go))
			{
				CDXRectangle rectTF = nag->TextFrame();
				rectTF.Offset(xShift, yShift);
				nag->TextFrame(rectTF);

				CDXRectangle rectGF = nag->GroupFrame();
				rectGF.Offset(xShift, yShift);
				nag->GroupFrame(rectGF);
			}
			else if(CDXPlateBase* plate = dynamic_cast<CDXPlateBase*>(go))
			{
				plate->m_topleft.Offset(xShift, yShift);
				plate->m_topright.Offset(xShift, yShift);
				plate->m_bottomleft.Offset(xShift, yShift);
				plate->m_bottomright.Offset(xShift, yShift);
			}
		}
		MoveObjects_Recurse(obj, xShift, yShift);
	}
}

/*
* GetObjectTags
* Returns all object tags for the specified object where the object tag name matches the specified name.
*/
std::vector<CDXObjectTag*> GetObjectTags(const CDXObject& object, const std::string& tagName)
{
	std::vector<CDXObjectTag*> objectTags;
	const CDXObjectsByTag* contents = object.GetContents();
	if (contents)
	{
		for (CDXObjectsByTag::const_iterator itContents = contents->begin(); itContents != contents->end(); ++itContents)
		{
			CDXObjectTag* tag = dynamic_cast<CDXObjectTag*>(itContents->second);
			if (tag)
			{
				// Return all tags if tag name is not specified
				if (tagName.empty() || (tag->m_Name == tagName))
				{
					objectTags.push_back(tag);
				}
			}
		}
	}
	return objectTags;
}

/*
* GetAllObjectTagsOnBasisObjects
* Returns all object tags on specified property's basis objects.
*/
std::vector<CDXObjectTag*> GetAllObjectTagsOnBasisObjects(const CDXChemicalProperty& prop, const CDXObject& cdxPage)
{
	std::vector<CDXObjectTag*> objectTags;
	for (std::vector<CDXObjectID>::const_iterator itBasis = prop.m_basisObjects.begin(); itBasis != prop.m_basisObjects.end(); ++itBasis)
	{
		CDXObject* object = cdxPage.FindByIDRecursive(*itBasis);
		if (object)
		{
			std::vector<CDXObjectTag*> tagsForObject = GetObjectTags(*object, "");
			if (!tagsForObject.empty())
			{
				objectTags.insert(objectTags.end(), tagsForObject.begin(), tagsForObject.end());
			}
		}
	}

	return objectTags;
}

const CDXPage* GetCDXPageForObject(const CDXObject& object)
{
	const CDXObject* parent = object.GetParent();
	while (parent)
	{
		if (dynamic_cast<const CDXPage*>(parent))
		{
			return (const CDXPage*)parent;
		}

		parent = parent->GetParent();
	}

	return NULL;
}
