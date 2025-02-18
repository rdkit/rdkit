// CommonCS/LibCommon/Hdr/CDXUtils.h
// Copyright © 1997-2007, CambridgeSoft Corp., All Rights Reserved

#pragma once

#include "CoreChemistryAPI.h"
#include "CDXStdObjects.h"
#include "CDMap.h"

// CDXIsValid performs a cursory check on whether the given object, and recursively its children, look sensible.
bool			CDXIsValid					(const CDXObject* pTopOfObjectTree, std::set<const CDXObject*>* pAlreadyChecked = NULL);	// typically pObj is a CDXDocument*

CDXRectangle	CDXReadOrComputeBoundingBox (const CDXObject* pObj, bool fIgnoreCachedInfo = false, bool fIgnoreObjectTags = false);
void			CDXBuildObjectMap			(const CDXObject& top, map<CDXObjectID,const CDXObject*> &objMap, bool* pDuplicateIdsPresent = NULL, bool fIgnoreIdZero = true);
CDXPoint2D		CDXGetTextExtent			(const CDXString& str);
double			CdxFindMedianBondLength		(const CDXObject &top, Dimensions2or3 dimension = kIn2D, double *pMean = NULL);

CORE_CHEMISTRY_API bool AnyChargeObjects(const CDXObject *cdxObj);

bool			AnyAbsRacRelMarkers			(const CDXObject *cdxObj, bool &hasAbsMarker, bool &hasRacMarker, bool &hasRelMarker);

CDXObjectID		GetIDOfExternalBond			(const CDXBond *pObjBond, const CDXBond **pObjExternalBond = NULL, const CDXNode *containingNodeOverride = NULL);
CDXNode *		GetInteriorConnectionForNode(map<const CDXNode *, vector<CDXNode *> > &mapAtToInteriorConnections, CDXNode *externalNode, const CDXBond *pObjBond);

CDXGroup*			BreakDiscontinuousSetsIntoFragments(
					CDXFragment* topLevelFragment,
					CDXObjectID& nextID);

void			BreakDiscontinuousSetsIntoFragmentsHelper(
					CDXNode* nextNode,
					CDXFragment* assignedFragment,
					CDMap<CDXObject*, CDXFragment*>& objectToFragmentMap);

void			FixNestingOfBrackets(CDXPage *cdxPage);

void			MoveObjects_Recurse(CDXObject* objRoot, CDXCoordinate xShift, CDXCoordinate yShift);

std::vector<CDXObjectTag*> GetObjectTags(const CDXObject& object, const std::string& tagName);
std::vector<CDXObjectTag*> GetAllObjectTagsOnBasisObjects(const CDXChemicalProperty& prop, const CDXObject& cdxPage);

CORE_CHEMISTRY_API const CDXPage* GetCDXPageForObject(const CDXObject& object);
