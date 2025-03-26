// CommonCS/LibCommon/Hdr/CDXUtils.h
// Copyright © 1997-2007, CambridgeSoft Corp., All Rights Reserved

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
