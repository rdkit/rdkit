// $Id$
//
//                 Copyright 2001-2006
//                   Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include "types.h"

namespace RDKit {
namespace common_properties {
const char * propnames[] = {
  "_Name",
  "_MolFileInfo",
  "_MolFileComments",
  "_2DConf",
  "_3DConf",
  "_doIsoSmiles",
  "extraRings",
  "_smilesAtomOutputOrder",
  "_StereochemDone",
  "_NeedsQueryScan",
  "_fragSMARTS",
  "maxAttachIdx",
  "origNoImplicit",
  "ringMembership",
  "_connectivityHKDeltas",
  "_connectivityNVals",
  "_crippenLogP",
  "_crippenLogPContribs",
  "_crippenMR",
  "_crippenMRContribs",
  "_labuteASA",
  "_labuteAtomContribs",
  "_labuteAtomHContrib",
  "_tpsa",
  "_tpsaAtomContribs",
  "numArom",
  "_MMFFSanitized",
  "_CrippenLogP",
  "_CrippenMR",
  "_BondsPotentialStereo",
  "_CIPCode",
  "_CIPRank",
  "_ChiralityPossible",
  "_UnknownStereo",
  "_ringStereoAtoms",
  "_ringStereochemCand",
  "_ringStereoWarning",
  "_SmilesStart",
  "_TraversalBondIndexOrder",
  "_TraversalRingClosureBond",
  "_TraversalStartPoint",
  "_queryRootAtom",
  "_hasMassQuery",
  "_protected",
  "_supplementalSmilesLabel",
  "_unspecifiedOrder",
  "_RingClosures",
  "molAtomMapNumber",
  "molFileAlias",
  "molFileValue",
  "molInversionFlag",
  "molParity",
  "molRxnComponent",
  "molRxnRole",
  "molTotValence",
  "_MolFileRLabel",
  "_MolFileChiralFlag",
  "dummyLabel",
  "_QueryFormalCharge",
  "_QueryHCount",
  "_QueryIsotope",
  "_QueryMass",
  "_ReactionDegreeChanged",
  "NullBond",
  "_rgroupAtomMaps",
  "_rgroupBonds",
  "_AtomID",
  "_starred",
  "_SLN_s",
  "_Unfinished_SLN_",
  "_brokenChirality",
  "isImplicit",
  "smilesSymbol",
  "_TriposAtomType",
  "2D",
  "BalabanJ",
  "BalanbanJ",
  "Discrims",
  "DistanceMatrix_Paths",
  //"__computedProps"
};


}  // end common_properties

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double EPS_DOUBLE = std::numeric_limits<double>::epsilon();
const double SMALL_DOUBLE = 1.0e-8;
const double MAX_INT = static_cast<double>(std::numeric_limits<int>::max());
const double MAX_LONGINT =
    static_cast<double>(std::numeric_limits<LONGINT>::max());

//  template <typename T>
//  T larger_of(T arg1,T arg2) { return arg1>arg2 ? arg1 : arg2; };

double round(double num) {
  double floorVal = floor(num);
  double ceilVal = ceil(num);
  return num - floorVal > ceilVal - num ? ceilVal : floorVal;
};

void Union(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
  res.resize(0);
  res = r1;
  INT_VECT_CI ri;
  for (ri = r2.begin(); ri != r2.end(); ri++) {
    if (std::find(res.begin(), res.end(), (*ri)) == res.end()) {
      res.push_back(*ri);
    }
  }
}

void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
  res.resize(0);
  INT_VECT_CI ri;
  for (ri = r1.begin(); ri != r1.end(); ri++) {
    if (std::find(r2.begin(), r2.end(), (*ri)) != r2.end()) {
      res.push_back(*ri);
    }
  }
}

void Union(const VECT_INT_VECT &rings, INT_VECT &res, const INT_VECT *exclude) {
  res.resize(0);
  INT_VECT ring;
  unsigned int id;
  unsigned int nrings = static_cast<unsigned int>(rings.size());
  INT_VECT_CI ri;

  for (id = 0; id < nrings; id++) {
    if (exclude) {
      if (std::find(exclude->begin(), exclude->end(), static_cast<int>(id)) !=
          exclude->end()) {
        continue;
      }
    }
    ring = rings[id];
    for (ri = ring.begin(); ri != ring.end(); ri++) {
      if (std::find(res.begin(), res.end(), (*ri)) == res.end()) {
        res.push_back(*ri);
      }
    }
  }
}

int nextCombination(INT_VECT &comb, int tot) {
  int nelem = static_cast<int>(comb.size());
  int celem = nelem - 1;

  while (comb[celem] == (tot - nelem + celem)) {
    celem--;
    if (celem < 0) {
      return -1;
    }
  }

  unsigned int i;
  comb[celem] += 1;
  for (i = celem + 1; i < comb.size(); i++) {
    comb[i] = comb[i - 1] + 1;
  }
  return celem;
}
}
