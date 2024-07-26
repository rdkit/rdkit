//
//  Copyright 2001-2021 Greg Landrum and other RDKit contributors
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
namespace detail {
const std::string computedPropName = "__computedProps";
}

namespace common_properties {
const std::string TWOD = "2D";
const std::string BalabanJ = "BalabanJ";
const std::string BalanbanJ = "BalanbanJ";
const std::string Discrims = "Discrims";
const std::string DistanceMatrix_Paths = "DistanceMatrix_Paths";
const std::string MolFileComments = "MolFileComments";
const std::string MolFileInfo = "MolFileInfo";
const std::string NullBond = "NullBond";
const std::string _2DConf = "_2DConf";
const std::string _3DConf = "_3DConf";
const std::string _AtomID = "_AtomID";
const std::string _BondsPotentialStereo = "_BondsPotentialStereo";
const std::string _ChiralAtomRank = "_chiralAtomRank";
const std::string _CIPCode = "_CIPCode";
const std::string _CIPRank = "_CIPRank";
const std::string _ChiralityPossible = "_ChiralityPossible";
const std::string _CrippenLogP = "_CrippenLogP";
const std::string _CrippenMR = "_CrippenMR";
const std::string _MMFFSanitized = "_MMFFSanitized";
const std::string _MolFileChiralFlag = "_MolFileChiralFlag";
const std::string MRV_SMA = "MRV SMA";
const std::string _MolFileRLabel = "_MolFileRLabel";
const std::string _MolFileAtomQuery = "_MolFileAtomQuery";
const std::string _MolFileBondQuery = "_MolFileBondQuery";
const std::string _MolFileBondEndPts = "_MolFileBondEndPts";
const std::string _MolFileBondAttach = "_MolFileBondAttach";
const std::string _MolFileBondType = "_MolFileBondType";
const std::string _MolFileBondStereo = "_MolFileBondStereo";
const std::string _MolFileBondCfg = "_MolFileBondCfg";

const std::string _Name = "_Name";
const std::string _NeedsQueryScan = "_NeedsQueryScan";
const std::string _NonExplicit3DChirality = "_NonExplicit3DChirality";
const std::string _QueryFormalCharge = "_QueryFormalCharge";
const std::string _QueryHCount = "_QueryHCount";
const std::string _QueryIsotope = "_QueryIsotope";
const std::string _QueryMass = "_QueryMass";
const std::string _ReactionDegreeChanged = "_ReactionDegreeChanged";
const std::string reactantAtomIdx = "react_atom_idx";
const std::string reactionMapNum = "old_mapno";

const std::string _RingClosures = "_RingClosures";
const std::string _SLN_s = "_SLN_s";
const std::string _SmilesStart = "_SmilesStart";
const std::string _StereochemDone = "_StereochemDone";
const std::string _TraversalBondIndexOrder = "_TraversalBondIndexOrder";
const std::string _TraversalRingClosureBond = "_TraversalRingClosureBond";
const std::string _TraversalStartPoint = "_TraversalStartPoint";
const std::string _TriposAtomType = "_TriposAtomType";
const std::string _Unfinished_SLN_ = "_Unfinished_SLN_";
const std::string _UnknownStereo = "_UnknownStereo";
const std::string _connectivityHKDeltas = "_connectivityHKDeltas";
const std::string _connectivityNVals = "_connectivityNVals";
const std::string _crippenLogP = "_crippenLogP";
const std::string _crippenLogPContribs = "_crippenLogPContribs";
const std::string _crippenMR = "_crippenMR";
const std::string _crippenMRContribs = "_crippenMRContribs";
const std::string _GasteigerCharge = "_GasteigerCharge";
const std::string _GasteigerHCharge = "_GasteigerHCharge";
const std::string _doIsoSmiles = "_doIsoSmiles";
const std::string _fragSMARTS = "_fragSMARTS";
const std::string _hasMassQuery = "_hasMassQuery";
const std::string _labuteASA = "_labuteASA";
const std::string _labuteAtomContribs = "_labuteAtomContribs";
const std::string _labuteAtomHContrib = "_labuteAtomHContrib";
const std::string _protected = "_protected";
const std::string _queryRootAtom = "_queryRootAtom";
const std::string _ringStereoAtoms = "_ringStereoAtoms";
const std::string _ringStereoWarning = "_ringStereoWarning";
const std::string _ringStereochemCand = "_ringStereochemCand";
const std::string _ringStereoOtherAtom = "_ringStereoOtherAtom";
const std::string _mesoOtherAtom = "_mesoOtherAtom";
const std::string _chiralPermutation = "_chiralPermutation";
const std::string _smilesAtomOutputOrder = "_smilesAtomOutputOrder";
const std::string _smilesBondOutputOrder = "_smilesBondOutputOrder";
const std::string _starred = "_starred";
const std::string _supplementalSmilesLabel = "_supplementalSmilesLabel";
const std::string _tpsa = "_tpsa";
const std::string _tpsaAtomContribs = "_tpsaAtomContribs";
const std::string _unspecifiedOrder = "_unspecifiedOrder";
const std::string _brokenChirality = "_brokenChirality";
const std::string _rgroupAtomMaps = "_rgroupAtomMaps";
const std::string _rgroupBonds = "_rgroupBonds";
const std::string dummyLabel = "dummyLabel";
const std::string extraRings = "extraRings";
const std::string isImplicit = "isImplicit";
const std::string maxAttachIdx = "maxAttachIdx";
const std::string molAtomMapNumber = "molAtomMapNumber";
const std::string molFileAlias = "molFileAlias";
const std::string molFileValue = "molFileValue";
const std::string molInversionFlag = "molInversionFlag";
const std::string molParity = "molParity";
const std::string molStereoCare = "molStereoCare";
const std::string molRxnComponent = "molRxnComponent";
const std::string molRxnRole = "molRxnRole";
const std::string molTotValence = "molTotValence";
const std::string molFileLinkNodes = "_molLinkNodes";
const std::string numArom = "numArom";
const std::string ringMembership = "ringMembership";
const std::string smilesSymbol = "smilesSymbol";
const std::string atomLabel = "atomLabel";
const std::string OxidationNumber = "OxidationNumber";
const std::string internalRgroupSmiles = "internalRgroupSmiles";
const std::string molRingBondCount = "molRingBondCount";
const std::string molSubstCount = "molSubstCount";
const std::string molAttachPoint = "molAttchpt";
const std::string molAttachOrder = "molAttchord";
const std::string molAtomClass = "molClass";
const std::string molAtomSeqId = "molSeqid";
const std::string molRxnExactChange = "molRxnExachg";
const std::string molReactStatus = "molReactStatus";
const std::string _fromAttachPoint = "_fromAttchpt";

const std::string molNote = "molNote";
const std::string atomNote = "atomNote";
const std::string bondNote = "bondNote";
const std::string _isotopicHs = "_isotopicHs";

const std::string _QueryAtomGenericLabel = "_QueryAtomGenericLabel";

// molecule drawing
const std::string _displayLabel = "_displayLabel";
const std::string _displayLabelW = "_displayLabelW";

}  // namespace common_properties

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double EPS_DOUBLE = std::numeric_limits<double>::epsilon();
const double SMALL_DOUBLE = 1.0e-8;
const double MAX_INT = static_cast<double>(std::numeric_limits<int>::max());
const double MAX_LONGINT =
    static_cast<double>(std::numeric_limits<LONGINT>::max());

//  template <typename T>
//  T larger_of(T arg1,T arg2) { return arg1>arg2 ? arg1 : arg2; };

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
  auto nrings = static_cast<unsigned int>(rings.size());
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
}  // namespace RDKit
