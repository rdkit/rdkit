//
//  Copyright 2001-2021 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_TYPES_H
#define RD_TYPES_H

#ifdef WIN32
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#define _DEFINED_USE_MATH_DEFINES
#endif
#endif
#include <cmath>
#ifdef _DEFINED_USE_MATH_DEFINES
#undef _DEFINED_USE_MATH_DEFINES
#undef _USE_MATH_DEFINES
#endif

#include "Invariant.h"
#include "Dict.h"

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <numeric>
#include <list>
#include <limits>

#include <cstring>
#include <any>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace detail {
// used in various places for computed properties
inline const std::string computedPropName = "__computedProps";
}  // namespace detail

namespace common_properties {
///////////////////////////////////////////////////////////////
// Molecule Props
inline const std::string TWOD = "2D";
inline const std::string BalabanJ = "BalabanJ";
inline const std::string BalanbanJ = "BalanbanJ";
inline const std::string Discrims = "Discrims";
inline const std::string DistanceMatrix_Paths = "DistanceMatrix_Paths";
inline const std::string MolFileComments = "MolFileComments";
inline const std::string MolFileInfo = "MolFileInfo";
inline const std::string NullBond = "NullBond";
inline const std::string _2DConf = "_2DConf";
inline const std::string _3DConf = "_3DConf";
inline const std::string _AtomID = "_AtomID";
inline const std::string _BondsPotentialStereo = "_BondsPotentialStereo";
inline const std::string _ChiralAtomRank = "_chiralAtomRank";
inline const std::string _CIPCode = "_CIPCode";
inline const std::string _CIPRank = "_CIPRank";
inline const std::string _CIPComputed = "_CIPComputed";
inline const std::string _CanonicalRankingNumber = "_CanonicalRankingNumber";
inline const std::string _ChiralityPossible = "_ChiralityPossible";
inline const std::string _CrippenLogP = "_CrippenLogP";
inline const std::string _CrippenMR = "_CrippenMR";
inline const std::string _MMFFSanitized = "_MMFFSanitized";
inline const std::string _MolFileChiralFlag = "_MolFileChiralFlag";
inline const std::string MRV_SMA = "MRV SMA";
inline const std::string _MolFileRLabel = "_MolFileRLabel";
inline const std::string _MolFileAtomQuery = "_MolFileAtomQuery";
inline const std::string _MolFileBondQuery = "_MolFileBondQuery";
inline const std::string _MolFileBondEndPts = "_MolFileBondEndPts";
inline const std::string _MolFileBondAttach = "_MolFileBondAttach";
inline const std::string _MolFileBondType = "_MolFileBondType";
inline const std::string _MolFileBondStereo = "_MolFileBondStereo";
inline const std::string _MolFileBondCfg = "_MolFileBondCfg";

inline const std::string _Name = "_Name";
inline const std::string _NeedsQueryScan = "_NeedsQueryScan";
inline const std::string _NonExplicit3DChirality = "_NonExplicit3DChirality";
inline const std::string _QueryFormalCharge = "_QueryFormalCharge";
inline const std::string _QueryHCount = "_QueryHCount";
inline const std::string _QueryIsotope = "_QueryIsotope";
inline const std::string _QueryMass = "_QueryMass";
inline const std::string _ReactionDegreeChanged = "_ReactionDegreeChanged";
inline const std::string reactantAtomIdx = "react_atom_idx";
inline const std::string reactionMapNum = "old_mapno";
inline const std::string reactantIdx = "react_idx";

inline const std::string _RingClosures = "_RingClosures";
inline const std::string _SLN_s = "_SLN_s";
inline const std::string _SmilesStart = "_SmilesStart";
inline const std::string _StereochemDone = "_StereochemDone";
inline const std::string _TraversalBondIndexOrder = "_TraversalBondIndexOrder";
inline const std::string _TraversalRingClosureBond =
    "_TraversalRingClosureBond";
inline const std::string _TraversalStartPoint = "_TraversalStartPoint";
inline const std::string _TriposAtomType = "_TriposAtomType";
inline const std::string _Unfinished_SLN_ = "_Unfinished_SLN_";
inline const std::string _UnknownStereo = "_UnknownStereo";
inline const std::string _connectivityHKDeltas = "_connectivityHKDeltas";
inline const std::string _connectivityNVals = "_connectivityNVals";
inline const std::string _crippenLogP = "_crippenLogP";
inline const std::string _crippenLogPContribs = "_crippenLogPContribs";
inline const std::string _crippenMR = "_crippenMR";
inline const std::string _crippenMRContribs = "_crippenMRContribs";
inline const std::string _GasteigerCharge = "_GasteigerCharge";
inline const std::string _GasteigerHCharge = "_GasteigerHCharge";
inline const std::string _doIsoSmiles = "_doIsoSmiles";
inline const std::string _fragSMARTS = "_fragSMARTS";
inline const std::string _hasMassQuery = "_hasMassQuery";
inline const std::string _labuteASA = "_labuteASA";
inline const std::string _labuteAtomContribs = "_labuteAtomContribs";
inline const std::string _labuteAtomHContrib = "_labuteAtomHContrib";
inline const std::string _protected = "_protected";
inline const std::string _queryRootAtom = "_queryRootAtom";
inline const std::string _ringStereoAtoms = "_ringStereoAtoms";
inline const std::string _ringStereoWarning = "_ringStereoWarning";
inline const std::string _ringStereochemCand = "_ringStereochemCand";
inline const std::string _ringStereoOtherAtom = "_ringStereoOtherAtom";
inline const std::string _mesoOtherAtom = "_mesoOtherAtom";
inline const std::string _chiralPermutation = "_chiralPermutation";
inline const std::string _smilesAtomOutputOrder = "_smilesAtomOutputOrder";
inline const std::string _smilesBondOutputOrder = "_smilesBondOutputOrder";
inline const std::string _starred = "_starred";
inline const std::string _supplementalSmilesLabel = "_supplementalSmilesLabel";
inline const std::string _tpsa = "_tpsa";
inline const std::string _tpsaAtomContribs = "_tpsaAtomContribs";
inline const std::string _unspecifiedOrder = "_unspecifiedOrder";
inline const std::string _brokenChirality = "_brokenChirality";
inline const std::string _rgroupAtomMaps = "_rgroupAtomMaps";
inline const std::string _rgroupBonds = "_rgroupBonds";
inline const std::string _rgroupTargetAtoms = "_rgroupTargetAtoms";
inline const std::string _rgroupTargetBonds = "_rgroupTargetBonds";
inline const std::string dummyLabel = "dummyLabel";
inline const std::string extraRings = "extraRings";
inline const std::string isImplicit = "isImplicit";
inline const std::string maxAttachIdx = "maxAttachIdx";
inline const std::string molAtomMapNumber = "molAtomMapNumber";
inline const std::string molFileAlias = "molFileAlias";
inline const std::string molFileValue = "molFileValue";
inline const std::string molInversionFlag = "molInversionFlag";
inline const std::string molParity = "molParity";
inline const std::string molStereoCare = "molStereoCare";
inline const std::string molRxnComponent = "molRxnComponent";
inline const std::string molRxnRole = "molRxnRole";
inline const std::string molTotValence = "molTotValence";
inline const std::string molFileLinkNodes = "_molLinkNodes";
inline const std::string numArom = "numArom";
inline const std::string ringMembership = "ringMembership";
inline const std::string smilesSymbol = "smilesSymbol";
inline const std::string atomLabel = "atomLabel";
inline const std::string OxidationNumber = "OxidationNumber";
inline const std::string internalRgroupSmiles = "internalRgroupSmiles";
inline const std::string molRingBondCount = "molRingBondCount";
inline const std::string molSubstCount = "molSubstCount";
inline const std::string molAttachPoint = "molAttchpt";
inline const std::string molAttachOrder = "molAttchord";
inline const std::string molAttachOrderTemplate = "molAttachOrderTemplate";
inline const std::string molAtomClass = "molClass";
inline const std::string molAtomSeqId = "molSeqid";
inline const std::string molRxnExactChange = "molRxnExachg";
inline const std::string molReactStatus = "molReactStatus";
inline const std::string _fromAttachPoint = "_fromAttchpt";
inline const std::string natReplace = "natReplace";
inline const std::string templateNames = "templateNames";

inline const std::string molNote = "molNote";
inline const std::string atomNote = "atomNote";
inline const std::string bondNote = "bondNote";
inline const std::string _isotopicHs = "_isotopicHs";

inline const std::string _QueryAtomGenericLabel = "_QueryAtomGenericLabel";

// molecule drawing
inline const std::string _displayLabel = "_displayLabel";
inline const std::string _displayLabelW = "_displayLabelW";

}  // namespace common_properties
#ifndef WIN32
typedef long long int LONGINT;
#else
typedef __int64 LONGINT;
#endif
#ifdef max
#undef max  // FUCK I hate this nonsense
#endif
#ifdef min
#undef min  // FUCK I hate this nonsense
#endif

RDKIT_RDGENERAL_EXPORT extern const double MAX_DOUBLE;
RDKIT_RDGENERAL_EXPORT extern const double EPS_DOUBLE;
RDKIT_RDGENERAL_EXPORT extern const double SMALL_DOUBLE;
RDKIT_RDGENERAL_EXPORT extern const double MAX_INT;
RDKIT_RDGENERAL_EXPORT extern const double MAX_LONGINT;

typedef unsigned int UINT;
typedef unsigned short USHORT;
typedef unsigned char UCHAR;

typedef std::vector<int> INT_VECT;
typedef INT_VECT::iterator INT_VECT_I;
typedef INT_VECT::const_iterator INT_VECT_CI;
typedef INT_VECT::reverse_iterator INT_VECT_RI;
typedef INT_VECT::const_reverse_iterator INT_VECT_CRI;

typedef std::list<int> INT_LIST;
typedef INT_LIST::iterator INT_LIST_I;
typedef INT_LIST::const_iterator INT_LIST_CI;

typedef std::list<INT_VECT> LIST_INT_VECT;
typedef LIST_INT_VECT::iterator LIST_INT_VECT_I;
typedef LIST_INT_VECT::const_iterator LIST_INT_VECT_CI;

typedef std::vector<INT_VECT> VECT_INT_VECT;
typedef VECT_INT_VECT::iterator VECT_INT_VECT_I;
typedef VECT_INT_VECT::const_iterator VECT_INT_VECT_CI;

typedef std::vector<UINT>::const_iterator UINT_VECT_CI;
typedef std::vector<UINT> UINT_VECT;

typedef std::vector<std::string>::const_iterator STR_VECT_CI;
typedef std::vector<std::string>::iterator STR_VECT_I;
typedef std::vector<std::string> STR_VECT;

typedef std::vector<double> DOUBLE_VECT;
typedef DOUBLE_VECT::iterator DOUBLE_VECT_I;
typedef DOUBLE_VECT::const_iterator DOUBLE_VECT_CI;
typedef std::vector<DOUBLE_VECT> VECT_DOUBLE_VECT;
typedef VECT_DOUBLE_VECT::iterator VECT_DOUBLE_VECT_I;
typedef VECT_DOUBLE_VECT::const_iterator VECT_DOUBLE_VECT_CI;

typedef std::map<std::string, UINT> STR_UINT_MAP;
typedef std::map<std::string, UINT>::const_iterator STR_UINT_MAP_CI;

typedef std::map<int, INT_VECT> INT_INT_VECT_MAP;
typedef INT_INT_VECT_MAP::const_iterator INT_INT_VECT_MAP_CI;

typedef std::map<int, int> INT_MAP_INT;
typedef INT_MAP_INT::iterator INT_MAP_INT_I;
typedef INT_MAP_INT::const_iterator INT_MAP_INT_CI;

typedef std::deque<int> INT_DEQUE;
typedef INT_DEQUE::iterator INT_DEQUE_I;
typedef INT_DEQUE::const_iterator INT_DEQUE_CI;

typedef std::map<int, INT_DEQUE> INT_INT_DEQ_MAP;
typedef INT_INT_DEQ_MAP::const_iterator INT_INT_DEQ_MAP_CI;

typedef std::set<int> INT_SET;
typedef INT_SET::iterator INT_SET_I;
typedef INT_SET::const_iterator INT_SET_CI;

//! functor to compare two doubles with a tolerance
struct RDKIT_RDGENERAL_EXPORT ltDouble {
 public:
  ltDouble() {}
  bool operator()(double d1, double d2) const {
    if (fabs(d1 - d2) < _tol) {
      return false;
    } else {
      return (d1 < d2);
    }
  }

 private:
  double _tol{1.0e-8};
};

//! std::map from double to integer.
typedef std::map<double, int, ltDouble> DOUBLE_INT_MAP;

//! functor for returning the larger of two values
template <typename T>
struct RDKIT_RDGENERAL_EXPORT larger_of {
  T operator()(T arg1, T arg2) { return arg1 > arg2 ? arg1 : arg2; }
};

//! functor for comparing two strings
struct RDKIT_RDGENERAL_EXPORT charptr_functor {
  bool operator()(const char *s1, const char *s2) const {
    return strcmp(s1, s2) < 0;
  }
};

//! \brief calculate the union of two INT_VECTs and put the results in a
//! third vector
RDKIT_RDGENERAL_EXPORT void Union(const INT_VECT &r1, const INT_VECT &r2,
                                  INT_VECT &res);

//! \brief calculate the intersection of two INT_VECTs and put the results in a
//! third vector
RDKIT_RDGENERAL_EXPORT void Intersect(const INT_VECT &r1, const INT_VECT &r2,
                                      INT_VECT &res);

//! calculating the union of the INT_VECT's in a VECT_INT_VECT
/*!
    \param rings   the INT_VECT's to consider
    \param res     used to return results
    \param exclude any values in this optional INT_VECT will be excluded
           from the union.
*/
RDKIT_RDGENERAL_EXPORT void Union(const VECT_INT_VECT &rings, INT_VECT &res,
                                  const INT_VECT *exclude = nullptr);

//! given a current combination of numbers change it to the next possible
// combination
/*!
  \param comb the <b>sorted</b> vector to consider
  \param tot the maximum number possible in the vector

  \return -1 on failure, the index of the last number changed on success.
  Example:
    for all combinations 3 of numbers between 0 and tot=5
    given (0,1,2) the function wil return (0,1,3) etc.


*/
RDKIT_RDGENERAL_EXPORT int nextCombination(INT_VECT &comb, int tot);
};  // namespace RDKit

#endif
