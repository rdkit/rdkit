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
#include <string_view>
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
inline constexpr std::string_view computedPropName = "__computedProps";
}  // namespace detail

namespace common_properties {
///////////////////////////////////////////////////////////////
// Molecule Props
inline constexpr std::string_view TWOD = "2D";
inline constexpr std::string_view BalabanJ = "BalabanJ";
inline constexpr std::string_view BalanbanJ = "BalanbanJ";
inline constexpr std::string_view Discrims = "Discrims";
inline constexpr std::string_view DistanceMatrix_Paths = "DistanceMatrix_Paths";
inline constexpr std::string_view MolFileComments = "MolFileComments";
inline constexpr std::string_view MolFileInfo = "MolFileInfo";
inline constexpr std::string_view NullBond = "NullBond";
inline constexpr std::string_view _2DConf = "_2DConf";
inline constexpr std::string_view _3DConf = "_3DConf";
inline constexpr std::string_view _AtomID = "_AtomID";
inline constexpr std::string_view _BondsPotentialStereo = "_BondsPotentialStereo";
inline constexpr std::string_view _ChiralAtomRank = "_chiralAtomRank";
inline constexpr std::string_view _CIPCode = "_CIPCode";
inline constexpr std::string_view _CIPRank = "_CIPRank";
inline constexpr std::string_view _CIPComputed = "_CIPComputed";
inline constexpr std::string_view _CanonicalRankingNumber = "_CanonicalRankingNumber";
inline constexpr std::string_view _ChiralityPossible = "_ChiralityPossible";
inline constexpr std::string_view _CrippenLogP = "_CrippenLogP";
inline constexpr std::string_view _CrippenMR = "_CrippenMR";
inline constexpr std::string_view _MMFFSanitized = "_MMFFSanitized";
inline constexpr std::string_view _MolFileChiralFlag = "_MolFileChiralFlag";
inline constexpr std::string_view MRV_SMA = "MRV SMA";
inline constexpr std::string_view _MolFileRLabel = "_MolFileRLabel";
inline constexpr std::string_view _MolFileAtomQuery = "_MolFileAtomQuery";
inline constexpr std::string_view _MolFileBondQuery = "_MolFileBondQuery";
inline constexpr std::string_view _MolFileBondEndPts = "_MolFileBondEndPts";
inline constexpr std::string_view _MolFileBondAttach = "_MolFileBondAttach";
inline constexpr std::string_view _MolFileBondType = "_MolFileBondType";
inline constexpr std::string_view _MolFileBondStereo = "_MolFileBondStereo";
inline constexpr std::string_view _MolFileBondCfg = "_MolFileBondCfg";

inline constexpr std::string_view _Name = "_Name";
inline constexpr std::string_view _NeedsQueryScan = "_NeedsQueryScan";
inline constexpr std::string_view _NonExplicit3DChirality = "_NonExplicit3DChirality";
inline constexpr std::string_view _QueryFormalCharge = "_QueryFormalCharge";
inline constexpr std::string_view _QueryHCount = "_QueryHCount";
inline constexpr std::string_view _QueryIsotope = "_QueryIsotope";
inline constexpr std::string_view _QueryMass = "_QueryMass";
inline constexpr std::string_view _ReactionDegreeChanged = "_ReactionDegreeChanged";
inline constexpr std::string_view reactantAtomIdx = "react_atom_idx";
inline constexpr std::string_view reactionMapNum = "old_mapno";
inline constexpr std::string_view reactantIdx = "react_idx";

inline constexpr std::string_view _RingClosures = "_RingClosures";
inline constexpr std::string_view _SLN_s = "_SLN_s";
inline constexpr std::string_view _SmilesStart = "_SmilesStart";
inline constexpr std::string_view _StereochemDone = "_StereochemDone";
inline constexpr std::string_view _TraversalBondIndexOrder = "_TraversalBondIndexOrder";
inline constexpr std::string_view _TraversalRingClosureBond =
    "_TraversalRingClosureBond";
inline constexpr std::string_view _TraversalStartPoint = "_TraversalStartPoint";
inline constexpr std::string_view _TriposAtomType = "_TriposAtomType";
inline constexpr std::string_view _Unfinished_SLN_ = "_Unfinished_SLN_";
inline constexpr std::string_view _UnknownStereo = "_UnknownStereo";
inline constexpr std::string_view _connectivityHKDeltas = "_connectivityHKDeltas";
inline constexpr std::string_view _connectivityNVals = "_connectivityNVals";
inline constexpr std::string_view _crippenLogP = "_crippenLogP";
inline constexpr std::string_view _crippenLogPContribs = "_crippenLogPContribs";
inline constexpr std::string_view _crippenMR = "_crippenMR";
inline constexpr std::string_view _crippenMRContribs = "_crippenMRContribs";
inline constexpr std::string_view _GasteigerCharge = "_GasteigerCharge";
inline constexpr std::string_view _GasteigerHCharge = "_GasteigerHCharge";
inline constexpr std::string_view _doIsoSmiles = "_doIsoSmiles";
inline constexpr std::string_view _fragSMARTS = "_fragSMARTS";
inline constexpr std::string_view _hasMassQuery = "_hasMassQuery";
inline constexpr std::string_view _labuteASA = "_labuteASA";
inline constexpr std::string_view _labuteAtomContribs = "_labuteAtomContribs";
inline constexpr std::string_view _labuteAtomHContrib = "_labuteAtomHContrib";
inline constexpr std::string_view _protected = "_protected";
inline constexpr std::string_view _queryRootAtom = "_queryRootAtom";
inline constexpr std::string_view _ringStereoAtoms = "_ringStereoAtoms";
inline constexpr std::string_view _ringStereoWarning = "_ringStereoWarning";
inline constexpr std::string_view _ringStereochemCand = "_ringStereochemCand";
inline constexpr std::string_view _ringStereoOtherAtom = "_ringStereoOtherAtom";
inline constexpr std::string_view _mesoOtherAtom = "_mesoOtherAtom";
inline constexpr std::string_view _chiralPermutation = "_chiralPermutation";
inline constexpr std::string_view _smilesAtomOutputOrder = "_smilesAtomOutputOrder";
inline constexpr std::string_view _smilesBondOutputOrder = "_smilesBondOutputOrder";
inline constexpr std::string_view _starred = "_starred";
inline constexpr std::string_view _supplementalSmilesLabel = "_supplementalSmilesLabel";
inline constexpr std::string_view _tpsa = "_tpsa";
inline constexpr std::string_view _tpsaAtomContribs = "_tpsaAtomContribs";
inline constexpr std::string_view _unspecifiedOrder = "_unspecifiedOrder";
inline constexpr std::string_view _brokenChirality = "_brokenChirality";
inline constexpr std::string_view _rgroupAtomMaps = "_rgroupAtomMaps";
inline constexpr std::string_view _rgroupBonds = "_rgroupBonds";
inline constexpr std::string_view _rgroupTargetAtoms = "_rgroupTargetAtoms";
inline constexpr std::string_view _rgroupTargetBonds = "_rgroupTargetBonds";
inline constexpr std::string_view dummyLabel = "dummyLabel";
inline constexpr std::string_view extraRings = "extraRings";
inline constexpr std::string_view isImplicit = "isImplicit";
inline constexpr std::string_view maxAttachIdx = "maxAttachIdx";
inline constexpr std::string_view molAtomMapNumber = "molAtomMapNumber";
inline constexpr std::string_view molFileAlias = "molFileAlias";
inline constexpr std::string_view molFileValue = "molFileValue";
inline constexpr std::string_view molInversionFlag = "molInversionFlag";
inline constexpr std::string_view molParity = "molParity";
inline constexpr std::string_view molStereoCare = "molStereoCare";
inline constexpr std::string_view molRxnComponent = "molRxnComponent";
inline constexpr std::string_view molRxnRole = "molRxnRole";
inline constexpr std::string_view molTotValence = "molTotValence";
inline constexpr std::string_view molFileLinkNodes = "_molLinkNodes";
inline constexpr std::string_view numArom = "numArom";
inline constexpr std::string_view ringMembership = "ringMembership";
inline constexpr std::string_view smilesSymbol = "smilesSymbol";
inline constexpr std::string_view atomLabel = "atomLabel";
inline constexpr std::string_view OxidationNumber = "OxidationNumber";
inline constexpr std::string_view internalRgroupSmiles = "internalRgroupSmiles";
inline constexpr std::string_view molRingBondCount = "molRingBondCount";
inline constexpr std::string_view molSubstCount = "molSubstCount";
inline constexpr std::string_view molAttachPoint = "molAttchpt";
inline constexpr std::string_view molAttachOrder = "molAttchord";
inline constexpr std::string_view molAttachOrderTemplate = "molAttachOrderTemplate";
inline constexpr std::string_view molAtomClass = "molClass";
inline constexpr std::string_view molAtomSeqId = "molSeqid";
inline constexpr std::string_view molAtomSeqName = "molSeqName";
inline constexpr std::string_view molRxnExactChange = "molRxnExachg";
inline constexpr std::string_view molReactStatus = "molReactStatus";
inline constexpr std::string_view _fromAttachPoint = "_fromAttchpt";
inline constexpr std::string_view natReplace = "natReplace";
inline constexpr std::string_view templateNames = "templateNames";

inline constexpr std::string_view molNote = "molNote";
inline constexpr std::string_view atomNote = "atomNote";
inline constexpr std::string_view bondNote = "bondNote";
inline constexpr std::string_view _isotopicHs = "_isotopicHs";

inline constexpr std::string_view _QueryAtomGenericLabel = "_QueryAtomGenericLabel";

// molecule drawing
inline constexpr std::string_view _displayLabel = "_displayLabel";
inline constexpr std::string_view _displayLabelW = "_displayLabelW";

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
