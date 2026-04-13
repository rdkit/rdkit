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
#include <memory>
#include <numeric>
#include <list>
#include <limits>

#include <cstring>
#include <any>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

class PropToken {
  struct PropTokenImpl {
    uint64_t hash;
    std::string text;
    PropTokenImpl(uint64_t initHash, std::string initText)
        : hash(initHash), text(std::move(initText)) {}
    PropTokenImpl(uint64_t initHash, const std::string_view &initText)
        : hash(initHash), text(initText) {}
  };
  std::shared_ptr<PropTokenImpl> impl;

 public:
  static uint64_t computeHash(const char *text) {
    size_t length = strlen(text);
    uint64_t hash = 0;
    while (length >= 8) {
      uint64_t v = *reinterpret_cast<const uint64_t*>(text);
      // The large magic number is the golden ratio minus 1, times 2^64,
      // which is good for a simple low discrepancy sequence.
      // 12345 is arbitrary.
      hash = 12345 * hash + v + 0x9E3779B97F4A7C15;
      length -= 8;
      text += 8;
    }
    if (length != 0) {
      uint64_t v = 0;
      for (size_t i = 0; i < length; ++i) {
        v |= uint64_t(uint8_t(*text)) << (i * 8);
        ++text;
      }
      hash = 12345 * hash + v + 0x9E3779B97F4A7C15;
    }
    return hash;
  }

  constexpr PropToken() = default;
  explicit PropToken(const char *text)
      : impl(std::make_shared<PropTokenImpl>(computeHash(text),
                                             std::string(text))) {}
  explicit PropToken(const std::string &text)
      : impl(std::make_shared<PropTokenImpl>(computeHash(text.c_str()), text)) {}
  explicit PropToken(const std::string_view &text)
      : impl(nullptr) {
    std::string textString(text);
    auto hash = computeHash(textString.c_str());
    impl.reset(new PropTokenImpl(hash, std::move(textString)));
  }
  PropToken(const PropToken &) = default;
  PropToken(PropToken &&) = default;
  PropToken &operator=(const PropToken &) = default;
  PropToken &operator=(PropToken &&) = default;
  ~PropToken() = default;
  bool operator==(const PropToken &other)
          const {
    const PropTokenImpl *p0 = impl.get();
    const PropTokenImpl *p1 = other.impl.get();
    if (p0 == p1) {
      return true;
    }
    if (p0->hash != p1->hash) {
      return false;
    }
    // For common tokens, this case should be relatively uncommon
    return (p0->text == p1->text);
  }
  bool operator!=(const PropToken &other) const {
    return !(*this == other);
  }
  bool operator==(const std::string &other) const {
    return (impl->text == other);
  }
  bool operator!=(const std::string &other) const {
    return !(*this == other);
  }
  bool operator==(const char *other) const {
    return (impl->text == other);
  }
  bool operator!=(const char *other) const {
    return !(*this == other);
  }
  const std::string &getString() const {
    return impl->text;
  }
  const char *getText() const {
    return impl->text.c_str();
  }
  bool isNull() const { return impl.get() == nullptr; }
  operator bool() const { return impl.get() != nullptr; }
  bool operator!() const { return impl.get() == nullptr; }
};

namespace detail {
// used in various places for computed properties
inline constexpr std::string_view computedPropName = "__computedProps";
RDKIT_RDGENERAL_EXPORT extern const PropToken computedPropNameToken;

// This wraps a vector or local buffer as a bit set,
// so that often, no allocation is needed, and if it is
// needed, it can be reused.
class BitSetWrapper {
  constexpr static size_t localBufferSize = 2;
  uint64_t localBuffer[2] = {0, 0};
  uint64_t *const buffer;
  const size_t n;

 public:
  BitSetWrapper(std::vector<uint64_t> &storage, size_t n_, bool value = false)
      : buffer((n_ > 64 * localBufferSize)
                   ? (storage.clear(), storage.resize((n_ + 63) / 64, value ? ~uint64_t(0) : 0),
                      storage.data())
                   : localBuffer),
        n(n_) {
    if (value && n_ <= 64 * localBufferSize) {
      for (size_t i = 0; i < localBufferSize; ++i) {
        localBuffer[i] = ~uint64_t(0);
      }
    }
  }
  BitSetWrapper(uint64_t *storage, size_t n_) : buffer(storage), n(n_) {}
  BitSetWrapper(uint64_t *storage, size_t n_, bool value)
      : buffer(storage), n(n_) {
    const size_t bufferSize = (n + 63) / 64;
    for (size_t i = 0; i < bufferSize; ++i) {
      buffer[i] = value ? ~uint64_t(0) : 0;
    }
  }

  bool get(size_t i) const {
    return (buffer[i / 64] & (uint64_t(1) << (i % 64))) != 0;
  }
  bool operator[](size_t i) const { return get(i); }
  void set(size_t i) {
    PRECONDITION(i < n, "BitSetWrapper::set index out of bounds");
    buffer[i / 64] |= (uint64_t(1) << (i % 64));
  }
  void reset(size_t i) {
    PRECONDITION(i < n, "BitSetWrapper::reset index out of bounds");
    buffer[i / 64] &= ~(uint64_t(1) << (i % 64));
  }
  void set(size_t i, bool value) {
    PRECONDITION(i < n, "BitSetWrapper::set index out of bounds");
    if (value) {
      buffer[i / 64] |= (uint64_t(1) << (i % 64));
    } else {
      buffer[i / 64] &= ~(uint64_t(1) << (i % 64));
    }
  }
  void set() {
    for (size_t i = 0; 64 * i < n; ++i) {
      buffer[i] = ~uint64_t(0);
    }
  }
  void reset() {
    for (size_t i = 0; 64 * i < n; ++i) {
      buffer[i] = 0;
    }
  }
  size_t size() const { return n; }

  const uint64_t *data() const { return buffer; }
  uint64_t *data() { return buffer; }
  size_t dataSize() const { return (n + 63) / 64; }

  bool empty() const {
    assert(n != 0);
    if (n == 0) {
      return true;
    }
    const size_t bufferSize = (n + 63) / 64;
    // Just OR together all the bits
    uint64_t bits = 0;
    for (size_t i = 0; i < bufferSize-1; ++i) {
      bits |= buffer[i];
    }
    // Only include the valid bits for the last buffer element.
    bits |= (buffer[bufferSize - 1] & (~uint64_t(0) >> ((-int64_t(n)) & 0x3F)));
    return (bits == 0);
  }
};
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
inline constexpr std::string_view _BondsPotentialStereo =
    "_BondsPotentialStereo";
inline constexpr std::string_view _ChiralAtomRank = "_chiralAtomRank";
inline constexpr std::string_view _CIPCode = "_CIPCode";
inline constexpr std::string_view _CIPRank = "_CIPRank";
inline constexpr std::string_view _CIPComputed = "_CIPComputed";
inline constexpr std::string_view _CIPNeighborOrder = "_CIPNeighborOrder";
inline constexpr std::string_view _CanonicalRankingNumber =
    "_CanonicalRankingNumber";
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
inline constexpr std::string_view _NonExplicit3DChirality =
    "_NonExplicit3DChirality";
inline constexpr std::string_view _QueryFormalCharge = "_QueryFormalCharge";
inline constexpr std::string_view _QueryHCount = "_QueryHCount";
inline constexpr std::string_view _QueryIsotope = "_QueryIsotope";
inline constexpr std::string_view _QueryMass = "_QueryMass";
inline constexpr std::string_view _ReactionDegreeChanged =
    "_ReactionDegreeChanged";
inline constexpr std::string_view reactantAtomIdx = "react_atom_idx";
inline constexpr std::string_view reactionMapNum = "old_mapno";
inline constexpr std::string_view reactantIdx = "react_idx";

inline constexpr std::string_view _RingClosures = "_RingClosures";
inline constexpr std::string_view _SLN_s = "_SLN_s";
inline constexpr std::string_view _SmilesStart = "_SmilesStart";
inline constexpr std::string_view _StereochemDone = "_StereochemDone";
inline constexpr std::string_view _TraversalBondIndexOrder =
    "_TraversalBondIndexOrder";
inline constexpr std::string_view _TraversalRingClosureBond =
    "_TraversalRingClosureBond";
inline constexpr std::string_view _TraversalStartPoint = "_TraversalStartPoint";
inline constexpr std::string_view _TriposAtomType = "_TriposAtomType";
inline constexpr std::string_view _Unfinished_SLN_ = "_Unfinished_SLN_";
inline constexpr std::string_view _UnknownStereo = "_UnknownStereo";
inline constexpr std::string_view _connectivityHKDeltas =
    "_connectivityHKDeltas";
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
inline constexpr std::string_view _smilesAtomOutputOrder =
    "_smilesAtomOutputOrder";
inline constexpr std::string_view _smilesBondOutputOrder =
    "_smilesBondOutputOrder";
inline constexpr std::string_view _starred = "_starred";
inline constexpr std::string_view _supplementalSmilesLabel =
    "_supplementalSmilesLabel";
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
inline constexpr std::string_view molAttachOrderTemplate =
    "molAttachOrderTemplate";
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

inline constexpr std::string_view _QueryAtomGenericLabel =
    "_QueryAtomGenericLabel";

// molecule drawing
inline constexpr std::string_view _displayLabel = "_displayLabel";
inline constexpr std::string_view _displayLabelW = "_displayLabelW";

///////////////////////////////////////////////////////////////
// misc props
RDKIT_RDGENERAL_EXPORT extern const PropToken _hasMassQueryToken;  // atom bool

RDKIT_RDGENERAL_EXPORT extern const PropToken _ChiralityPossibleToken;  // bool
RDKIT_RDGENERAL_EXPORT extern const PropToken _chiralPermutationToken;  // uint
RDKIT_RDGENERAL_EXPORT extern const PropToken _CIPCodeToken;            // char
RDKIT_RDGENERAL_EXPORT extern const PropToken _CIPRankToken;            // uint
RDKIT_RDGENERAL_EXPORT extern const PropToken _isotopicHsToken;         // uint64
RDKIT_RDGENERAL_EXPORT extern const PropToken _MolFileBondEndPtsToken;  // string
RDKIT_RDGENERAL_EXPORT extern const PropToken _MolFileBondAttachToken;  // string
RDKIT_RDGENERAL_EXPORT extern const PropToken _MolFileRLabelToken;      // uint
RDKIT_RDGENERAL_EXPORT extern const PropToken _ringStereoAtomsAllToken; // int
RDKIT_RDGENERAL_EXPORT extern const PropToken _ringStereoAtomsBeginsToken;// uint
RDKIT_RDGENERAL_EXPORT extern const PropToken _ringStereoGroupToken;    // int
RDKIT_RDGENERAL_EXPORT extern const PropToken _supplementalSmilesLabelToken;// string
RDKIT_RDGENERAL_EXPORT extern const PropToken _UnknownStereoToken;      // bool
RDKIT_RDGENERAL_EXPORT extern const PropToken dummyLabelToken;          // token
RDKIT_RDGENERAL_EXPORT extern const PropToken isImplicitToken;          // bool
RDKIT_RDGENERAL_EXPORT extern const PropToken molAtomMapNumberToken;    // int
RDKIT_RDGENERAL_EXPORT extern const PropToken molFileAliasToken;        // string
RDKIT_RDGENERAL_EXPORT extern const PropToken molFileValueToken;        // string

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
