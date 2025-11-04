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
RDKIT_RDGENERAL_EXPORT extern const std::string computedPropName;
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
RDKIT_RDGENERAL_EXPORT extern const std::string _Name;            // string
RDKIT_RDGENERAL_EXPORT extern const std::string MolFileInfo;      // string
RDKIT_RDGENERAL_EXPORT extern const std::string MolFileComments;  // string
RDKIT_RDGENERAL_EXPORT extern const std::string
    _2DConf;  // int (combine into dimension?)
RDKIT_RDGENERAL_EXPORT extern const std::string _3DConf;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _doIsoSmiles;  // int (should probably be removed)
RDKIT_RDGENERAL_EXPORT extern const std::string extraRings;  // vec<vec<int> >
RDKIT_RDGENERAL_EXPORT extern const std::string
    _smilesAtomOutputOrder;  // vec<int> computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _smilesBondOutputOrder;  // vec<int> computed
RDKIT_RDGENERAL_EXPORT extern const std::string _StereochemDone;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string _NeedsQueryScan;  // int (bool)
RDKIT_RDGENERAL_EXPORT extern const std::string _fragSMARTS;      // std::string
RDKIT_RDGENERAL_EXPORT extern const std::string
    maxAttachIdx;  // int TemplEnumTools.cpp
RDKIT_RDGENERAL_EXPORT extern const std::string
    ringMembership;  //? unused (molopstest.cpp)

// Computed Values
// ConnectivityDescriptors
RDKIT_RDGENERAL_EXPORT extern const std::string
    _connectivityHKDeltas;  // std::vector<double> computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _connectivityNVals;  // std::vector<double> computed

RDKIT_RDGENERAL_EXPORT extern const std::string
    _crippenLogP;  // double computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _crippenLogPContribs;  // std::vector<double> computed

RDKIT_RDGENERAL_EXPORT extern const std::string _crippenMR;  // double computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _crippenMRContribs;  // std::vector<double> computed

RDKIT_RDGENERAL_EXPORT extern const std::string _labuteASA;  // double computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _labuteAtomContribs;  // vec<double> computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _labuteAtomHContrib;  // double computed

RDKIT_RDGENERAL_EXPORT extern const std::string _tpsa;  // double computed
RDKIT_RDGENERAL_EXPORT extern const std::string
    _tpsaAtomContribs;  // vec<double> computed

RDKIT_RDGENERAL_EXPORT extern const std::string
    numArom;  // int computed (only uses in tests?)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _MMFFSanitized;  // int (bool) computed

RDKIT_RDGENERAL_EXPORT extern const std::string
    _CrippenLogP;  // Unused (in the basement)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _CrippenMR;  // Unused (in the basement)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _GasteigerCharge;  // used to hold partial charges
RDKIT_RDGENERAL_EXPORT extern const std::string
    _GasteigerHCharge;  // used to hold partial charges from implicit Hs

///////////////////////////////////////////////////////////////
// Atom Props

// Chirality stuff
RDKIT_RDGENERAL_EXPORT extern const std::string
    _BondsPotentialStereo;  // int (or bool) COMPUTED
RDKIT_RDGENERAL_EXPORT extern const std::string
    _CIPCode;  // std::string COMPUTED
RDKIT_RDGENERAL_EXPORT extern const std::string _CIPRank;  // int COMPUTED
RDKIT_RDGENERAL_EXPORT extern const std::string
    _CIPComputed;  // int (bool) COMPUTED
RDKIT_RDGENERAL_EXPORT extern const std::string
    _CanonicalRankingNumber;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string _ChiralityPossible;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _UnknownStereo;  // int (bool) AddHs/Chirality
RDKIT_RDGENERAL_EXPORT extern const std::string
    _ringStereoAtoms;  // int vect Canon/Chiral/MolHash/MolOps//Renumber//RWmol
RDKIT_RDGENERAL_EXPORT extern const std::string
    _ringStereochemCand;  // chirality bool COMPUTED
RDKIT_RDGENERAL_EXPORT extern const std::string
    _ringStereoWarning;  // obsolete ?
RDKIT_RDGENERAL_EXPORT extern const std::string _chiralPermutation;    // int
RDKIT_RDGENERAL_EXPORT extern const std::string _ringStereoOtherAtom;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string _mesoOtherAtom;        // int

// Smiles parsing
RDKIT_RDGENERAL_EXPORT extern const std::string _SmilesStart;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _TraversalBondIndexOrder;  // ? unused
RDKIT_RDGENERAL_EXPORT extern const std::string
    _TraversalRingClosureBond;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string _TraversalStartPoint;  // bool
RDKIT_RDGENERAL_EXPORT extern const std::string
    _queryRootAtom;  // int SLNParse/SubstructMatch
RDKIT_RDGENERAL_EXPORT extern const std::string _hasMassQuery;  // atom bool
RDKIT_RDGENERAL_EXPORT extern const std::string _protected;  // atom int (bool)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _ChiralAtomRank;  // atom rank (unsigned int)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _supplementalSmilesLabel;  // atom string (SmilesWrite)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _unspecifiedOrder;  // atom int (bool) smarts/smiles
RDKIT_RDGENERAL_EXPORT extern const std::string
    _RingClosures;  // INT_VECT smarts/smiles/canon
RDKIT_RDGENERAL_EXPORT extern const std::string
    atomLabel;  // atom string from CXSMILES
RDKIT_RDGENERAL_EXPORT extern const std::string OxidationNumber;  // int

// MDL Style Properties (MolFileParser)
RDKIT_RDGENERAL_EXPORT extern const std::string molAtomMapNumber;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string molFileAlias;      // string
RDKIT_RDGENERAL_EXPORT extern const std::string molFileValue;      // string
RDKIT_RDGENERAL_EXPORT extern const std::string molInversionFlag;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string molParity;         // int
RDKIT_RDGENERAL_EXPORT extern const std::string molStereoCare;     // int
RDKIT_RDGENERAL_EXPORT extern const std::string molRxnComponent;   // int
RDKIT_RDGENERAL_EXPORT extern const std::string molRxnRole;        // int
RDKIT_RDGENERAL_EXPORT extern const std::string molTotValence;     // int
RDKIT_RDGENERAL_EXPORT extern const std::string molRingBondCount;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string molSubstCount;     // int
RDKIT_RDGENERAL_EXPORT extern const std::string molAttachPoint;    // int
RDKIT_RDGENERAL_EXPORT extern const std::string molAttachOrder;    // int
RDKIT_RDGENERAL_EXPORT extern const std::string
    molAttachOrderTemplate;  // std::vector<AtomAttchOrd>

RDKIT_RDGENERAL_EXPORT extern const std::string molAtomClass;  // string
RDKIT_RDGENERAL_EXPORT extern const std::string natReplace;    // string
RDKIT_RDGENERAL_EXPORT extern const std::string
    templateNames;  // vector of strings
RDKIT_RDGENERAL_EXPORT extern const std::string molAtomSeqId;       // int
RDKIT_RDGENERAL_EXPORT extern const std::string molRxnExactChange;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string molReactStatus;     // int
RDKIT_RDGENERAL_EXPORT extern const std::string molFileLinkNodes;   // string
RDKIT_RDGENERAL_EXPORT extern const std::string _fromAttachPoint;   // int

RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileRLabel;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileChiralFlag;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileAtomQuery;   // int
RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileBondQuery;   // int
RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileBondEndPts;  // string
RDKIT_RDGENERAL_EXPORT extern const std::string _MolFileBondAttach;  // string
RDKIT_RDGENERAL_EXPORT extern const std::string
    _MolFileBondType;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _MolFileBondStereo;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _MolFileBondCfg;  // unsigned int
RDKIT_RDGENERAL_EXPORT extern const std::string
    MRV_SMA;  // smarts string from Marvin

// flag indicating that the chirality wasn't specified in the input,
// but was calculated from 3D coordinates in the input
RDKIT_RDGENERAL_EXPORT extern const std::string _NonExplicit3DChirality;  // int
RDKIT_RDGENERAL_EXPORT extern const std::string dummyLabel;  // atom string

RDKIT_RDGENERAL_EXPORT extern const std::string
    _QueryAtomGenericLabel;  // string

// Reaction Information (Reactions.cpp)
RDKIT_RDGENERAL_EXPORT extern const std::string _QueryFormalCharge;  //  int
RDKIT_RDGENERAL_EXPORT extern const std::string _QueryHCount;        // int
RDKIT_RDGENERAL_EXPORT extern const std::string _QueryIsotope;       // int
RDKIT_RDGENERAL_EXPORT extern const std::string
    _QueryMass;  // int = round(float * 1000)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _ReactionDegreeChanged;                                // int (bool)
RDKIT_RDGENERAL_EXPORT extern const std::string NullBond;  // int (bool)
RDKIT_RDGENERAL_EXPORT extern const std::string _rgroupAtomMaps;
RDKIT_RDGENERAL_EXPORT extern const std::string _rgroupBonds;
RDKIT_RDGENERAL_EXPORT extern const std::string _rgroupTargetAtoms;
RDKIT_RDGENERAL_EXPORT extern const std::string _rgroupTargetBonds;
RDKIT_RDGENERAL_EXPORT extern const std::string reactantAtomIdx;
RDKIT_RDGENERAL_EXPORT extern const std::string reactionMapNum;
RDKIT_RDGENERAL_EXPORT extern const std::string reactantIdx;

// SLN
RDKIT_RDGENERAL_EXPORT extern const std::string
    _AtomID;  // unsigned int SLNParser
RDKIT_RDGENERAL_EXPORT extern const std::string
    _starred;  // atom int COMPUTED (SLN)
RDKIT_RDGENERAL_EXPORT extern const std::string
    _SLN_s;  // string SLNAttribs (chiral info)
RDKIT_RDGENERAL_EXPORT extern const std::string _Unfinished_SLN_;  // int (bool)

// Smarts Smiles
RDKIT_RDGENERAL_EXPORT extern const std::string _brokenChirality;  // atom bool
RDKIT_RDGENERAL_EXPORT extern const std::string isImplicit;  // atom int (bool)
RDKIT_RDGENERAL_EXPORT extern const std::string
    smilesSymbol;  // atom string (only used in test?)

// Tripos
RDKIT_RDGENERAL_EXPORT extern const std::string
    _TriposAtomType;  // string Mol2FileParser
// missing defs for _TriposAtomName//_TriposPartialCharge...

// molecule drawing
RDKIT_RDGENERAL_EXPORT extern const std::string _displayLabel;   // string
RDKIT_RDGENERAL_EXPORT extern const std::string _displayLabelW;  // string

///////////////////////////////////////////////////////////////
// misc props
RDKIT_RDGENERAL_EXPORT extern const std::string
    TWOD;  // need THREED -> confusing using in TDTMol supplier
           //  converge with _2DConf?
RDKIT_RDGENERAL_EXPORT extern const std::string BalabanJ;   // mol double
RDKIT_RDGENERAL_EXPORT extern const std::string BalanbanJ;  // typo!! fix...

RDKIT_RDGENERAL_EXPORT extern const std::string Discrims;  // FragCatalog Entry
// Subgraphs::DiscrimTuple (uint32,uint32,uint32)
RDKIT_RDGENERAL_EXPORT extern const std::string
    DistanceMatrix_Paths;  // boost::shared_array<double>
//  - note, confusing creation of names in
//  - getDistanceMat
RDKIT_RDGENERAL_EXPORT extern const std::string internalRgroupSmiles;
RDKIT_RDGENERAL_EXPORT extern const std::string molNote;
RDKIT_RDGENERAL_EXPORT extern const std::string atomNote;
RDKIT_RDGENERAL_EXPORT extern const std::string bondNote;
RDKIT_RDGENERAL_EXPORT extern const std::string _isotopicHs;

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
