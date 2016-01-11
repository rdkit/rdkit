//
// Copyright (C) 2001-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef __RD_TYPES_H__
#define __RD_TYPES_H__

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include <cmath>

#include <RDGeneral/Invariant.h>
#include "Dict.h"

namespace detail {
// used in various places for computed properties
const std::string computedPropName = "__computedProps";
}

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

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

namespace RDKit {
namespace common_properties {
extern const std::string TWOD;
extern const std::string BalabanJ;
extern const std::string BalanbanJ;
extern const std::string Discrims;
extern const std::string DistanceMatrix_Paths;
extern const std::string MolFileComments;
extern const std::string MolFileInfo;
extern const std::string NullBond;
extern const std::string _2DConf;
extern const std::string _3DConf;
extern const std::string _AtomID;
extern const std::string _BondsPotentialStereo;
extern const std::string _CIPCode;
extern const std::string _CIPRank;
extern const std::string _ChiralityPossible;
extern const std::string _CrippenLogP;
extern const std::string _CrippenMR;
extern const std::string _MMFFSanitized;
extern const std::string _MolFileChiralFlag;
extern const std::string _MolFileRLabel;
extern const std::string _Name;
extern const std::string _NeedsQueryScan;
extern const std::string _QueryFormalCharge;
extern const std::string _QueryHCount;
extern const std::string _QueryIsotope;
extern const std::string _QueryMass;
extern const std::string _ReactionDegreeChanged;
extern const std::string _RingClosures;
extern const std::string _SLN_s;
extern const std::string _SmilesStart;
extern const std::string _StereochemDone;
extern const std::string _TraversalBondIndexOrder;
extern const std::string _TraversalRingClosureBond;
extern const std::string _TriposAtomType;
extern const std::string _Unfinished_SLN_;
extern const std::string _UnknownStereo;
extern const std::string _connectivityHKDeltas;
extern const std::string _connectivityNVals;
extern const std::string _crippenLogP;
extern const std::string _crippenLogPContribs;
extern const std::string _crippenMR;
extern const std::string _crippenMRContribs;
extern const std::string _doIsoSmiles;
extern const std::string _fragSMARTS;
extern const std::string _hasMassQuery;
extern const std::string _labuteASA;
extern const std::string _labuteAtomContribs;
extern const std::string _labuteAtomHContrib;
extern const std::string _protected;
extern const std::string _queryRootAtom;
extern const std::string _ringStereoAtoms;
extern const std::string _ringStereoWarning;
extern const std::string _ringStereochemCand;
extern const std::string _smilesAtomOutputOrder;
extern const std::string _starred;
extern const std::string _supplementalSmilesLabel;
extern const std::string _tpsa;
extern const std::string _tpsaAtomContribs;
extern const std::string _unspecifiedOrder;
extern const std::string _brokenChirality;
extern const std::string _rgroupAtomMaps;
extern const std::string _rgroupBonds;
extern const std::string dummyLabel;
extern const std::string extraRings;
extern const std::string isImplicit;
extern const std::string maxAttachIdx;
extern const std::string molAtomMapNumber;
extern const std::string molFileAlias;
extern const std::string molFileValue;
extern const std::string molInversionFlag;
extern const std::string molParity;
extern const std::string molRxnComponent;
extern const std::string molRxnRole;
extern const std::string molTotValence;
extern const std::string numArom;
extern const std::string origNoImplicit;
extern const std::string ringMembership;
extern const std::string smilesSymbol;
}  // end common_properties
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

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double EPS_DOUBLE = std::numeric_limits<double>::epsilon();
const double SMALL_DOUBLE = 1.0e-8;
const double MAX_INT = static_cast<double>(std::numeric_limits<int>::max());
const double MAX_LONGINT =
    static_cast<double>(std::numeric_limits<LONGINT>::max());

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
struct ltDouble {
 public:
  ltDouble() : _tol(1.0e-8){};
  bool operator()(double d1, double d2) const {
    if (fabs(d1 - d2) < _tol) {
      return false;
    } else {
      return (d1 < d2);
    }
  }

 private:
  double _tol;
};

//! std::map from double to integer.
typedef std::map<double, int, ltDouble> DOUBLE_INT_MAP;

//! functor for returning the larger of two values
template <typename T>
struct larger_of {
  T operator()(T arg1, T arg2) { return arg1 > arg2 ? arg1 : arg2; };
};

//! functor for comparing two strings
struct charptr_functor {
  bool operator()(const char *s1, const char *s2) const {
    // std::cout << s1 << " " << s2 << " " << strcmp(s1, s2) << "\n";

    return strcmp(s1, s2) < 0;
  };
};

//! \brief calculate the union of two INT_VECTs and put the results in a
//! third vector
void Union(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res);

//! \brief calculate the intersection of two INT_VECTs and put the results in a
//! third vector
void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res);

//! calculating the union of the INT_VECT's in a VECT_INT_VECT
/*!
    \param rings   the INT_VECT's to consider
    \param res     used to return results
    \param exclude any values in this optional INT_VECT will be excluded
           from the union.
*/
void Union(const VECT_INT_VECT &rings, INT_VECT &res,
           const INT_VECT *exclude = NULL);

//! given a current combination of numbers change it to the next possible
//combination
/*!
  \param comb the <b>sorted</b> vector to consider
  \param tot the maximum number possible in the vector

  \return -1 on failure, the index of the last number changed on success.
  Example:
    for all combinations 3 of numbers between 0 and tot=5
    given (0,1,2) the function wil return (0,1,3) etc.


*/
int nextCombination(INT_VECT &comb, int tot);

//! rounds a value to the closest int
double round(double v);

};  // end of namespace

#endif
