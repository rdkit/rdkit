//
// Copyright (C) 2001-2016 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_TYPES_H
#define RD_TYPES_H

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include <cmath>

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
#include <RDGeneral/BoostStartInclude.h>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

  namespace detail {
  // used in various places for computed properties
  extern const std::string computedPropName;
  }

namespace common_properties {
///////////////////////////////////////////////////////////////
// Molecule Props
extern const std::string _Name;           // string
extern const std::string MolFileInfo;     // string
extern const std::string MolFileComments; // string
extern const std::string _2DConf;         // int (combine into dimension?)
extern const std::string _3DConf;         // int
extern const std::string _doIsoSmiles;    // int (should probably be removed)
extern const std::string extraRings;      // vec<vec<int> >
extern const std::string _smilesAtomOutputOrder; // vec<int> computed
extern const std::string _StereochemDone; // int
extern const std::string _NeedsQueryScan; // int (bool)
extern const std::string _fragSMARTS;     // std::string
extern const std::string maxAttachIdx;    // int TemplEnumTools.cpp
extern const std::string origNoImplicit;  // int (bool)
extern const std::string ringMembership;  //? unused (molopstest.cpp)

// Computed Values
// ConnectivityDescriptors
extern const std::string _connectivityHKDeltas;// std::vector<double> computed
extern const std::string _connectivityNVals;   // std::vector<double> computed

extern const std::string _crippenLogP;         // double computed
extern const std::string _crippenLogPContribs; // std::vector<double> computed

extern const std::string _crippenMR;           // double computed
extern const std::string _crippenMRContribs;   // std::vector<double> computed

extern const std::string _labuteASA;           // double computed
extern const std::string _labuteAtomContribs;  // vec<double> computed
extern const std::string _labuteAtomHContrib;  // double computed

extern const std::string _tpsa;                // double computed
extern const std::string _tpsaAtomContribs;    // vec<double> computed

extern const std::string numArom;              // int computed (only uses in tests?)
extern const std::string _MMFFSanitized;       // int (bool) computed

extern const std::string _CrippenLogP; // Unused (in the basement)
extern const std::string _CrippenMR;   // Unused (in the basement)

///////////////////////////////////////////////////////////////
// Atom Props

// Chirality stuff
extern const std::string _BondsPotentialStereo; // int (or bool) COMPUTED
extern const std::string _CIPCode; // std::string COMPUTED
extern const std::string _CIPRank; // int COMPUTED
extern const std::string _ChiralityPossible; // int
extern const std::string _UnknownStereo; // int (bool) AddHs/Chirality
extern const std::string _ringStereoAtoms; // int vect Canon/Chiral/MolHash/MolOps//Renumber//RWmol
extern const std::string _ringStereochemCand; // chirality bool COMPUTED
extern const std::string _ringStereoWarning; // obsolete ?

// Smiles parsing
extern const std::string _SmilesStart; // int
extern const std::string _TraversalBondIndexOrder; // ? unused
extern const std::string _TraversalRingClosureBond; // unsigned int
extern const std::string _TraversalStartPoint; // bool
extern const std::string _queryRootAtom; // int SLNParse/SubstructMatch
extern const std::string _hasMassQuery; // atom bool
extern const std::string _protected; // atom int (bool)
extern const std::string _supplementalSmilesLabel; // atom string (SmilesWrite)
extern const std::string _unspecifiedOrder;// atom int (bool) smarts/smiles
extern const std::string _RingClosures; // INT_VECT smarts/smiles/canon
extern const std::string atomLabel; // atom string from CXSMILES

// MDL Style Properties (MolFileParser)
extern const std::string molAtomMapNumber; // int
extern const std::string molFileAlias;  // string
extern const std::string molFileValue;  // string
extern const std::string molInversionFlag; // int
extern const std::string molParity;     // int
extern const std::string molRxnComponent; // int
extern const std::string molRxnRole;    // int
extern const std::string molTotValence; // int
extern const std::string _MolFileRLabel; // int
extern const std::string _MolFileChiralFlag; // int

extern const std::string dummyLabel; // atom string

// Reaction Information (Reactions.cpp)
extern const std::string _QueryFormalCharge; //  int
extern const std::string _QueryHCount; // int
extern const std::string _QueryIsotope; // int
extern const std::string _QueryMass; // int = round(float * 1000)
extern const std::string _ReactionDegreeChanged; // int (bool)
extern const std::string NullBond; // int (bool)
extern const std::string _rgroupAtomMaps;
extern const std::string _rgroupBonds;

// SLN
extern const std::string _AtomID; // unsigned int SLNParser
extern const std::string _starred; // atom int COMPUTED (SLN)
extern const std::string _SLN_s; // string SLNAttribs (chiral info)
extern const std::string _Unfinished_SLN_; // int (bool)

// Smarts Smiles
extern const std::string _brokenChirality; // atom bool
extern const std::string isImplicit;  // atom int (bool)
extern const std::string smilesSymbol; // atom string (only used in test?)

// Tripos
extern const std::string _TriposAtomType; // string Mol2FileParser
// missing defs for _TriposAtomName//_TriposPartialCharge...


///////////////////////////////////////////////////////////////
// misc props
extern const std::string TWOD; // need THREED -> confusing using in TDTMol supplier
                               //  converge with _2DConf?
extern const std::string BalabanJ; // mol double
extern const std::string BalanbanJ; // typo!! fix...

extern const std::string Discrims; // FragCatalog Entry
                                   // Subgraphs::DiscrimTuple (uint32,uint32,uint32)
extern const std::string DistanceMatrix_Paths; // boost::shared_array<double>
                                               //  - note, confusing creation of names in
                                               //  - getDistanceMat

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

extern const double MAX_DOUBLE;
extern const double EPS_DOUBLE;
extern const double SMALL_DOUBLE;
extern const double MAX_INT;
extern const double MAX_LONGINT;

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
// combination
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
