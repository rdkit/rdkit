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
#include "Exceptions.h"

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
///////////////////////////////////////////////////////////////
// Molecule Props
const int _Name                  = 0;   // string
const int MolFileInfo            = 1;   // string
const int MolFileComments        = 2;   // string
const int _2DConf                = 3;   // int (combine into dimension?)
const int _3DConf                = 4;   // int
const int _doIsoSmiles           = 5;   // int (should probably be removed)
const int extraRings             = 6;   // vec<vec<int> >
const int _smilesAtomOutputOrder = 7;   // vec<int> computed
const int _StereochemDone        = 8;   // int
const int _NeedsQueryScan        = 9;   // int (bool)
const int _fragSMARTS            = 10;  // std::string
const int maxAttachIdx           = 11;  // int TemplEnumTools.cpp
const int origNoImplicit         = 12;  // int (bool)
const int ringMembership         = 13;  //? unused (molopstest.cpp)

// Computed Values
// ConnectivityDescriptors
const int _connectivityHKDeltas = 14;   // std::vector<double> computed
const int _connectivityNVals    = 15;   // std::vector<double> computed
const int _crippenLogP          = 16;   // double computed
const int _crippenLogPContribs  = 17;   // std::vector<double> computed
const int _crippenMR            = 18;   // double computed
const int _crippenMRContribs    = 19;   // std::vector<double> computed
const int _labuteASA            = 20;   // double computed
const int _labuteAtomContribs   = 21;   // vec<double> computed
const int _labuteAtomHContrib   = 22;   // double computed
const int _tpsa                 = 23;   // double computed
const int _tpsaAtomContribs     = 24;   // vec<double> computed
const int numArom               = 25;   // int computed (only uses in tests?)
const int _MMFFSanitized        = 26;   // int (bool) computed
const int _CrippenLogP          = 27;   // Unused (in the basement)
const int _CrippenMR            = 28;   // Unused (in the basement)

///////////////////////////////////////////////////////////////
// Atom Props

// Chirality stuff
const int _BondsPotentialStereo = 29;   // int (or bool) COMPUTED
const int _CIPCode              = 30;   // std::string COMPUTED
const int _CIPRank              = 31;   // int COMPUTED
const int _ChiralityPossible    = 32;   // int
const int _UnknownStereo        = 33;   // int (bool) AddHs/Chirality
const int _ringStereoAtoms      = 34;   // int vect Canon/Chiral/MolHash/MolOps//Renumber
                                        //RWmol
const int _ringStereochemCand   = 35;   // chirality bool COMPUTED
const int _ringStereoWarning    = 36;   // obsolete ?

// Smiles parsing
const int _SmilesStart              = 37; // int
const int _TraversalBondIndexOrder  = 38; // ? unused
const int _TraversalRingClosureBond = 39; // unsigned int
const int _TraversalStartPoint      = 40; // bool
const int _queryRootAtom            = 41; // int SLNParse/SubstructMatch
const int _hasMassQuery             = 42; // atom bool
const int _protected                = 43; // atom int (bool)
const int _supplementalSmilesLabel  = 44; // atom string (SmilesWrite)
const int _unspecifiedOrder         = 45; // atom int (bool) smarts/smiles
const int _RingClosures             = 46; // INT_VECT smarts/smiles/canon

// MDL Style Properties (MolFileParser)
const int molAtomMapNumber   = 47;      // int
const int molFileAlias       = 48;      // string
const int molFileValue       = 49;      // string
const int molInversionFlag   = 50;      // int
const int molParity          = 51;      // int
const int molRxnComponent    = 52;      // int
const int molRxnRole         = 53;      // int
const int molTotValence      = 54;      // int
const int _MolFileRLabel     = 55;      // int
const int _MolFileChiralFlag = 56;      // int
const int dummyLabel         = 57;       // atom string
    
// Reaction Information (Reactions.cpp)
const int _QueryFormalCharge     = 58;  //  int
const int _QueryHCount           = 59;  // int
const int _QueryIsotope          = 60;  // int
const int _QueryMass             = 61;  // int = round(float * 1000)
const int _ReactionDegreeChanged = 62;  // int (bool)
const int NullBond               = 63;  // int (bool)
const int _rgroupAtomMaps        = 64;  // int
const int _rgroupBonds           = 65;  // int

// SLN
const int _AtomID          = 66;        // unsigned int SLNParser
const int _starred         = 67;        // atom int COMPUTED (SLN)
const int _SLN_s           = 68;        // string SLNAttribs (chiral info)
const int _Unfinished_SLN_ = 69;        // int (bool)

// Smarts Smiles
const int _brokenChirality = 70;        // atom bool
const int isImplicit       = 71;        // atom int (bool)
const int smilesSymbol     = 72;        // atom string (only used in test?)

// Tripos
const int _TriposAtomType = 73; // string Mol2FileParser
// missing defs for _TriposAtomName//_TriposPartialCharge...


///////////////////////////////////////////////////////////////
// misc props
const int TWOD                 = 74;    // need THREED -> confusing using in TDTMol supplier                                           //  converge with _2DConf?
const int BalabanJ             = 75;    // mol double
const int BalanbanJ            = 76;    // typo!! fix...
const int Discrims             = 77;    // FragCatalog Entry
                                        // Subgraphs::DiscrimTuple (uint32,uint32,uint32)
const int DistanceMatrix_Paths = 78;    // boost::shared_array<double>
                                        //  - note, confusing creation of names in
                                        //  - getDistanceMat
const int MAX                  = 78;   // reserved spaces
const int RESERVED             = 256;
///////////////////////////////////////////////////////////////
// Name Lookup (see Dict.cpp)
extern const char * propnames[];
const char *getPropName(int);

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
