//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBSTRUCT_UTILS_H_
#define _RD_SUBSTRUCT_UTILS_H_

#include "SubstructMatch.h"
#include <boost/smart_ptr.hpp>

namespace RDKit{
  class ROMol;
  class Atom;
  class Bond;
  typedef boost::shared_ptr<Atom>    ATOM_SPTR;
  typedef boost::shared_ptr<Bond>    BOND_SPTR;
  
  double toPrime(const MatchVectType &v);
  void removeDuplicates(std::vector<MatchVectType> &v,unsigned int nAtoms);
  bool atomCompat(const ATOM_SPTR a1,const ATOM_SPTR a2);
  bool chiralAtomCompat(const ATOM_SPTR a1,const ATOM_SPTR a2);
  bool bondCompat(const BOND_SPTR b1,const BOND_SPTR b2);
}


#endif
