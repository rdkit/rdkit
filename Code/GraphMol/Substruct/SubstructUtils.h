//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
  bool atomCompat(const ATOM_SPTR &a1,const ATOM_SPTR &a2,bool useQueryQueryMatches=false);
  bool chiralAtomCompat(const ATOM_SPTR &a1,const ATOM_SPTR &a2);
  bool bondCompat(const BOND_SPTR &b1,const BOND_SPTR &b2,bool useQueryQueryMatches=false);
}


#endif
