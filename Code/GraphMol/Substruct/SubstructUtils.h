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

namespace RDKit {
class ROMol;
class Atom;
class Bond;

double toPrime(const MatchVectType &v);
void removeDuplicates(std::vector<MatchVectType> &v, unsigned int nAtoms);
bool atomCompat(const Atom* a1, const Atom* a2,
                bool useQueryQueryMatches = false);
bool chiralAtomCompat(const Atom* a1, const Atom* a2);
bool bondCompat(const Bond* b1, const Bond* b2,
                bool useQueryQueryMatches = false);
}

#endif
