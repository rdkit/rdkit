//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SUBSTRUCT_UTILS_H_
#define _RD_SUBSTRUCT_UTILS_H_

#include "SubstructMatch.h"
#include <boost/smart_ptr.hpp>

namespace RDKit {
class ROMol;
class Atom;
class Bond;

RDKIT_SUBSTRUCTMATCH_EXPORT double toPrime(const MatchVectType &v);
RDKIT_SUBSTRUCTMATCH_EXPORT void removeDuplicates(std::vector<MatchVectType> &v, unsigned int nAtoms);
RDKIT_SUBSTRUCTMATCH_EXPORT bool atomCompat(const Atom* a1, const Atom* a2,
                bool useQueryQueryMatches = false);
RDKIT_SUBSTRUCTMATCH_EXPORT bool chiralAtomCompat(const Atom* a1, const Atom* a2);
RDKIT_SUBSTRUCTMATCH_EXPORT bool bondCompat(const Bond* b1, const Bond* b2,
                bool useQueryQueryMatches = false);
}

#endif
