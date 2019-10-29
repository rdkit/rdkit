//
//  Copyright (C) 2003-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SUBSTRUCT_UTILS_H
#define RD_SUBSTRUCT_UTILS_H

#include "SubstructMatch.h"

namespace RDKit {
class ROMol;
class Atom;
class Bond;

RDKIT_SUBSTRUCTMATCH_EXPORT double toPrime(const MatchVectType& v);
RDKIT_SUBSTRUCTMATCH_EXPORT void removeDuplicates(std::vector<MatchVectType>& v,
                                                  unsigned int nAtoms);
RDKIT_SUBSTRUCTMATCH_EXPORT bool atomCompat(const Atom* a1, const Atom* a2,
                                            const SubstructMatchParameters& ps);
RDKIT_SUBSTRUCTMATCH_EXPORT bool chiralAtomCompat(const Atom* a1,
                                                  const Atom* a2);
RDKIT_SUBSTRUCTMATCH_EXPORT bool bondCompat(const Bond* b1, const Bond* b2,
                                            const SubstructMatchParameters& ps);
}  // namespace RDKit

#endif
