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
class RDProps;
class Atom;
class Bond;

RDKIT_SUBSTRUCTMATCH_EXPORT double toPrime(const MatchVectType& v);
RDKIT_SUBSTRUCTMATCH_EXPORT void removeDuplicates(std::vector<MatchVectType>& v,
                                                  unsigned int nAtoms);
RDKIT_SUBSTRUCTMATCH_EXPORT bool propertyCompat(const RDProps* r1, const RDProps* r2,
                                                const std::vector<std::string>& properties);
RDKIT_SUBSTRUCTMATCH_EXPORT bool atomCompat(const Atom* a1, const Atom* a2,
                                            const SubstructMatchParameters& ps);
RDKIT_SUBSTRUCTMATCH_EXPORT bool chiralAtomCompat(const Atom* a1,
                                                  const Atom* a2);
RDKIT_SUBSTRUCTMATCH_EXPORT bool bondCompat(const Bond* b1, const Bond* b2,
                                            const SubstructMatchParameters& ps);
//! This postprocesses the passed substruct matches and returns
//! the match that has the largest number of non-hydrogen atoms
//! in correspondence of terminal dummy atoms
RDKIT_SUBSTRUCTMATCH_EXPORT const MatchVectType& getMostSubstitutedCoreMatch(
    const ROMol& mol, const ROMol& core,
    const std::vector<MatchVectType>& matches);
//! This returns a copy of the passed substruct matches sorted by decreasing
//! number of non-hydrogen atoms in correspondence of terminal dummy atoms
RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType>
sortMatchesByDegreeOfCoreSubstitution(
    const ROMol& mol, const ROMol& core,
    const std::vector<MatchVectType>& matches);
RDKIT_SUBSTRUCTMATCH_EXPORT bool isAtomTerminalRGroupOrQueryHydrogen(
    const Atom* atom);

}  // namespace RDKit

#endif
