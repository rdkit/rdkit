//
//  Copyright (C) 2001-2015 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SUBSTRUCTMATCH_H__
#define _RD_SUBSTRUCTMATCH_H__

// std bits
#include <vector>

namespace RDKit {
class ROMol;
class Atom;
class Bond;
class ResonanceMolSupplier;

//! \brief used to return matches from substructure searching,
//!   The format is (queryAtomIdx, molAtomIdx)
typedef std::vector<std::pair<int, int> > MatchVectType;

//! Find a substructure match for a query in a molecule
/*!
    \param mol       The ROMol to be searched
    \param query     The query ROMol
    \param matchVect Used to return the match
                     (pre-existing contents will be deleted)
    \param recursionPossible  flags whether or not recursive matches are allowed
    \param useChirality  use atomic CIP codes as part of the comparison
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching

    \return whether or not a match was found

*/
bool SubstructMatch(const ROMol &mol, const ROMol &query,
                    MatchVectType &matchVect, bool recursionPossible = true,
                    bool useChirality = false,
                    bool useQueryQueryMatches = false);

//! Find a substructure match for a query in a ResonanceMolSupplier object
/*!
    \param resMolSuppl The ResonanceMolSupplier object to be searched
    \param query       The query ROMol
    \param matchVect   Used to return the match
                       (pre-existing contents will be deleted)
    \param recursionPossible  flags whether or not recursive matches are allowed
    \param useChirality  use atomic CIP codes as part of the comparison
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching

    \return whether or not a match was found

*/
bool SubstructMatch(ResonanceMolSupplier &resMolSuppl, const ROMol &query,
                    MatchVectType &matchVect, bool recursionPossible = true,
                    bool useChirality = false,
                    bool useQueryQueryMatches = false);

//! Find all substructure matches for a query in a molecule
/*!
    \param mol       The ROMol to be searched
    \param query     The query ROMol
    \param matchVect Used to return the matches
                     (pre-existing contents will be deleted)
    \param uniquify  Toggles uniquification (by atom index) of the results
    \param recursionPossible  flags whether or not recursive matches are allowed
    \param useChirality  use atomic CIP codes as part of the comparison
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
    \param maxMatches  The maximum number of matches that will be returned.
                       In high-symmetry cases with medium-sized molecules, it is
   very
                       easy to end up with a combinatorial explosion in the
   number of
                       possible matches. This argument prevents that from having
                       unintended consequences

    \return the number of matches found

*/
unsigned int SubstructMatch(const ROMol &mol, const ROMol &query,
                            std::vector<MatchVectType> &matchVect,
                            bool uniquify = true, bool recursionPossible = true,
                            bool useChirality = false,
                            bool useQueryQueryMatches = false,
                            unsigned int maxMatches = 1000);

//! Find all substructure matches for a query in a ResonanceMolSupplier object
/*!
    \param resMolSuppl The ResonanceMolSupplier object to be searched
    \param query       The query ROMol
    \param matchVect   Used to return the matches
                       (pre-existing contents will be deleted)
    \param uniquify  Toggles uniquification (by atom index) of the results
    \param recursionPossible  flags whether or not recursive matches are allowed
    \param useChirality  use atomic CIP codes as part of the comparison
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
    \param maxMatches  The maximum number of matches that will be returned.
                       In high-symmetry cases with medium-sized molecules, it is
   very
                       easy to end up with a combinatorial explosion in the
   number of
                       possible matches. This argument prevents that from having
                       unintended consequences
    \param numThreads  The number of threads used during the search
                       (defaults to 1; 0 selects the number of
                       concurrent threads supported by the hardware;
                       negative values are added to the number of
                       concurrent threads supported by the
                       hardware)

    \return the number of matches found

*/
unsigned int SubstructMatch(ResonanceMolSupplier &resMolSuppl,
                            const ROMol &query,
                            std::vector<MatchVectType> &matchVect,
                            bool uniquify = false,
                            bool recursionPossible = true,
                            bool useChirality = false,
                            bool useQueryQueryMatches = false,
                            unsigned int maxMatches = 1000, int numThreads = 1);
}

#endif
