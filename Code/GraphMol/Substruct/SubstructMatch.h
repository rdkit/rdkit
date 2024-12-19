//
//  Copyright (C) 2001-2020 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SUBSTRUCTMATCH_H
#define RD_SUBSTRUCTMATCH_H

// std bits
#include <vector>

#include <unordered_set>
#include <functional>
#include <unordered_map>
#include <cstdint>
#include <string>

#include <boost/dynamic_bitset.hpp>
#if BOOST_VERSION >= 107100
#define RDK_INTERNAL_BITSET_HAS_HASH
#endif

#include <GraphMol/StereoGroup.h>

namespace RDKit {
class ROMol;
class Atom;
class Bond;
class ResonanceMolSupplier;
class MolBundle;

//! \brief used to return matches from substructure searching,
//!   The format is (queryAtomIdx, molAtomIdx)
typedef std::vector<std::pair<int, int>> MatchVectType;

struct RDKIT_SUBSTRUCTMATCH_EXPORT SubstructMatchParameters {
  bool useChirality = false;  //!< Use chirality in determining whether or not
                              //!< atoms/bonds match
  bool useEnhancedStereo = false;  //!< Use enhanced stereochemistry in
                                   //!< determining whether atoms/bonds match
  bool aromaticMatchesConjugated = false;  //!< Aromatic and conjugated bonds
                                           //!< match each other
  bool useQueryQueryMatches = false;  //!< Consider query-query matches, not
                                      //!< just simple matches
  bool useGenericMatchers = false;    //!< Looks for generic atoms in the query
                                      //!< and uses them as part of the matching
  bool recursionPossible = true;      //!< Allow recursive queries
  bool uniquify = true;            //!< uniquify (by atom index) match results
  unsigned int maxMatches = 1000;  //!< maximum number of matches to return
  int numThreads = 1;  //!< number of threads to use when multi-threading
                       //!< is possible. 0 selects the number of
                       //!< concurrent threads supported by the hardware
                       //!< negative values are added to the number of
                       //!< concurrent threads supported by the hardware
  std::vector<std::string> atomProperties;  //!< atom properties that must be
                                            //!< equivalent in order to match
  std::vector<std::string> bondProperties;  //!< bond properties that must be
                                            //!< equivalent in order to match
  std::function<bool(const ROMol &mol,
                     const std::vector<unsigned int> &match)>
      extraFinalCheck;  //!< a function to be called at the end to validate a
                        //!< match
  unsigned int maxRecursiveMatches =
      1000;  //!< maximum number of matches that the recursive substructure
  //!< matching should return
  bool specifiedStereoQueryMatchesUnspecified =
      false;  //!< If set, query atoms and bonds with specified stereochemistry
              //!< will match atoms and bonds with unspecified stereochemistry
  SubstructMatchParameters() {}
};

RDKIT_SUBSTRUCTMATCH_EXPORT void updateSubstructMatchParamsFromJSON(
    SubstructMatchParameters &params, const std::string &json);
RDKIT_SUBSTRUCTMATCH_EXPORT std::string substructMatchParamsToJSON(
    const SubstructMatchParameters &params);

//! Find a substructure match for a query in a molecule
/*!
    \param mol         The ROMol to be searched
    \param query       The query ROMol
    \param matchParams Parameters controlling the matching

    \return The matches, if any

*/
RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ROMol &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());

//! Find all substructure matches for a query in a ResonanceMolSupplier object
/*!
    \param resMolSuppl The ResonanceMolSupplier object to be searched
    \param query       The query ROMol
    \param matchParams Parameters controlling the matching

    \return The matches, if any

*/
RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType> SubstructMatch(
    ResonanceMolSupplier &resMolSuppl, const ROMol &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());

RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType> SubstructMatch(
    const MolBundle &bundle, const ROMol &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());
RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const MolBundle &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());
RDKIT_SUBSTRUCTMATCH_EXPORT std::vector<MatchVectType> SubstructMatch(
    const MolBundle &bundle, const MolBundle &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());

//! Find a substructure match for a query
/*!
    \param mol       The object to be searched
    \param query     The query
    \param matchVect Used to return the match
                     (pre-existing contents will be deleted)
    \param recursionPossible  flags whether or not recursive matches are allowed
    \param useChirality  use atomic CIP codes as part of the comparison
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching

    \return whether or not a match was found

*/
template <typename T1, typename T2>
bool SubstructMatch(T1 &mol, const T2 &query, MatchVectType &matchVect,
                    bool recursionPossible = true, bool useChirality = false,
                    bool useQueryQueryMatches = false) {
  SubstructMatchParameters params;
  params.recursionPossible = recursionPossible;
  params.useChirality = useChirality;
  params.useQueryQueryMatches = useQueryQueryMatches;
  params.maxMatches = 1;
  std::vector<MatchVectType> matchVects = SubstructMatch(mol, query, params);
  if (matchVects.size()) {
    matchVect = matchVects.front();
  } else {
    matchVect.clear();
  }
  return matchVect.size() != 0;
};

//! Find all substructure matches for a query
/*!
    \param mol       The object to be searched
    \param query     The query
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
template <typename T1, typename T2>
unsigned int SubstructMatch(T1 &mol, const T2 &query,
                            std::vector<MatchVectType> &matchVect,
                            bool uniquify = true, bool recursionPossible = true,
                            bool useChirality = false,
                            bool useQueryQueryMatches = false,
                            unsigned int maxMatches = 1000,
                            int numThreads = 1) {
  SubstructMatchParameters params;
  params.uniquify = uniquify;
  params.recursionPossible = recursionPossible;
  params.useChirality = useChirality;
  params.useQueryQueryMatches = useQueryQueryMatches;
  params.maxMatches = maxMatches;
  params.numThreads = numThreads;
  matchVect = SubstructMatch(mol, query, params);
  return static_cast<unsigned int>(matchVect.size());
};

// ----------------------------------------------
//
// find one match in ResonanceMolSupplier object
//
template <>
inline bool SubstructMatch(ResonanceMolSupplier &resMolSupplier,
                           const ROMol &query, MatchVectType &matchVect,
                           bool recursionPossible, bool useChirality,
                           bool useQueryQueryMatches) {
  SubstructMatchParameters params;
  params.recursionPossible = recursionPossible;
  params.useChirality = useChirality;
  params.useQueryQueryMatches = useQueryQueryMatches;
  params.maxMatches = 1;
  std::vector<MatchVectType> matchVects =
      SubstructMatch(resMolSupplier, query, params);
  if (matchVects.size()) {
    matchVect = matchVects.front();
  } else {
    matchVect.clear();
  }
  return matchVect.size() != 0;
}

template <>
inline unsigned int SubstructMatch(ResonanceMolSupplier &resMolSupplier,
                                   const ROMol &query,
                                   std::vector<MatchVectType> &matchVect,
                                   bool uniquify, bool recursionPossible,
                                   bool useChirality, bool useQueryQueryMatches,
                                   unsigned int maxMatches, int numThreads) {
  SubstructMatchParameters params;
  params.uniquify = uniquify;
  params.recursionPossible = recursionPossible;
  params.useChirality = useChirality;
  params.useQueryQueryMatches = useQueryQueryMatches;
  params.maxMatches = maxMatches;
  params.numThreads = numThreads;
  matchVect = SubstructMatch(resMolSupplier, query, params);
  return static_cast<unsigned int>(matchVect.size());
};

//! Class used as a final step to confirm whether or not a given atom->atom
//! mapping is a valid substructure match.
class RDKIT_SUBSTRUCTMATCH_EXPORT MolMatchFinalCheckFunctor {
 public:
  MolMatchFinalCheckFunctor(const ROMol &query, const ROMol &mol,
                            const SubstructMatchParameters &ps);

  bool operator()(const std::uint32_t q_c[], const std::uint32_t m_c[]);

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  const SubstructMatchParameters &d_params;
  std::unordered_map<unsigned int, StereoGroup const *> d_molStereoGroups;
#ifdef RDK_INTERNAL_BITSET_HAS_HASH
  // Boost 1.71 added support for std::hash with dynamic_bitset.
  using HashedStorageType = boost::dynamic_bitset<>;
#else
  // otherwise we use a less elegant solution
  using HashedStorageType = std::string;
#endif
  std::unordered_set<HashedStorageType> matchesSeen;
};

}  // namespace RDKit

#endif
