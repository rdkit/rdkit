//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SUBSTRUCTGENERICS_H
#define RD_SUBSTRUCTGENERICS_H

#include <vector>
#include <functional>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
class ROMol;
class Atom;
class Bond;

namespace SubstructSearch {

// We'd like to be able to correctly interpret what's written by Marvin and
// MarvinJS, so the conditions for these are adapted from the ChemAxon
// documentation for homology groups
// (https://docs.chemaxon.com/display/docs/homology-groups.md)
//
// If I had questions about what the queries should do, I ran example in Reaxys
// with MarvinJS as the sketcher to see what that returns.
//
// I've tried to document deviations or surprises
namespace Generics {

//! Matches alkyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - at least one carbon is present
    - all bonds are single
    - no ring bonds

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool AlkylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches alkenyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C=C
    - no ring bonds

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool AlkenylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches alkynyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C#C
    - no ring bonds

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool AlkynylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches carbocyclic side chains
/*!

  Note: this is Reaxys query type CBC and matches carbocycles

  Conditions:
    - atom is in at least one ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarbocyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches cycloalkyl side chains
/*!

  Note: this is Reaxys query type CAL and is directly equivalent to alkyl,
  except the immediate atom needs to be in a ring.


  Conditions:
    - atom is in at least one ring composed entirely of carbon and connected
      with single bonds
    - atoms in the ring do not have unsaturations (including exocyclic)
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules (i.e. all
      single bonds)


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarbocycloalkylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches cycloalkenyl side chains
/*!

  Note: this is Reaxys query type CEL and matches carbocycles which have at
  least one double or aromatic bond.

  Conditions:
    - atom is in at least one ring composed entirely of carbon and with at least
      one double or aromatic bond
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules (including
      that each ring must have at least one double or aromatic bond)


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarbocycloalkenylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heterocyclic side chains
/*!

  Note: this is Reaxys query type CHC and matches heterocycles

  Conditions:
    - atom is in at least one fused ring with a heteroatom


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool HeterocyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches aryl side chains
/*!

  Note: this is Reaxys query type ARY and matches carbocycles which are aromatic

  Conditions:
    - atom is in at least one aromatic ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarboarylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heteroaryl side chains
/*!

  Note: this is Reaxys query type HAR and matches aromatic heterocycles

  Conditions:
    - atom is in at least one fused aromatic sytem with a heteroatom


*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool HeteroarylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches cyclic side chains
/*!

  Note: this is Reaxys query type CYC and matches cycles

  Conditions:
    - atom is in at least one ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains
/*!

  Note: this is Reaxys query type ACY and matches sidechains with no cycles

  Conditions:
    - no atom in the sidechain is in a ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool AcyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches all-carbon acyclic side chains
/*!

  Note: this is Reaxys query type ABC and matches all-carbon sidechains with no
  cycles

  Conditions:
    - all atoms in the sidechain are carbon
    - no atom in the sidechain is in a ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarboacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches all-carbon acyclic side chains
/*!

  Note: this is Reaxys query type ABC and matches all-carbon sidechains with no
  cycles

  Conditions:
    - all atoms in the sidechain are carbon
    - no atom in the sidechain is in a ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool CarboacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains with at least one heteroatom
/*!

  Note: this is Reaxys query type AHC and matches sidechains with no cycles and
  at least one heteroatom

  Conditions:
    - at least one non-carbon, non-hydrogen atom is in the sidechain
    - no atom in the sidechain is in a ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool HeteroacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches acyclic alkoxy side chains
/*!

  Note: this is Reaxys query type AOX and matches alkoxy sidechains

  Conditions:
    - first atom is an O
    - all other atoms are C
    - all single bonds
    - no atom in the sidechain is in a ring

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool AlkoxyacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches rings without carbon
/*!

  Note: this is Reaxys query type CXX and matches rings which contain no carbon

  Conditions:
    - a ring is present
    - none of the atoms in the fused ring system are carbon

*/
RDKIT_SUBSTRUCTMATCH_EXPORT bool NoCarbonRingAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

const static std::map<
    std::string,
    std::function<bool(const ROMol &, const Atom &, boost::dynamic_bitset<>)>>
    genericMatchers = {
        {"Alkyl", AlkylAtomMatcher},
        {"ALK", AlkylAtomMatcher},
        {"Alkenyl", AlkenylAtomMatcher},
        {"AEL", AlkenylAtomMatcher},
        {"Alkynyl", AlkynylAtomMatcher},
        {"AYL", AlkynylAtomMatcher},
        {"Carbocyclic", CarbocyclicAtomMatcher},
        {"CBC", CarbocyclicAtomMatcher},
        {"Carbocycloalkyl", CarbocycloalkylAtomMatcher},
        {"CAL", CarbocycloalkylAtomMatcher},
        {"Carbocycloalkenyl", CarbocycloalkenylAtomMatcher},
        {"CEL", CarbocycloalkenylAtomMatcher},
        {"Carboaryl", CarboarylAtomMatcher},
        {"ARY", CarboarylAtomMatcher},
        {"Cyclic", CyclicAtomMatcher},
        {"CYC", CyclicAtomMatcher},
        {"Acyclic", AcyclicAtomMatcher},
        {"ACY", AcyclicAtomMatcher},
        {"Carboacyclic", CarboacyclicAtomMatcher},
        {"ABC", CarboacyclicAtomMatcher},
        {"Heteroacyclic", HeteroacyclicAtomMatcher},
        {"AHC", HeteroacyclicAtomMatcher},
        {"Alkoxy", AlkoxyacyclicAtomMatcher},
        {"AOX", AlkoxyacyclicAtomMatcher},
        {"Heterocyclic", HeterocyclicAtomMatcher},
        {"CHC", HeterocyclicAtomMatcher},
        {"Heteroaryl", HeteroarylAtomMatcher},
        {"HAR", HeteroarylAtomMatcher},
        {"NoCarbonRing", NoCarbonRingAtomMatcher},
        {"CXX", NoCarbonRingAtomMatcher},
};
}  // namespace Generics
//! returns false if any of the molecule's generic atoms are not satisfied in
/// the current match
RDKIT_SUBSTRUCTMATCH_EXPORT bool GenericAtomMatcher(
    const ROMol &mol, const ROMol &query,
    const std::vector<unsigned int> &match);
//! sets the apropriate generic query tags based on atom labels and/or SGroups
/*

- Generic query tags found in the atom labels/SGroups will be overwrite existing
generic query tags (if there are any present).
- Any atom labels or SGroups which are converted will be removed
- If the both atom labels and SGroups are being used and an atom has generic
query tags in both, the one from the SGroup will be used.
- Generic query tags not found in Generics::genericMatchers will be ignored

*/
RDKIT_SUBSTRUCTMATCH_EXPORT void SetGenericQueriesFromProperties(
    ROMol &mol, bool useAtomLabels = true, bool useSGroups = true);
}  // namespace SubstructSearch
}  // namespace RDKit

#endif
