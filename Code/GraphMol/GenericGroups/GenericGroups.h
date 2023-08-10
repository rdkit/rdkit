//
//  Copyright (C) 2021-2023 Greg Landrum and RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#ifndef RD_GENERICGROUPS_H
#define RD_GENERICGROUPS_H

#include <vector>
#include <functional>
#include <map>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
class ROMol;
class Atom;
class Bond;

namespace GenericGroups {
// We'd like to be able to correctly interpret what's written by Marvin and
// MarvinJS, so the conditions for these are adapted from the ChemAxon
// documentation for homology groups
// (https://docs.chemaxon.com/display/docs/homology-groups.md)
//
// If I had questions about what the queries should do, I ran example in Reaxys
// with MarvinJS as the sketcher to see what that returns.
//
// I've tried to document deviations or surprises

namespace Matchers {

//! Matches any group as a side chain
/*!

  Note: this is Reaxys query type G and matches any sidechain

  Conditions:
    - at least one non-hydrogen atom is in the sidechain

*/
RDKIT_GENERICGROUPS_EXPORT bool GroupAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches any group as a side chain including just an H atom
/*!

  Note: this is Reaxys query type GH and matches any sidechain

  Conditions:
    - none

*/
RDKIT_GENERICGROUPS_EXPORT bool GroupHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches any group as a side chain
/*!

  Note: this is Reaxys query type G* and matches any sidechain that has a ring
  closure

  Conditions:
    - at least one non-hydrogen atom is in the sidechain
    - at least one ring closure

*/
RDKIT_GENERICGROUPS_EXPORT bool GroupStarAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches any group as a side chain that has a ring closure or just an H atom
/*!

  Note: this is Reaxys query type GH* and matches any sidechains

  Conditions:
    - at least one ring closure
    - OR
    - the entire group is just an H atom


*/
RDKIT_GENERICGROUPS_EXPORT bool GroupStarHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches alkyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - at least one carbon is present
    - all bonds are single
    - no ring bonds

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches alkyl side chains or an H atom
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - all bonds are single
    - no ring bonds

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches alkenyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C=C
    - no ring bonds

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkenylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches alkenyl side chains  or an H
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C=C
    - no ring bonds

*/

RDKIT_GENERICGROUPS_EXPORT bool AlkenylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches alkynyl side chains
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C#C
    - no ring bonds

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkynylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches alkynyl side chains or an H
/*!

  Conditions:
    - side chain consists entirely of carbon or hydrogen
    - contains at least one C#C
    - no ring bonds
    - OR
    - the whole group is an H atom

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkynylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches carbocyclic side chains
/*!

  Note: this is Reaxys query type CBC and matches carbocycles

  Conditions:
    - atom is in at least one ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules


*/
RDKIT_GENERICGROUPS_EXPORT bool CarbocyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches carbocyclic side chains or an H atom
/*!

  Note: this is Reaxys query type CBC and matches carbocycles

  Conditions:
    - atom is in at least one ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules

    - OR the entire group is just an H atom


*/
RDKIT_GENERICGROUPS_EXPORT bool CarbocyclicHAtomMatcher(
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
RDKIT_GENERICGROUPS_EXPORT bool CarbocycloalkylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches cycloalkyl side chains or an H atom
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
    - OR
    - the whole group is an H atom

*/
RDKIT_GENERICGROUPS_EXPORT bool CarbocycloalkylHAtomMatcher(
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
RDKIT_GENERICGROUPS_EXPORT bool CarbocycloalkenylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
//! Matches cycloalkenyl side chains or an H atom
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
RDKIT_GENERICGROUPS_EXPORT bool CarbocycloalkenylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heterocyclic side chains
/*!

  Note: this is Reaxys query type CHC and matches heterocycles

  Conditions:
    - atom is in at least one fused ring with a heteroatom


*/
RDKIT_GENERICGROUPS_EXPORT bool HeterocyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heterocyclic side chains or an H atom
/*!

  Note: this is Reaxys query type CHH and matches heterocycles or an H atom

  Conditions:
    - atom is in at least one fused ring with a heteroatom
    - or the entire group is a single H atom


*/
RDKIT_GENERICGROUPS_EXPORT bool HeterocyclicHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches aryl side chains
/*!

  Note: this is Reaxys query type ARY and matches carbocycles which are aromatic

  Conditions:
    - atom is in at least one aromatic ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules


*/
RDKIT_GENERICGROUPS_EXPORT bool CarboarylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches aryl side chains or an H atom
/*!

  Note: this is Reaxys query type ARH and matches carbocycles which are aromatic
  or an H atom

  Conditions:
    - atom is in at least one aromatic ring composed entirely of carbon
    - atom is not in any rings not compatible with the above conditions
    - additional fused rings in the system must obey the same rules


*/
RDKIT_GENERICGROUPS_EXPORT bool CarboarylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heteroaryl side chains
/*!

  Note: this is Reaxys query type HAR and matches aromatic heterocycles

  Conditions:
    - atom is in at least one fused aromatic sytem with a heteroatom


*/
RDKIT_GENERICGROUPS_EXPORT bool HeteroarylAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches heteroaryl side chains or an H atom
/*!

  Note: this is Reaxys query type HAR and matches aromatic heterocycles

  Conditions:
    - atom is in at least one fused aromatic sytem with a heteroatom
    - or the entire group is an H atom


*/
RDKIT_GENERICGROUPS_EXPORT bool HeteroarylHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches cyclic side chains
/*!

  Note: this is Reaxys query type CYC and matches cycles

  Conditions:
    - atom is in at least one ring

*/
RDKIT_GENERICGROUPS_EXPORT bool CyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches cyclic side chains or an H atom
/*!

  Note: this is Reaxys query type CYH and matches cycles

  Conditions:
    - atom is in at least one ring
    - or the entire group is just an H atom

*/
RDKIT_GENERICGROUPS_EXPORT bool CyclicHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains
/*!

  Note: this is Reaxys query type ACY and matches sidechains with no cycles

  Conditions:
    - no atom in the sidechain is in a ring and the group is NOT just an H atom

*/
RDKIT_GENERICGROUPS_EXPORT bool AcyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains or an H atom
/*!

  Note: this is Reaxys query type ACY and matches sidechains with no cycles

  Conditions:
    - no atom in the sidechain is in a ring

*/
RDKIT_GENERICGROUPS_EXPORT bool AcyclicHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches all-carbon acyclic side chains
/*!

  Note: this is Reaxys query type ABC and matches all-carbon sidechains with no
  cycles

  Conditions:
    - all atoms in the sidechain are carbon
    - no atom in the sidechain is in a ring

*/
RDKIT_GENERICGROUPS_EXPORT bool CarboacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches all-carbon acyclic side chainsor or an H atom
/*!

  Note: this is Reaxys query type ABH and matches all-carbon sidechains with no
  cycles or just a H atom

  Conditions:
    - all atoms in the sidechain are carbon or H
    - no atom in the sidechain is in a ring

*/
RDKIT_GENERICGROUPS_EXPORT bool CarboacyclicHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains with at least one heteroatom
/*!

  Note: this is Reaxys query type AHC and matches sidechains with no cycles and
  at least one heteroatom

  Conditions:
    - at least one non-carbon, non-hydrogen atom is in the sidechain
    - no atom in the sidechain is in a ring

*/
RDKIT_GENERICGROUPS_EXPORT bool HeteroacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches acyclic side chains with at least one heteroatom or an H atom
/*!

  Note: this is Reaxys query type AHC and matches sidechains with no cycles and
  at least one heteroatom

  Conditions:
    - at least one non-carbon, non-hydrogen atom is in the sidechain
    - no atom in the sidechain is in a ring

*/
RDKIT_GENERICGROUPS_EXPORT bool HeteroacyclicHAtomMatcher(
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
RDKIT_GENERICGROUPS_EXPORT bool AlkoxyacyclicAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);
/*!

Note: this is Reaxys query type AOH and matches alkoxy sidechains or a hydrogen

Conditions:
- first atom is an O
- all other atoms are C
- all single bonds
- no atom in the sidechain is in a ring
- OR
- the whole group is just an H atom

*/
RDKIT_GENERICGROUPS_EXPORT bool AlkoxyacyclicHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches rings without carbon
/*!

  Note: this is Reaxys query type CXX and matches rings which contain no carbon

  Conditions:
    - a ring is present
    - none of the atoms in the fused ring system are carbon

*/
RDKIT_GENERICGROUPS_EXPORT bool NoCarbonRingAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

//! Matches rings without carbon or just an H
/*!

  Note: this is Reaxys query type CXH and matches rings which contain no carbon

  Conditions:
    - a ring is present
    - none of the atoms in the fused ring system are carbon
    - OR
    - the entire group is just an H atom

*/

RDKIT_GENERICGROUPS_EXPORT bool NoCarbonRingHAtomMatcher(
    const ROMol &mol, const Atom &atom, boost::dynamic_bitset<> ignore);

}  // namespace Matchers
const static std::map<
    std::string,
    std::function<bool(const ROMol &, const Atom &, boost::dynamic_bitset<>)>>
    genericMatchers = {
        {"Group", Matchers::GroupAtomMatcher},
        {"G", Matchers::GroupAtomMatcher},
        {"GroupH", Matchers::GroupHAtomMatcher},
        {"GH", Matchers::GroupHAtomMatcher},
        {"Group*", Matchers::GroupStarAtomMatcher},
        {"G*", Matchers::GroupStarAtomMatcher},
        {"GroupH*", Matchers::GroupStarHAtomMatcher},
        {"GH*", Matchers::GroupStarHAtomMatcher},
        {"Alkyl", Matchers::AlkylAtomMatcher},
        {"ALK", Matchers::AlkylAtomMatcher},
        {"AlkylH", Matchers::AlkylHAtomMatcher},
        {"ALH", Matchers::AlkylHAtomMatcher},
        {"Alkenyl", Matchers::AlkenylAtomMatcher},
        {"AEL", Matchers::AlkenylAtomMatcher},
        {"AlkenylH", Matchers::AlkenylHAtomMatcher},
        {"AEH", Matchers::AlkenylHAtomMatcher},
        {"Alkynyl", Matchers::AlkynylAtomMatcher},
        {"AYL", Matchers::AlkynylAtomMatcher},
        {"AlkynylH", Matchers::AlkynylHAtomMatcher},
        {"AYH", Matchers::AlkynylHAtomMatcher},
        {"Carbocyclic", Matchers::CarbocyclicAtomMatcher},
        {"CBC", Matchers::CarbocyclicAtomMatcher},
        {"CarbocyclicH", Matchers::CarbocyclicHAtomMatcher},
        {"CBH", Matchers::CarbocyclicHAtomMatcher},
        {"Carbocycloalkyl", Matchers::CarbocycloalkylAtomMatcher},
        {"CAL", Matchers::CarbocycloalkylAtomMatcher},
        {"CarbocycloalkylH", Matchers::CarbocycloalkylHAtomMatcher},
        {"CAH", Matchers::CarbocycloalkylHAtomMatcher},
        {"Carbocycloalkenyl", Matchers::CarbocycloalkenylAtomMatcher},
        {"CEL", Matchers::CarbocycloalkenylAtomMatcher},
        {"CarbocycloalkenylH", Matchers::CarbocycloalkenylHAtomMatcher},
        {"CEH", Matchers::CarbocycloalkenylHAtomMatcher},
        {"Carboaryl", Matchers::CarboarylAtomMatcher},
        {"ARY", Matchers::CarboarylAtomMatcher},
        {"CarboarylH", Matchers::CarboarylHAtomMatcher},
        {"ARH", Matchers::CarboarylHAtomMatcher},
        {"Cyclic", Matchers::CyclicAtomMatcher},
        {"CYC", Matchers::CyclicAtomMatcher},
        {"CyclicH", Matchers::CyclicHAtomMatcher},
        {"CYH", Matchers::CyclicHAtomMatcher},
        {"Acyclic", Matchers::AcyclicAtomMatcher},
        {"ACY", Matchers::AcyclicAtomMatcher},
        {"AcyclicH", Matchers::AcyclicHAtomMatcher},
        {"ACH", Matchers::AcyclicHAtomMatcher},
        {"Carboacyclic", Matchers::CarboacyclicAtomMatcher},
        {"ABC", Matchers::CarboacyclicAtomMatcher},
        {"CarboacyclicH", Matchers::CarboacyclicHAtomMatcher},
        {"ABH", Matchers::CarboacyclicHAtomMatcher},
        {"Heteroacyclic", Matchers::HeteroacyclicAtomMatcher},
        {"AHC", Matchers::HeteroacyclicAtomMatcher},
        {"HeteroacyclicH", Matchers::HeteroacyclicHAtomMatcher},
        {"AHH", Matchers::HeteroacyclicHAtomMatcher},
        {"Alkoxy", Matchers::AlkoxyacyclicAtomMatcher},
        {"AOX", Matchers::AlkoxyacyclicAtomMatcher},
        {"AlkoxyH", Matchers::AlkoxyacyclicHAtomMatcher},
        {"AOH", Matchers::AlkoxyacyclicHAtomMatcher},
        {"Heterocyclic", Matchers::HeterocyclicAtomMatcher},
        {"CHC", Matchers::HeterocyclicAtomMatcher},
        {"HeterocyclicH", Matchers::HeterocyclicHAtomMatcher},
        {"CHH", Matchers::HeterocyclicHAtomMatcher},
        {"Heteroaryl", Matchers::HeteroarylAtomMatcher},
        {"HAR", Matchers::HeteroarylAtomMatcher},
        {"HeteroarylH", Matchers::HeteroarylHAtomMatcher},
        {"HAH", Matchers::HeteroarylHAtomMatcher},
        {"NoCarbonRing", Matchers::NoCarbonRingAtomMatcher},
        {"CXX", Matchers::NoCarbonRingAtomMatcher},
        {"NoCarbonRingH", Matchers::NoCarbonRingHAtomMatcher},
        {"CXH", Matchers::NoCarbonRingHAtomMatcher}
};

// This is an extension of adjustQueryProperties from GraphMol that allows the search of generic groups 
RDKIT_GENERICGROUPS_EXPORT ROMol *adjustQueryPropertiesWithGenericGroups(
    const ROMol &mol,
    const MolOps::AdjustQueryParameters *inParams=nullptr);

//! returns false if any of the molecule's generic atoms are not satisfied in
/// the current match
RDKIT_GENERICGROUPS_EXPORT bool genericAtomMatcher(
    const ROMol &mol, const ROMol &query,
    const std::vector<unsigned int> &match);
//! sets the apropriate generic query tags based on atom labels and/or SGroups
/*

- Generic query tags found in the atom labels/SGroups will be overwrite existing
generic query tags (if there are any present).
- only SUP SGroups are considered
- Any atom labels or SGroups which are converted will be removed
- If both atom labels and SGroups are being used and an atom has generic
query tags in both, the one from the SGroup will be used.
- Generic query tags not found in GenericGroups::genericMatchers will be ignored

*/
RDKIT_GENERICGROUPS_EXPORT void setGenericQueriesFromProperties(
    ROMol &mol, bool useAtomLabels = true, bool useSGroups = true);
RDKIT_GENERICGROUPS_EXPORT void convertGenericQueriesToSubstanceGroups(
    ROMol &mol);
}  // namespace GenericGroups
}  // namespace RDKit

#endif
