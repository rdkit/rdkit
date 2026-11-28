//
//
//  Copyright (C) 2018-2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstanceGroup.h"
#include "ROMol.h"
#include "RWMol.h"
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

namespace RDKit {

namespace {

template <class T>
void remove_element(std::vector<T> &container, unsigned int element) {
  auto pos = std::find(container.begin(), container.end(), element);
  if (pos != container.end()) {
    container.erase(pos);
  }
}
}  // namespace

SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
    : RDProps(), dp_mol(owning_mol) {
  PRECONDITION(owning_mol, "supplied owning molecule is bad");

  // TYPE is required to be set , as other properties will depend on it.
  setProp<std::string>("TYPE", type);
}

void SubstanceGroup::setOwningMol(ROMol *mol) {
  PRECONDITION(mol, "owning molecule is nullptr");

  dp_mol = mol;
}

unsigned int SubstanceGroup::getIndexInMol() const {
  PRECONDITION(dp_mol, "SubstanceGroup is not owned by any molecule");

  const auto &sgroups = getSubstanceGroups(*dp_mol);
  CHECK_INVARIANT(!sgroups.empty(),
                  "No SubstanceGroups found on owning molecule");

  auto match_sgroup = [&](const SubstanceGroup &sg) { return this == &sg; };
  auto sgroupItr = std::find_if(sgroups.begin(), sgroups.end(), match_sgroup);

  if (sgroupItr == sgroups.end()) {
    std::ostringstream errout;
    errout << "Unable to find own index in owning mol SubstanceGroup collection"
           << std::endl;
    throw SubstanceGroupException(errout.str());
  }

  return sgroupItr - sgroups.begin();
}

void SubstanceGroup::setAtoms(std::vector<unsigned int> atoms) {
  PRECONDITION(dp_mol, "bad mol");
  auto natoms = dp_mol->getNumAtoms();
  for (auto atomIdx : atoms) {
    if (atomIdx >= natoms) {
      std::ostringstream errout;
      errout << "Atom index " << atomIdx << " is out of range";
      throw ValueErrorException(errout.str());
    }
  }
  d_atoms = std::move(atoms);
}
void SubstanceGroup::setParentAtoms(std::vector<unsigned int> patoms) {
  PRECONDITION(dp_mol, "bad mol");
  auto natoms = dp_mol->getNumAtoms();
  for (auto atomIdx : patoms) {
    if (atomIdx >= natoms) {
      std::ostringstream errout;
      errout << "Atom index " << atomIdx << " is out of range";
      throw ValueErrorException(errout.str());
    }
  }
  d_patoms = std::move(patoms);
}
void SubstanceGroup::setBonds(std::vector<unsigned int> bonds) {
  PRECONDITION(dp_mol, "bad mol");
  auto nbonds = dp_mol->getNumBonds();
  for (auto bondIdx : bonds) {
    if (bondIdx >= nbonds) {
      std::ostringstream errout;
      errout << "Bond index " << bondIdx << " is out of range";
      throw ValueErrorException(errout.str());
    }
  }
  d_bonds = std::move(bonds);
}

void SubstanceGroup::addAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  if (idx >= dp_mol->getNumAtoms()) {
    throw ValueErrorException("Atom index out of range");
  }
  d_atoms.push_back(idx);
}

void SubstanceGroup::addAtomWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");
  Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
  PRECONDITION(atom, "atom not found");
  d_atoms.push_back(atom->getIdx());
}

void SubstanceGroup::addParentAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  if (idx >= dp_mol->getNumAtoms()) {
    throw ValueErrorException("Atom index out of range");
  }
  if (std::find(d_atoms.begin(), d_atoms.end(), idx) == d_atoms.end()) {
    std::ostringstream errout;
    errout << "Atom " << idx << " is not a member of current SubstanceGroup";
    throw SubstanceGroupException(errout.str());
  }

  d_patoms.push_back(idx);
}

void SubstanceGroup::addParentAtomWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");

  Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
  unsigned int idx = atom->getIdx();
  if (std::find(d_atoms.begin(), d_atoms.end(), idx) == d_atoms.end()) {
    std::ostringstream errout;
    errout << "Atom with bookmark " << mark
           << " is not a member of current SubstanceGroup ";
    throw SubstanceGroupException(errout.str());
  }

  d_patoms.push_back(idx);
}

void SubstanceGroup::addBondWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  if (idx >= dp_mol->getNumBonds()) {
    throw ValueErrorException("Bond index out of range");
  }

  d_bonds.push_back(idx);
}

void SubstanceGroup::addBondWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");
  Bond *bond = dp_mol->getUniqueBondWithBookmark(mark);
  d_bonds.push_back(bond->getIdx());
}

void SubstanceGroup::removeAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  remove_element(d_atoms, idx);
}

void SubstanceGroup::removeParentAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  remove_element(d_patoms, idx);
}

void SubstanceGroup::removeBondWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  remove_element(d_bonds, idx);
}

void SubstanceGroup::addBracket(const SubstanceGroup::Bracket &bracket) {
  d_brackets.push_back(bracket);
}

void SubstanceGroup::addCState(unsigned int bondIdx,
                               const RDGeom::Point3D &vector) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_bonds.empty(), "no bonds");

  if (getBondType(bondIdx) != SubstanceGroup::BondType::XBOND) {
    std::ostringstream errout;
    errout << "Bond with index " << bondIdx
           << " is not an XBOND for current SubstanceGroup";
    throw SubstanceGroupException(errout.str());
  }
  d_cstates.push_back({bondIdx, vector});
}

void SubstanceGroup::addAttachPoint(unsigned int aIdx, int lvIdx,
                                    const std::string &idStr) {
  d_saps.push_back({aIdx, lvIdx, idStr});
}

//! check if the bond is SubstanceGroup XBOND or CBOND
SubstanceGroup::BondType SubstanceGroup::getBondType(
    unsigned int bondIdx) const {
  PRECONDITION(
      std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end(),
      "bond is not part of the SubstanceGroup")

  auto bond = dp_mol->getBondWithIdx(bondIdx);
  bool begin_atom_in_sgroup =
      std::find(d_atoms.begin(), d_atoms.end(), bond->getBeginAtomIdx()) !=
      d_atoms.end();
  bool end_atom_in_sgroup = std::find(d_atoms.begin(), d_atoms.end(),
                                      bond->getEndAtomIdx()) != d_atoms.end();

  if (begin_atom_in_sgroup && end_atom_in_sgroup) {
    return SubstanceGroup::BondType::CBOND;
  } else if (begin_atom_in_sgroup || end_atom_in_sgroup) {
    return SubstanceGroup::BondType::XBOND;
  } else {
    std::ostringstream errout;
    errout << "Neither beginning nor ending atoms of bond " << bond->getIdx()
           << " is in this SubstanceGroup.";
    throw SubstanceGroupException(errout.str());
  }
}

bool SubstanceGroup::adjustToRemovedAtom(unsigned int atomIdx) {
  bool res = false;
  for (auto &aid : d_atoms) {
    if (aid == atomIdx) {
      throw SubstanceGroupException(
          "adjustToRemovedAtom() called on SubstanceGroup which contains the "
          "atom");
    }
    if (aid > atomIdx) {
      res = true;
      --aid;
    }
  }
  for (auto &aid : d_patoms) {
    if (aid == atomIdx) {
      throw SubstanceGroupException(
          "adjustToRemovedAtom() called on SubstanceGroup which contains the "
          "atom");
    }
    if (aid > atomIdx) {
      res = true;
      --aid;
    }
  }
  for (auto &ap : d_saps) {
    if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
      throw SubstanceGroupException(
          "adjustToRemovedAtom() called on SubstanceGroup which contains the "
          "atom");
    }
    if (ap.aIdx > atomIdx) {
      res = true;
      --ap.aIdx;
    }
    if (ap.lvIdx > rdcast<int>(atomIdx)) {
      res = true;
      --ap.lvIdx;
    }
  }

  return res;
}

bool SubstanceGroup::adjustToRemovedBond(unsigned int bondIdx) {
  bool res = false;
  for (auto &bid : d_bonds) {
    if (bid == bondIdx) {
      throw SubstanceGroupException(
          "adjustToRemovedBond() called on SubstanceGroup which contains the "
          "bond");
    }
    if (bid > bondIdx) {
      res = true;
      --bid;
    }
  }
  for (auto &cs : d_cstates) {
    if (cs.bondIdx == bondIdx) {
      throw SubstanceGroupException(
          "adjustToRemovedBond() called on SubstanceGroup which contains the "
          "bond");
    }
    if (cs.bondIdx > bondIdx) {
      res = true;
      --cs.bondIdx;
    }
  }
  return res;
}

bool SubstanceGroup::includesAtom(unsigned int atomIdx) const {
  if (std::find(d_atoms.begin(), d_atoms.end(), atomIdx) != d_atoms.end()) {
    return true;
  }
  if (std::find(d_patoms.begin(), d_patoms.end(), atomIdx) != d_patoms.end()) {
    return true;
  }
  for (const auto &ap : d_saps) {
    if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
      return true;
    }
  }
  return false;
}

bool SubstanceGroup::includesBond(unsigned int bondIdx) const {
  if (std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end()) {
    return true;
  }
  for (const auto &cs : d_cstates) {
    if (cs.bondIdx == bondIdx) {
      return true;
    }
  }
  return false;
}

namespace {
bool containsIdx(const std::vector<unsigned int> &sortedIndices,
                 unsigned int idx) {
  return std::binary_search(sortedIndices.begin(), sortedIndices.end(), idx);
}
// the index that idx becomes after removing the (sorted) indices
unsigned int adjustedIdx(const std::vector<unsigned int> &sortedIndices,
                         unsigned int idx) {
  const auto numBefore =
      std::lower_bound(sortedIndices.begin(), sortedIndices.end(), idx) -
      sortedIndices.begin();
  return idx - static_cast<unsigned int>(numBefore);
}
}  // namespace

bool SubstanceGroup::adjustToRemovedAtoms(
    const std::vector<unsigned int> &atomIndices) {
  if (atomIndices.empty()) {
    return false;
  }
  bool res = false;
  auto adjust = [&atomIndices, &res](unsigned int &aid) {
    if (containsIdx(atomIndices, aid)) {
      throw SubstanceGroupException(
          "adjustToRemovedAtoms() called on SubstanceGroup which contains "
          "one of the atoms");
    }
    const auto newIdx = adjustedIdx(atomIndices, aid);
    if (newIdx != aid) {
      res = true;
      aid = newIdx;
    }
  };
  for (auto &aid : d_atoms) {
    adjust(aid);
  }
  for (auto &aid : d_patoms) {
    adjust(aid);
  }
  for (auto &ap : d_saps) {
    adjust(ap.aIdx);
    if (ap.lvIdx >= 0) {
      auto lvIdx = rdcast<unsigned int>(ap.lvIdx);
      adjust(lvIdx);
      ap.lvIdx = rdcast<int>(lvIdx);
    }
  }
  return res;
}

bool SubstanceGroup::includesAnyAtom(
    const std::vector<unsigned int> &atomIndices) const {
  if (atomIndices.empty()) {
    return false;
  }
  for (auto aid : d_atoms) {
    if (containsIdx(atomIndices, aid)) {
      return true;
    }
  }
  for (auto aid : d_patoms) {
    if (containsIdx(atomIndices, aid)) {
      return true;
    }
  }
  for (const auto &ap : d_saps) {
    if (containsIdx(atomIndices, ap.aIdx) ||
        (ap.lvIdx >= 0 &&
         containsIdx(atomIndices, rdcast<unsigned int>(ap.lvIdx)))) {
      return true;
    }
  }
  return false;
}

bool SubstanceGroup::adjustToRemovedBonds(
    const std::vector<unsigned int> &bondIndices) {
  if (bondIndices.empty()) {
    return false;
  }
  bool res = false;
  auto adjust = [&bondIndices, &res](unsigned int &bid) {
    if (containsIdx(bondIndices, bid)) {
      throw SubstanceGroupException(
          "adjustToRemovedBonds() called on SubstanceGroup which contains "
          "one of the bonds");
    }
    const auto newIdx = adjustedIdx(bondIndices, bid);
    if (newIdx != bid) {
      res = true;
      bid = newIdx;
    }
  };
  for (auto &bid : d_bonds) {
    adjust(bid);
  }
  for (auto &cs : d_cstates) {
    adjust(cs.bondIdx);
  }
  return res;
}

bool SubstanceGroup::includesAnyBond(
    const std::vector<unsigned int> &bondIndices) const {
  if (bondIndices.empty()) {
    return false;
  }
  for (auto bid : d_bonds) {
    if (containsIdx(bondIndices, bid)) {
      return true;
    }
  }
  for (const auto &cs : d_cstates) {
    if (containsIdx(bondIndices, cs.bondIdx)) {
      return true;
    }
  }
  return false;
}

bool SubstanceGroupChecks::isValidType(const std::string &type) {
  return std::find(SubstanceGroupChecks::sGroupTypes.begin(),
                   SubstanceGroupChecks::sGroupTypes.end(),
                   type) != SubstanceGroupChecks::sGroupTypes.end();
}

bool SubstanceGroupChecks::isValidSubType(const std::string &type) {
  return std::find(SubstanceGroupChecks::sGroupSubtypes.begin(),
                   SubstanceGroupChecks::sGroupSubtypes.end(),
                   type) != SubstanceGroupChecks::sGroupSubtypes.end();
}

bool SubstanceGroupChecks::isValidConnectType(const std::string &type) {
  return std::find(SubstanceGroupChecks::sGroupConnectTypes.begin(),
                   SubstanceGroupChecks::sGroupConnectTypes.end(),
                   type) != SubstanceGroupChecks::sGroupConnectTypes.end();
}

bool SubstanceGroupChecks::isValidClass(const std::string &sgroupClass) {
  return std::find(SubstanceGroupChecks::sGroupClasses.begin(),
                   SubstanceGroupChecks::sGroupClasses.end(),
                   sgroupClass) != SubstanceGroupChecks::sGroupClasses.end();
}

bool SubstanceGroupChecks::isSubstanceGroupIdFree(const ROMol &mol,
                                                  unsigned int id) {
  auto match_sgroup = [id](const SubstanceGroup &sg) {
    unsigned int storedId;
    return sg.getPropIfPresent("ID", storedId) && id == storedId;
  };

  const auto &sgroups = getSubstanceGroups(mol);
  return std::find_if(sgroups.begin(), sgroups.end(), match_sgroup) ==
         sgroups.end();
}

std::vector<SubstanceGroup> &getSubstanceGroups(ROMol &mol) {
  return mol.d_sgroups;
}
const std::vector<SubstanceGroup> &getSubstanceGroups(const ROMol &mol) {
  return mol.d_sgroups;
}

unsigned int addSubstanceGroup(ROMol &mol, SubstanceGroup sgroup) {
  sgroup.setOwningMol(&mol);

  auto &&sgroups = getSubstanceGroups(mol);
  unsigned int id = sgroups.size();

  sgroups.push_back(std::move(sgroup));

  return id;
}

namespace {
bool removedParentInHierarchy(
    unsigned int idx, const std::vector<SubstanceGroup> &sgs,
    const boost::dynamic_bitset<> &toRemove,
    const std::map<unsigned int, unsigned int> &indexLookup) {
  PRECONDITION(idx < sgs.size(), "cannot find SubstanceGroup");
  if (toRemove[idx]) {
    return true;
  }

  unsigned int parent;
  if (sgs[idx].getPropIfPresent("PARENT", parent)) {
    auto piter = indexLookup.find(parent);
    if (piter != indexLookup.end()) {
      return removedParentInHierarchy(piter->second, sgs, toRemove,
                                      indexLookup);
    }
  }
  return false;
}

// Removes the SubstanceGroups for which includes(sg) is true, along with their
// descendants, and calls adjust(sg) on the others.
template <typename IncludesFunc, typename AdjustFunc>
void removeSubstanceGroupsReferencing(RWMol &mol, IncludesFunc includes,
                                      AdjustFunc adjust) {
  auto &sgs = getSubstanceGroups(mol);
  if (!sgs.empty()) {
    // first collect the ones that should be removed
    boost::dynamic_bitset<> toRemove(sgs.size());
    unsigned int nRemoved = 0;
    bool parentsPresent = false;
    for (unsigned int i = 0; i < sgs.size(); ++i) {
      if (!parentsPresent && sgs[i].hasProp("PARENT")) {
        parentsPresent = true;
      }
      if (includes(sgs[i])) {
        toRemove.set(i);
        ++nRemoved;
      }
    }

    // if we're going to be removing anything and there are PARENTS present,
    // we need to build a lookup map between index->position in original array
    std::map<unsigned int, unsigned int> indexLookup;
    if (parentsPresent && nRemoved) {
      for (unsigned int i = 0; i < sgs.size(); ++i) {
        unsigned int index;
        if (sgs[i].getPropIfPresent("index", index)) {
          indexLookup[index] = i;
        }
      }
    }
    // now go through and keep everything that shouldn't be removed
    // and who doesn't have a PARENT that should be removed in their hierarchy
    std::vector<SubstanceGroup> newsgs;
    newsgs.reserve(sgs.size() - nRemoved);
    unsigned int i = 0;
    for (auto &&sg : sgs) {
      if (!toRemove[i]) {
        // we might be keeping it. Check the parent
        if (!parentsPresent || !sg.hasProp("PARENT")) {
          adjust(sg);
          newsgs.push_back(std::move(sg));
        } else if (parentsPresent) {
          unsigned int parent;
          // has our parent been removed?
          if (sg.getPropIfPresent("PARENT", parent)) {
            auto piter = indexLookup.find(parent);
            bool keepIt = false;
            if (piter == indexLookup.end()) {
              // our parent isn't around, so it isn't being removed
              // note: this is an odd case and probably shouldn't happen, but
              // this isn't the place to enforce that
              keepIt = true;
            } else if (!toRemove[piter->second]) {
              // our parent isn't being removed, recursively check up through
              // parents to see if we find any that are being removed:
              if (!removedParentInHierarchy(piter->second, sgs, toRemove,
                                            indexLookup)) {
                keepIt = true;
              }
            }
            if (keepIt) {
              adjust(sg);
              newsgs.push_back(std::move(sg));
            }
          }
        }
      }
      ++i;
    }
    sgs = std::move(newsgs);
  }
}
}  // namespace
void removeSubstanceGroupsReferencingAtom(RWMol &mol, unsigned int idx) {
  // Delete substance groups containing this atom. It could be that it's ok to
  // keep it, but we just don't know
  removeSubstanceGroupsReferencing(
      mol, [idx](const SubstanceGroup &sg) { return sg.includesAtom(idx); },
      [idx](SubstanceGroup &sg) { sg.adjustToRemovedAtom(idx); });
}

void removeSubstanceGroupsReferencingBond(RWMol &mol, unsigned int idx) {
  // Delete substance groups containing this bond. It could be that it's ok to
  // keep it, but we just don't know
  removeSubstanceGroupsReferencing(
      mol, [idx](const SubstanceGroup &sg) { return sg.includesBond(idx); },
      [idx](SubstanceGroup &sg) { sg.adjustToRemovedBond(idx); });
}

void removeSubstanceGroupsReferencingAtoms(
    RWMol &mol, const std::vector<unsigned int> &atomIndices) {
  PRECONDITION(std::is_sorted(atomIndices.begin(), atomIndices.end()),
               "atom indices must be sorted");
  if (atomIndices.empty()) {
    return;
  }
  removeSubstanceGroupsReferencing(
      mol,
      [&atomIndices](const SubstanceGroup &sg) {
        return sg.includesAnyAtom(atomIndices);
      },
      [&atomIndices](SubstanceGroup &sg) {
        sg.adjustToRemovedAtoms(atomIndices);
      });
}

void removeSubstanceGroupsReferencingBonds(
    RWMol &mol, const std::vector<unsigned int> &bondIndices) {
  PRECONDITION(std::is_sorted(bondIndices.begin(), bondIndices.end()),
               "bond indices must be sorted");
  if (bondIndices.empty()) {
    return;
  }
  removeSubstanceGroupsReferencing(
      mol,
      [&bondIndices](const SubstanceGroup &sg) {
        return sg.includesAnyBond(bondIndices);
      },
      [&bondIndices](SubstanceGroup &sg) {
        sg.adjustToRemovedBonds(bondIndices);
      });
}

}  // namespace RDKit

std::ostream &operator<<(std::ostream &target,
                         const RDKit::SubstanceGroup &sgroup) {
  target << sgroup.getIndexInMol() << ' '
         << sgroup.getProp<std::string>("TYPE");

  auto brackets = sgroup.getBrackets();
  if (!brackets.empty()) {
    target << " Brk: " << brackets.size();
  }

  auto cstates = sgroup.getCStates();
  if (!cstates.empty()) {
    target << " CSt: " << cstates.size();
  }

  auto attachpts = sgroup.getAttachPoints();
  if (!attachpts.empty()) {
    target << " AtPt: " << attachpts.size();
  }

  auto atoms = sgroup.getAtoms();
  if (!atoms.empty()) {
    target << " Atoms: { ";
    for (auto atom_idx : atoms) {
      target << atom_idx << ' ';
    }
    target << '}';
  }

  auto patoms = sgroup.getParentAtoms();
  if (!patoms.empty()) {
    target << " PAtoms: { ";
    for (auto atom_idx : patoms) {
      target << atom_idx << ' ';
    }
    target << '}';
  }

  auto bonds = sgroup.getBonds();
  if (!bonds.empty()) {
    target << " Bonds: { ";
    for (auto bond_idx : bonds) {
      target << bond_idx << ' ';
    }
    target << '}';
  }

  return target;
}
