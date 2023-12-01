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

void SubstanceGroup::addAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(dp_mol->getAtomWithIdx(idx), "wrong atom index");

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
  PRECONDITION(dp_mol->getBondWithIdx(idx), "wrong bond index");

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

namespace SubstanceGroupChecks {
const std::vector<const char *> sGroupTypes = {
    // polymer sgroups:
    "SRU", "MON", "COP", "CRO", "GRA", "MOD", "MER", "ANY",
    // formulations/mixtures:
    "COM", "MIX", "FOR",
    // other
    "SUP", "MUL", "DAT", "GEN"};

bool isValidType(const std::string &type) {
  return std::find(sGroupTypes.begin(), sGroupTypes.end(), type) !=
         sGroupTypes.end();
}

bool isValidSubType(const std::string &type) {
  return type == "ALT" || type == "RAN" || type == "BLO";
}

bool isValidConnectType(const std::string &type) {
  return type == "HH" || type == "HT" || type == "EU";
}

bool isSubstanceGroupIdFree(const ROMol &mol, unsigned int id) {
  auto match_sgroup = [id](const SubstanceGroup &sg) {
    unsigned int storedId;
    return sg.getPropIfPresent("ID", storedId) && id == storedId;
  };

  const auto &sgroups = getSubstanceGroups(mol);
  return std::find_if(sgroups.begin(), sgroups.end(), match_sgroup) ==
         sgroups.end();
}
}  // namespace SubstanceGroupChecks

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
bool includesBond(SubstanceGroup &sg, unsigned int idx) {
  return sg.includesBond(idx);
}
bool includesAtom(SubstanceGroup &sg, unsigned int idx) {
  return sg.includesAtom(idx);
}
void removedBond(SubstanceGroup &sg, unsigned int idx) {
  sg.adjustToRemovedBond(idx);
}
void removedAtom(SubstanceGroup &sg, unsigned int idx) {
  sg.adjustToRemovedAtom(idx);
}

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

template <bool INCLUDES_METHOD(SubstanceGroup &, unsigned int),
          void ADJUST_METHOD(SubstanceGroup &, unsigned int)>
void removeSubstanceGroupsReferencing(RWMol &mol, unsigned int idx) {
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
      if (INCLUDES_METHOD(sgs[i], idx)) {
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
          ADJUST_METHOD(sg, idx);
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
              ADJUST_METHOD(sg, idx);
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
  removeSubstanceGroupsReferencing<includesAtom, removedAtom>(mol, idx);
}

void removeSubstanceGroupsReferencingBond(RWMol &mol, unsigned int idx) {
  // Delete substance groups containing this bond. It could be that it's ok to
  // keep it, but we just don't know
  removeSubstanceGroupsReferencing<includesBond, removedBond>(mol, idx);
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
