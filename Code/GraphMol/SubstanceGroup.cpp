//
//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstanceGroup.h"
#include "ROMol.h"

namespace RDKit {

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

bool SubstanceGroupChecks::isSubstanceGroupIdFree(const ROMol &mol,
                                                  unsigned int id) {
  auto match_sgroup = [&](const SubstanceGroup &sg) {
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

}  // namespace RDKit

std::ostream &operator<<(std::ostream &target,
                         const RDKit::SubstanceGroup &sgroup) {
  target << sgroup.getIndexInMol() << ' '
         << sgroup.getProp<std::string>("TYPE");
  return target;
}
