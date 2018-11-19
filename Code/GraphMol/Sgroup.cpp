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
#include "Sgroup.h"
#include "ROMol.h"

namespace RDKit {

static const std::vector<std::string> sGroupTypes = {
    // polymer sgroups:
    "SRU", "MON", "COP", "CRO", "GRA", "MOD", "MER", "ANY",
    // formulations/mixtures:
    "COM", "MIX", "FOR",
    // other
    "SUP", "MUL", "DAT", "GEN"};

static const std::vector<std::string> sGroupSubtypes = {"ALT", "RAN", "BLO"};
static const std::vector<std::string> sGroupConnectTypes = {"HH", "HT", "EU"};

SGroup::SGroup(ROMol *owning_mol, const std::string &_type)
    : d_id(0), d_type(_type), dp_mol(owning_mol) {
  PRECONDITION(owning_mol, "owning molecule is nullptr");
}

void SGroup::setOwningMol(ROMol *mol) {
  PRECONDITION(mol, "owning molecule is nullptr");

  dp_mol = mol;
}

const std::string &SGroup::getProp(const std::string &prop) const {
  try {
    return d_prop.at(prop);
  } catch (const std::out_of_range &) {
    std::ostringstream errout;
    errout << "String Property '" << prop << "' not found in SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
}

//! check if SGroup has the given property set
bool SGroup::hasProp(const std::string &prop) const {
  return d_prop.find(prop) != d_prop.end();
}

namespace {
  bool isSGroupIdFree(const ROMol &mol,unsigned int id){
    auto sgroups = getMolSGroups(mol);
    CHECK_INVARIANT(sgroups,"no sgroups on molecule owning an sgroup");
	auto match_sgroup = [&](const SGROUP_SPTR &sg) { return id == sg->getId(); };
    return std::find_if(sgroups->begin(),sgroups->end(),match_sgroup) == sgroups->end();
  }

unsigned int getNextFreeSGroupId(const ROMol &mol) {
  std::set<unsigned int> ids;
  CHECK_INVARIANT(getMolSGroups(mol),"no sgroups on molecule");
  for (const auto &sgroup : *getMolSGroups(mol)) {
    ids.insert(sgroup->getId());
  }
  unsigned int nexId = 1;  // smallest possible ID
  for (const auto &id : ids) {
    if (id > nexId) {
      break;
    }
    ++nexId;
  }
  return nexId;
}

}

void SGroup::setId(unsigned int id) {
  PRECONDITION(dp_mol, "SGroup is not owned by any molecule");

  if (id == 0) {
    d_id = getNextFreeSGroupId(*dp_mol);
    return;
  } else if (!isSGroupIdFree(*dp_mol,id)) {
    std::ostringstream errout;
    errout << "ID " << id
           << " is already assigned to a SGroup on the same molecule";
    throw SGroupException(errout.str());
  }

  d_id = id;
}

unsigned int SGroup::getIndexInMol() const {
  if (!dp_mol) {
    std::ostringstream errout;
    errout << "SGroup has no owning molecule" << std::endl;
    throw SGroupException(errout.str());
  }

  auto match_sgroup = [&](const SGROUP_SPTR &sg) { return this == sg.get(); };
  auto sgroups = getMolSGroups(*dp_mol);
  TEST_ASSERT(sgroups);
	  auto sgroupItr =
      std::find_if(sgroups->begin(), sgroups->end(), match_sgroup);

  if (sgroupItr == sgroups->end()) {
    std::ostringstream errout;
    errout << "Unable to find own index in owning mol SGroup collection"
           << std::endl;
    throw SGroupException(errout.str());
  }

  return sgroupItr - sgroups->begin();
}

void SGroup::addAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");

  Atom *atom = dp_mol->getAtomWithIdx(idx);
  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot find Atom " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
  d_atoms.push_back(atom);
}

void SGroup::addAtomWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");
  Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
  d_atoms.push_back(atom);
}

void SGroup::addPAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_atoms.empty(), "no atoms");

  Atom *atom = dp_mol->getAtomWithIdx(idx);
  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot find Atom " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }

  if (std::find(d_atoms.begin(), d_atoms.end(), atom) == d_atoms.end()) {
    std::ostringstream errout;
    errout << "Atom " << idx << " is not a member of SGroup " << getId();
    throw SGroupException(errout.str());
  }

  d_patoms.push_back(atom);
}

void SGroup::addPAtomWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_atoms.empty(), "no atoms");

  Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
  if (std::find(d_atoms.begin(), d_atoms.end(), atom) == d_atoms.end()) {
    std::ostringstream errout;
    errout << "Atom with bookmark " << mark << " is not a member of SGroup "
           << getId();
    throw SGroupException(errout.str());
  }

  d_patoms.push_back(atom);
}

void SGroup::addBondWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");

  Bond *bond = dp_mol->getBondWithIdx(idx);
  if (!bond) {
    std::ostringstream errout;
    errout << "Cannot find Bond " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
  d_bonds.push_back(bond);
}

void SGroup::addBondWithBookmark(int mark) {
  PRECONDITION(dp_mol, "bad mol");
  Bond *bond = dp_mol->getUniqueBondWithBookmark(mark);
  d_bonds.push_back(bond);
}

void SGroup::addBracket(const SGroup::Bracket &bracket) {
  d_brackets.push_back(bracket);
}

void SGroup::addCState(Bond *bond, RDGeom::Point3D *vectorPtr) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_bonds.empty(), "no bonds");

  if (getBondType(bond) != SGroup::BondType::XBOND) {
    std::ostringstream errout;
    errout << "Bond with index " << bond->getIdx()
           << " is not an XBOND for SGroup " << getId();
    throw SGroupException(errout.str());
  }
  auto vectorSPtr = boost::shared_ptr<RDGeom::Point3D>(vectorPtr);
  d_cstates.push_back({bond, vectorSPtr});
}

void SGroup::addDataField(const std::string &data) {
  d_dataFields.push_back(data);
}

void SGroup::updateOwningMol(ROMol *other_mol) {
  for (auto &&atom : d_atoms) {
    unsigned int idx = atom->getIdx();
    atom = other_mol->getAtomWithIdx(idx);
  }

  for (auto &&patom : d_patoms) {
    unsigned int idx = patom->getIdx();
    patom = other_mol->getAtomWithIdx(idx);
  }

  for (auto &&bond : d_bonds) {
    unsigned int idx = bond->getIdx();
    bond = other_mol->getBondWithIdx(idx);
  }

  for (auto &&cstate : d_cstates) {
    unsigned int idx = cstate.bond->getIdx();
    cstate.bond = other_mol->getBondWithIdx(idx);
  }

  for (auto &&sap : d_saps) {
    unsigned int aIdx = sap.aAtom->getIdx();
    sap.aAtom = other_mol->getAtomWithIdx(aIdx);

    if (sap.lvAtom) {
      unsigned int lvIdx = sap.lvAtom->getIdx();
      sap.lvAtom = other_mol->getAtomWithIdx(lvIdx);
    }
  }

  if (d_parent) {
    auto parent_idx = d_parent->getIndexInMol();
    d_parent = getMolSGroup(*other_mol,parent_idx);
  }
}

void SGroup::addAttachPoint(Atom *aAtomPtr, Atom *lvAtomPtr,
                            std::string idStr) {
  d_saps.push_back({aAtomPtr, lvAtomPtr, idStr});
}

//! check if the bond is SGroup XBOND or CBOND
SGroup::BondType SGroup::getBondType(Bond *bond) const {
  PRECONDITION(std::find(d_bonds.begin(), d_bonds.end(), bond) != d_bonds.end(),
               "bond is not part of the SGroup")
  Atom *atom1 = bond->getBeginAtom();
  Atom *atom2 = bond->getEndAtom();

  bool atom1_in_sgroup =
      std::find(d_atoms.begin(), d_atoms.end(), atom1) != d_atoms.end();
  bool atom2_in_sgroup =
      std::find(d_atoms.begin(), d_atoms.end(), atom2) != d_atoms.end();

  if (atom1_in_sgroup && atom2_in_sgroup) {
    return SGroup::BondType::CBOND;
  } else if (atom1_in_sgroup || atom2_in_sgroup) {
    return SGroup::BondType::XBOND;
  } else {
    std::ostringstream errout;
    errout << "Neither beginning nor ending atoms of bond " << bond->getIdx()
           << " is in this SGroup.";
    throw SGroupException(errout.str());
  }
}

bool SGroupTypeOK(std::string typ) {
  return std::find(sGroupTypes.begin(), sGroupTypes.end(), typ) !=
         sGroupTypes.end();
}

bool SGroupSubTypeOK(std::string typ) {
  return std::find(sGroupSubtypes.begin(), sGroupSubtypes.end(), typ) !=
         sGroupSubtypes.end();
}

bool SGroupConnectTypeOK(std::string typ) {
  return std::find(sGroupConnectTypes.begin(), sGroupConnectTypes.end(), typ) !=
         sGroupConnectTypes.end();
}

unsigned int getMolNumSGroups(const ROMol &mol) {
  if(!getMolSGroups(mol)) return 0;
  return rdcast<unsigned int>(getMolSGroups(mol)->size());
}

std::vector<boost::shared_ptr<SGroup>> *getMolSGroups(ROMol &mol) {
	return &mol.d_sgroups;
}
const std::vector<boost::shared_ptr<SGroup>> *getMolSGroups(const ROMol &mol) {
	return &mol.d_sgroups;
}

SGroup *getMolSGroup(ROMol &mol,unsigned int idx) {
  // make sure we have more than one sgroup
  auto sgroups = getMolSGroups(mol);
  if (!sgroups || sgroups->empty()){
    throw SGroupException("No SGroups available on the molecule");
  }
  return sgroups->at(idx).get();
}
const SGroup *getMolSGroup(const ROMol &mol,unsigned int idx) {
  // make sure we have more than one sgroup
  auto sgroups = getMolSGroups(mol);
  if (!sgroups || sgroups->empty()){
    throw SGroupException("No SGroups available on the molecule");
  }
  return sgroups->at(idx).get();
}

unsigned int addMolSGroup(ROMol &mol, SGroup *sgroup) {
  sgroup->setOwningMol(&mol);
  SGROUP_SPTR nSGroup(sgroup);
  auto sgroups = getMolSGroups(mol);
  CHECK_INVARIANT(sgroups,"no s group container for molecule");
  unsigned int id = sgroups->size();
  sgroups->push_back(nSGroup);
  return id;
}


}  // namespace RDKit

// TO DO: finish this!
std::ostream &operator<<(std::ostream &target, const RDKit::SGroup &sgroup) {
  target << sgroup.getId() << " " << sgroup.getType();
  return target;
}
