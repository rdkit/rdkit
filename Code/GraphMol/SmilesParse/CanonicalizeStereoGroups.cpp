
//  Copyright (C) 2001-2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/new_canon.h>

#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/utils.h>
#include <algorithm>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>

namespace RDKit {

namespace {

void buildTree(int atomIndexToAdd, const ROMol *mol,
               std::vector<unsigned int> &chosenOrder,
               std::vector<int> &reverseOrder,
               std::vector<unsigned int> &ranks) {
  // build a tree of the atoms in the molecule
  // starting with the atom at atomIndexToAdd
  // and using the ranks to determine
  // the order of the neighbors in each atom in the tree.
  //
  // The chosenOrder is the list of old atom numbers in the order they are
  // chosen reverseOrder is the reference for each old atom number to its new
  // place in the new chosenOrder.
  //
  // the tree is built by adding, recursively, the neighbor atoms in order of
  // rank
  //

  PRECONDITION(mol, "bad molecule");

  chosenOrder.push_back(atomIndexToAdd);
  reverseOrder[atomIndexToAdd] = chosenOrder.size() - 1;

  auto atomToAdd = mol->getAtomWithIdx(atomIndexToAdd);

  std::vector<std::pair<unsigned int, unsigned int>> nbrRanks;
  nbrRanks.reserve(mol->getAtomDegree(atomToAdd));

  for (const auto nbr : mol->atomNeighbors(atomToAdd)) {
    nbrRanks.push_back(std::make_pair(ranks[nbr->getIdx()], nbr->getIdx()));
  }
  std::sort(nbrRanks.begin(), nbrRanks.end());
  for (const auto &pr : nbrRanks) {
    if (reverseOrder[pr.second] == -1) {
      buildTree(pr.second, mol, chosenOrder, reverseOrder, ranks);
    }
  }
}

class ChiralAtomItem {
 private:
  unsigned int atomId;
  RDKit::Atom::ChiralType chiralType;

 public:
  ChiralAtomItem() = delete;
  ChiralAtomItem(const RDKit::Atom *atomInit,
                 const std::vector<unsigned int> &atomsToInvert)
      : atomId(atomInit->getIdx()), chiralType(atomInit->getChiralTag()) {
    if (std::find(atomsToInvert.begin(), atomsToInvert.end(), atomId) !=
        atomsToInvert.end()) {
      if (chiralType == RDKit::Atom::CHI_TETRAHEDRAL_CW) {
        chiralType = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
      } else if (chiralType == RDKit::Atom::CHI_TETRAHEDRAL_CCW) {
        chiralType = RDKit::Atom::CHI_TETRAHEDRAL_CW;
      }
    }
  }

  unsigned int getAtomId() const { return atomId; }
  RDKit::Atom::ChiralType getChiralType() const { return chiralType; }

  bool operator<(const ChiralAtomItem &other) const {
    if (atomId < other.atomId) {
      return true;
    } else if (atomId > other.atomId) {
      return false;
    }
    // note:  CCW is considered less that CW
    if (chiralType < other.chiralType) {
      return true;
    } else if (chiralType > other.chiralType) {
      return false;
    }
    return false;
  }

  bool operator==(const ChiralAtomItem &other) const {
    if (atomId != other.atomId) {
      return false;
    }

    if (chiralType != other.chiralType) {
      return false;
    }
    return true;
  }

  bool operator!=(const ChiralAtomItem &other) const {
    return !(*this == other);
  }
};

class ChiralBondItem {
 private:
  RDKit::Bond::BondStereo stereoType = RDKit::Bond::BondStereo::STEREONONE;
  unsigned int bondId;
  unsigned int atomId1;
  unsigned int atomId2;

 public:
  unsigned int getBondId() const { return bondId; }
  unsigned int getAtomId1() const { return atomId1; }
  unsigned int getAtomId2() const { return atomId2; }

  RDKit::Bond::BondStereo getStereoType() const { return stereoType; }

  ChiralBondItem() = delete;
  ChiralBondItem(const RDKit::Bond *bondInit)
      : stereoType(bondInit->getStereo()),
        bondId(bondInit->getIdx()),
        atomId1(bondInit->getBeginAtomIdx()),
        atomId2(bondInit->getEndAtomIdx()) {
    if (atomId1 > atomId2) {
      std::swap(atomId1, atomId2);
    }
  }

  bool operator<(const ChiralBondItem &other) const {
    if (atomId1 < other.atomId1) {
      return true;
    } else if (atomId1 > other.atomId1) {
      return false;
    } else if (atomId2 < other.atomId2) {
      return true;
    } else if (atomId2 > other.atomId2) {
      return false;
    }

    if (stereoType < other.stereoType) {
      return true;
    } else if (stereoType > other.stereoType) {
      return false;
    }

    return false;
  }

  bool operator==(const ChiralBondItem &other) const {
    if (atomId1 != other.atomId1 || atomId2 != other.atomId2 ||
        stereoType != other.stereoType) {
      return false;
    }

    return true;
  }

  bool operator!=(const ChiralBondItem &other) const {
    return !(*this == other);
  }
};

class RankedValue {
 private:
  std::vector<ChiralAtomItem> chiralAtomItems;
  mutable std::vector<ChiralBondItem> chiralBondItems;
  mutable bool bondsSorted = false;

 public:
  void AddAtom(RDKit::Atom *atom,
               const std::vector<unsigned int> &atomsToInvert) {
    chiralAtomItems.emplace_back(atom, atomsToInvert);
  }

  void AddBond(RDKit::Bond *bond) {
    chiralBondItems.emplace_back(bond);
    bondsSorted = false;
  }

  unsigned int getNumChiralAtoms() const { return chiralAtomItems.size(); }
  unsigned int getNumChiralBonds() const { return chiralBondItems.size(); }

  const std::vector<ChiralAtomItem> &getChiralAtoms() const {
    return chiralAtomItems;
  }

  const std::vector<ChiralBondItem> &getChiralBonds() const {
    if (!bondsSorted) {
      if (chiralBondItems.size() > 1) {
        std::sort(chiralBondItems.begin(), chiralBondItems.end());
      }
      bondsSorted = true;
    }

    return chiralBondItems;
  }

  bool operator<(const RankedValue &other) const {
    if (chiralAtomItems.size() < other.chiralAtomItems.size()) {
      return true;
    } else if (chiralAtomItems.size() > other.chiralAtomItems.size()) {
      return false;
    }

    if (chiralBondItems.size() < other.chiralBondItems.size()) {
      return true;
    } else if (chiralBondItems.size() > other.chiralBondItems.size()) {
      return false;
    }

    for (auto it = chiralAtomItems.begin(), it2 = other.chiralAtomItems.begin();
         it != chiralAtomItems.end(); ++it, ++it2) {
      if (*it < *it2) {
        return true;
      } else if (*it2 < *it) {
        return false;
      }
    }

    for (auto it = chiralBondItems.begin(), it2 = other.chiralBondItems.begin();
         it != chiralBondItems.end(); ++it, ++it2) {
      if (*it < *it2) {
        return true;
      } else if (*it2 < *it) {
        return false;
      }
    }

    return false;
  }

  bool equivalentTo(const RankedValue &other) const {
    if (chiralAtomItems.size() != other.chiralAtomItems.size()) {
      return false;
    }

    if (chiralBondItems.size() != other.chiralBondItems.size()) {
      return false;
    }

    for (auto it = chiralAtomItems.begin(), it2 = other.chiralAtomItems.begin();
         it != chiralAtomItems.end(); ++it, ++it2) {
      if ((*it).getAtomId() != (*it2).getAtomId()) {
        return false;
      }
    }

    for (auto it = chiralBondItems.begin(), it2 = other.chiralBondItems.begin();
         it != chiralBondItems.end(); ++it, ++it2) {
      if ((*it).getAtomId1() != (*it2).getAtomId1() ||
          (*it).getAtomId2() != (*it2).getAtomId2()) {
        return false;
      }
    }

    return true;
  }

  bool operator==(const RankedValue &other) const {
    if (chiralAtomItems.size() != other.chiralAtomItems.size()) {
      return false;
    }

    if (chiralBondItems.size() != other.chiralBondItems.size()) {
      return false;
    }

    for (auto it = chiralAtomItems.begin(), it2 = other.chiralAtomItems.begin();
         it != chiralAtomItems.end(); ++it, ++it2) {
      if (*it != *it2) {
        return false;
      }
    }

    for (auto it = chiralBondItems.begin(), it2 = other.chiralBondItems.begin();
         it != chiralBondItems.end(); ++it, ++it2) {
      if (*it != *it2) {
        return false;
      }
    }

    return true;
  }
};

bool doesAtomChiralityVary(const std::set<RankedValue> &allRankedValues,
                           unsigned int index) {
  PRECONDITION(!allRankedValues.empty(), "bad allRankedValues size");
  PRECONDITION(allRankedValues.begin()->getNumChiralAtoms() > index,
               "index out of range");

  const auto firstChiralVal =
      allRankedValues.begin()->getChiralAtoms()[index].getChiralType();

  for (const auto &rankedValue : allRankedValues) {
    if (rankedValue.getChiralAtoms()[index].getChiralType() != firstChiralVal) {
      return true;
    }
  }

  return false;
}

bool doesBondStereoVary(const std::set<RankedValue> &allRankedValues,
                        unsigned int index) {
  PRECONDITION(!allRankedValues.empty(), "bad allRankedValues size");
  PRECONDITION(allRankedValues.begin()->getNumChiralBonds() > index,
               "index out of range");

  const auto firstStereoVal =
      allRankedValues.begin()->getChiralBonds()[index].getStereoType();

  for (const auto &rankedValue : allRankedValues) {
    if (rankedValue.getChiralBonds()[index].getStereoType() != firstStereoVal) {
      return true;
    }
  }

  return false;
}

bool doTwoAtomsVaryTheSame(const std::set<RankedValue> &allRankedValues,
                           unsigned int index1, unsigned int index2) {
  PRECONDITION(!allRankedValues.empty(), "bad allRankedValues size");
  PRECONDITION(allRankedValues.begin()->getNumChiralAtoms() > index1,
               "index1 out of range");
  PRECONDITION(allRankedValues.begin()->getNumChiralAtoms() > index2,
               "index2 out of range");

  const auto firstChiralVal1 =
      allRankedValues.begin()->getChiralAtoms()[index1].getChiralType();
  const auto firstChiralVal2 =
      allRankedValues.begin()->getChiralAtoms()[index2].getChiralType();

  for (const auto &rankedValue : allRankedValues) {
    if ((firstChiralVal1 ==
         rankedValue.getChiralAtoms()[index1].getChiralType()) !=
        (firstChiralVal2 ==
         rankedValue.getChiralAtoms()[index2].getChiralType())) {
      return false;
    }
  }

  return true;
}

bool doAtomAndBondVaryTheSame(const std::set<RankedValue> &allRankedValues,
                              unsigned int atomIndex1,
                              unsigned int bondIndex2) {
  PRECONDITION(!allRankedValues.empty(), "bad allMols size");
  PRECONDITION(allRankedValues.begin()->getNumChiralAtoms() > atomIndex1,
               "atomIndex1 out of range");
  PRECONDITION(allRankedValues.begin()->getNumChiralBonds() > bondIndex2,
               "bondIndex2 out of range");

  const auto firstChiralVal1 =
      allRankedValues.begin()->getChiralAtoms()[atomIndex1].getChiralType();
  const auto firstStereoVal2 =
      allRankedValues.begin()->getChiralBonds()[bondIndex2].getStereoType();

  for (const auto &rankedValue : allRankedValues) {
    if ((firstChiralVal1 ==
         rankedValue.getChiralAtoms()[atomIndex1].getChiralType()) !=
        (firstStereoVal2 ==
         rankedValue.getChiralBonds()[bondIndex2].getStereoType())) {
      return false;
    }
  }

  return true;
}

bool doTwoBondsVaryTheSame(const std::set<RankedValue> &allRankedValues,
                           unsigned int bondIndex1, unsigned int bondIndex2) {
  PRECONDITION(!allRankedValues.empty(), "bad allMols size");
  PRECONDITION(allRankedValues.begin()->getNumChiralBonds() > bondIndex1,
               "atomIndex1 out of range");
  PRECONDITION(allRankedValues.begin()->getNumChiralBonds() > bondIndex2,
               "bondIndex2 out of range");

  const auto firstStereoVal1 =
      allRankedValues.begin()->getChiralBonds()[bondIndex1].getStereoType();
  const auto firstStereoVal2 =
      allRankedValues.begin()->getChiralBonds()[bondIndex2].getStereoType();

  for (const auto &rankedValue : allRankedValues) {
    if ((firstStereoVal1 ==
         rankedValue.getChiralBonds()[bondIndex1].getStereoType()) !=
        (firstStereoVal2 ==
         rankedValue.getChiralBonds()[bondIndex2].getStereoType())) {
      return false;
    }
  }

  return true;
}

unsigned int countSwaps(std::vector<unsigned int> &nbrs) {
  unsigned int swaps = 0;
  for (unsigned int i = 0; i < nbrs.size(); ++i) {
    for (unsigned int j = i + 1; j < nbrs.size(); ++j) {
      if (nbrs[i] > nbrs[j]) {
        ++swaps;
      }
    }
  }
  return swaps;
}

}  // namespace

// the call to renumberAtoms will NOT invert chiral atoms.
// it does return the atoms in the order handed to it, which could be the order
// of atoms for a possible smiles string.
//
// RDKit internally bases the chiral atoms on the order of the bonds
// to an atom, and the renumber function below does not change the order of the
// bonds to an atom.  So, the chiral atoms are still correct and are unchanged
// by renumber.
//
//  this routine determines which atoms would be inverted in an actual
// smiles were it to be generated.  This allows processing of a mol in a
// smiles-like order without actually generating a smiles string.

void getAtomsToInvert2(const RDKit::ROMol &mol,
                       const std::vector<unsigned int> &newOrder,
                       const std::vector<int> &reversedOrder,
                       std::vector<unsigned int> &atomsToInvert) {
  unsigned int nAts = mol.getNumAtoms();
  PRECONDITION(newOrder.size() == nAts, "bad newOrder size");

  // copy over the atoms:
  for (unsigned int nIdx = 0; nIdx < nAts; ++nIdx) {
    unsigned int oIdx = newOrder[nIdx];
    const RDKit::Atom *oAtom = mol.getAtomWithIdx(oIdx);

    if (oAtom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED) {
      // get the neighbors in the new order

      std::vector<unsigned int> nbrs;
      nbrs.reserve(oAtom->getDegree());

      for (const auto &nbr : mol.atomNeighbors(oAtom)) {
        nbrs.push_back(reversedOrder[nbr->getIdx()]);
      }
      if (RDKit::countSwaps(nbrs) % 2) {
        atomsToInvert.push_back(nIdx);
      }
    }
  }

  return;
}

void addSingleAbsGroup(ROMol &mol) {
  // all chiral centers are added to an abs group
  // if there are not chiral centers, no group is added

  std::vector<StereoGroup> sgs;
  std::vector<Atom *> chiralAtoms;
  std::vector<Bond *> chiralBonds;
  for (auto &atom : mol.atoms()) {
    if (atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
        atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
      chiralAtoms.push_back(atom);
    }
  }
  for (auto &bond : mol.bonds()) {
    if (bond->getStereo() == Bond::BondStereo::STEREOATROPCW ||
        bond->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
      chiralBonds.push_back(bond);
    }
  }

  if (!chiralAtoms.empty() || !chiralBonds.empty()) {
    sgs.emplace_back(StereoGroupType::STEREO_ABSOLUTE, chiralAtoms,
                     chiralBonds);
  }
  mol.setStereoGroups(sgs);  // could be empty, or have one abs group
}

void clearStereoGroups(ROMol &mol) {
  // all chiral centers are added to an abs group
  // if there are not chiral centers, no group is added
  std::vector<StereoGroup> sgs;
  mol.setStereoGroups(sgs);
}

void canonicalizeStereoGroups_internal(
    std::unique_ptr<RDKit::ROMol> &mol, RDKit::StereoGroupType stereoGroupType,
    RDKit::StereoGroupAbsOptions outputAbsoluteGroups) {
  // this expands a mol with stereo groups to a vector of values that are the
  // result of expanding the stereo groups, then determines the stereo groups
  // from that set

  std::set<RDKit::RankedValue> allRankedValues;

  std::vector<RDKit::StereoGroup> groupsToProcess;
  std::vector<RDKit::StereoGroup> andGroupsToKeep;

  for (auto &grp : mol->getStereoGroups()) {
    if (stereoGroupType == grp.getGroupType()) {
      groupsToProcess.push_back(grp);
    } else if (stereoGroupType == RDKit::StereoGroupType::STEREO_OR &&
               grp.getGroupType() == RDKit::StereoGroupType::STEREO_AND) {
      andGroupsToKeep.push_back(grp);
    }
  }

  mol->setStereoGroups(
      andGroupsToKeep);  // these groups might be empty, especially if we
                         // are PROCESSING AND groups
  std::unique_ptr<RDKit::ROMol> bestNewMol;
  auto newMolCount = std::pow(2, groupsToProcess.size());

  for (unsigned int molIndex = 0; molIndex < newMolCount; ++molIndex) {
    auto newMol = std::unique_ptr<RDKit::ROMol>(new RDKit::RWMol(*(mol.get())));

    for (unsigned int grpIndex = 0; grpIndex < groupsToProcess.size();
         ++grpIndex) {
      if (molIndex & (1 << grpIndex)) {
        for (auto atomPtr : groupsToProcess[grpIndex].getAtoms()) {
          if (atomPtr->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW) {
            newMol->getAtomWithIdx(atomPtr->getIdx())
                ->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
          } else if (atomPtr->getChiralTag() ==
                     RDKit::Atom::CHI_TETRAHEDRAL_CCW) {
            newMol->getAtomWithIdx(atomPtr->getIdx())
                ->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
          }
        }
        // do any  atropisomer bonds in this stereo group

        for (auto bond : groupsToProcess[grpIndex].getBonds()) {
          if (bond->getStereo() == RDKit::Bond::STEREOATROPCW) {
            newMol->getBondWithIdx(bond->getIdx())
                ->setStereo(RDKit::Bond::STEREOATROPCCW);
          } else if (bond->getStereo() == RDKit::Bond::STEREOATROPCCW) {
            newMol->getBondWithIdx(bond->getIdx())
                ->setStereo(RDKit::Bond::STEREOATROPCW);
          }
        }
      }
    }

    if (!andGroupsToKeep.empty()) {
      canonicalizeStereoGroups_internal(
          newMol, RDKit::StereoGroupType::STEREO_AND,
          RDKit::StereoGroupAbsOptions::NeverInclude);
    }
    std::vector<unsigned int> ranks(mol->getNumAtoms());

    const bool breakTies = true;
    const bool includeChirality = true;
    const bool includeIsotopes = false;
    const bool includeAtomMaps = true;
    const bool useNonStereoRanks = true;
    const bool includeChiralPresence = true;
    const bool includeStereoGroups = true;

    RDKit::Canon::rankMolAtoms(*newMol, ranks, breakTies, includeChirality,
                               includeIsotopes, includeAtomMaps,
                               includeChiralPresence, includeStereoGroups,
                               useNonStereoRanks);

    // create an atoms ordering as if a smiles, but do this for the entire
    // mol - not fragments - it really is NOT the same order as generating a
    // smiles

    std::vector<unsigned int> chosenOrder;
    std::vector<int> reversedOrder(newMol->getNumAtoms(), -1);
    while (true) {
      int startingAtomIndex = -1;
      unsigned int lowestRank = UINT_MAX;

      for (unsigned int i = 0; i < newMol->getNumAtoms(); ++i) {
        if (reversedOrder[i] != -1) {
          continue;
        }
        if (ranks[i] < lowestRank) {
          lowestRank = ranks[i];
          startingAtomIndex = i;
        }
      }
      if (startingAtomIndex == -1) {
        break;  // all atoms are done
      }

      RDKit::buildTree(startingAtomIndex, newMol.get(), chosenOrder,
                       reversedOrder, ranks);
    }
    if (newMol->getNumAtoms() != chosenOrder.size()) {
      throw ValueErrorException("atomOrdering size mismatch");
    }

    // the call to renumberAtoms will NOT invert chiral atoms as would
    // happen if a smiles string were to be generated.
    //
    // the call to getAtomsToInvert will determine which atoms would be
    // inverted in a smiles string

    std::vector<unsigned int> atomsToInvert;
    RDKit::getAtomsToInvert2(*newMol.get(), chosenOrder, reversedOrder,
                             atomsToInvert);

    newMol.reset((RDKit::RWMol *)RDKit::MolOps::renumberAtoms(*newMol.get(),
                                                              chosenOrder));

    RDKit::RankedValue newRankedValue;

    boost::dynamic_bitset<> atomIndicesInStereoGroups(newMol->getNumAtoms());
    boost::dynamic_bitset<> bondIndicesInStereoGroups(newMol->getNumBonds());

    for (auto grp : newMol->getStereoGroups()) {
      for (auto atomPtr : grp.getAtoms()) {
        atomIndicesInStereoGroups.set(atomPtr->getIdx());
      }
      for (auto bondPtr : grp.getBonds()) {
        bondIndicesInStereoGroups.set(bondPtr->getIdx());
      }
    }

    // now get all chiral centers and atrop bonds that are not in the stereo
    // groups

    for (auto atom : newMol->atoms()) {
      if ((atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW ||
           atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW) &&
          !atomIndicesInStereoGroups[atom->getIdx()]) {
        newRankedValue.AddAtom(atom, atomsToInvert);
      }
    }

    for (auto bond : newMol->bonds()) {
      if ((bond->getStereo() == RDKit::Bond::BondStereo::STEREOATROPCCW ||
           bond->getStereo() == RDKit::Bond::BondStereo::STEREOATROPCW) &&
          !bondIndicesInStereoGroups[bond->getIdx()]) {
        newRankedValue.AddBond(bond);
      }
    }

    atomIndicesInStereoGroups.clear();  // not needed past here
    bondIndicesInStereoGroups.clear();  // not needed past here

    if (!allRankedValues.empty() &&
        !newRankedValue.equivalentTo(*allRankedValues.begin())) {
      throw RDKit::RigorousEnhancedStereoException(
          "ranked items  are not equivalent");
    }

    auto insertResult = allRankedValues.insert(newRankedValue);

    if (insertResult.second && newRankedValue == *allRankedValues.begin()) {
      bestNewMol = std::move(newMol);
    }
  }

  // now figure out the stereo groups to create

  std::vector<bool> atomsDone(allRankedValues.begin()->getNumChiralAtoms(),
                              false);
  std::vector<bool> bondsDone(allRankedValues.begin()->getNumChiralBonds(),
                              false);

  std::vector<RDKit::Atom *> absGroupAtoms;
  std::vector<RDKit::Bond *> absGroupBonds;

  // if there is only one smiles, then there is no variation and
  // add stereo groups are actual abs (and there is only one group)

  std::vector<RDKit::StereoGroup> newGroups;

  if (allRankedValues.size() == 1) {
    if (outputAbsoluteGroups == RDKit::StereoGroupAbsOptions::NeverInclude) {
      mol.swap(bestNewMol);
      return;
    }

    for (const auto &chiralAtom : allRankedValues.begin()->getChiralAtoms()) {
      absGroupAtoms.push_back(
          bestNewMol->getAtomWithIdx(chiralAtom.getAtomId()));
    }
    for (const auto &chiralBond : allRankedValues.begin()->getChiralBonds()) {
      absGroupBonds.push_back(
          bestNewMol->getBondWithIdx(chiralBond.getBondId()));
    }

  } else {
    // Now make the new stereo-enhanced mol

    unsigned int groupCount = 0;

    for (unsigned int index1 = 0;
         index1 < allRankedValues.begin()->getNumChiralAtoms(); ++index1) {
      if (atomsDone[index1]) {
        continue;
      }

      unsigned int atomIndex1 =
          allRankedValues.begin()->getChiralAtoms()[index1].getAtomId();

      if (!doesAtomChiralityVary(allRankedValues, index1)) {
        if (outputAbsoluteGroups !=
            RDKit::StereoGroupAbsOptions::NeverInclude) {
          absGroupAtoms.push_back(bestNewMol->getAtomWithIdx(atomIndex1));
        }
        atomsDone[index1] = true;
        continue;
      }

      std::vector<RDKit::Atom *> atomsToAdd;
      atomsToAdd.push_back(bestNewMol->getAtomWithIdx(atomIndex1));
      atomsDone[index1] = true;

      // now look through all other possible atoms and bonds to see if they
      // vary the same way as the first one in the group

      for (unsigned int index2 = index1 + 1;
           index2 < allRankedValues.begin()->getNumChiralAtoms(); ++index2) {
        if (atomsDone[index2]) {
          continue;
        }
        unsigned int atomIndex2 =
            allRankedValues.begin()->getChiralAtoms()[index2].getAtomId();

        if (doTwoAtomsVaryTheSame(allRankedValues, index1, index2)) {
          atomsToAdd.push_back(bestNewMol->getAtomWithIdx(atomIndex2));
          atomsDone[index2] = true;
        }
      }

      std::vector<RDKit::Bond *> bondsToAdd;
      for (unsigned int index2 = 0;
           index2 < allRankedValues.begin()->getNumChiralBonds(); ++index2) {
        if (bondsDone[index2]) {
          continue;
        }

        auto bondIndex2 =
            allRankedValues.begin()->getChiralBonds()[index2].getBondId();

        if (doAtomAndBondVaryTheSame(allRankedValues, index1, index2)) {
          bondsToAdd.push_back(bestNewMol->getBondWithIdx(bondIndex2));
          bondsDone[index2] = true;
        }
      }

      std::sort(atomsToAdd.begin(), atomsToAdd.end(),
                [](const RDKit::Atom *a, const RDKit::Atom *b) {
                  return a->getIdx() < b->getIdx();
                });
      std::sort(bondsToAdd.begin(), bondsToAdd.end(),
                [](const RDKit::Bond *a, const RDKit::Bond *b) {
                  return a->getIdx() < b->getIdx();
                });
      newGroups.emplace_back(stereoGroupType, atomsToAdd, bondsToAdd,
                             ++groupCount);
    }

    // now any groups that only involve bonds

    for (unsigned int index1 = 0;
         index1 < allRankedValues.begin()->getNumChiralBonds(); ++index1) {
      if (bondsDone[index1]) {
        continue;
      }

      unsigned int bondIndex1 =
          allRankedValues.begin()->getChiralBonds()[index1].getBondId();

      if (!doesBondStereoVary(allRankedValues, index1)) {
        if (outputAbsoluteGroups !=
            RDKit::StereoGroupAbsOptions::NeverInclude) {
          absGroupBonds.push_back(bestNewMol->getBondWithIdx(bondIndex1));
        }
        bondsDone[index1] = true;
        continue;
      }

      std::vector<RDKit::Bond *> bondsToAdd;
      bondsToAdd.push_back(bestNewMol->getBondWithIdx(bondIndex1));
      bondsDone[index1] = true;

      // now look through all other possible bonds to see if they vary
      // the same way as the first one in the group

      for (unsigned int index2 = index1 + 1;
           index2 < allRankedValues.begin()->getNumChiralBonds(); ++index2) {
        if (bondsDone[index2]) {
          continue;
        }
        unsigned int bondIndex2 =
            allRankedValues.begin()->getChiralBonds()[index2].getBondId();

        if (doTwoBondsVaryTheSame(allRankedValues, bondIndex1, bondIndex2)) {
          bondsToAdd.push_back(bestNewMol->getBondWithIdx(bondIndex2));
          bondsDone[index2] = true;
        }
      }

      std::sort(bondsToAdd.begin(), bondsToAdd.end(),
                [](const RDKit::Bond *a, const RDKit::Bond *b) {
                  return a->getIdx() < b->getIdx();
                });
      std::vector<RDKit::Atom *> atomsToAdd;  // nothing added to this one here
      newGroups.emplace_back(stereoGroupType, atomsToAdd, bondsToAdd,
                             ++groupCount);
    }
  }

  // keep the groups from the best mol, if it had them from a call to this
  // routine for the other kind of stereo groups.

  if (!bestNewMol->getStereoGroups().empty()) {
    for (auto grp : bestNewMol->getStereoGroups()) {
      newGroups.push_back(grp);
    }
  }

  // if the abs group is not empty, add it
  if ((outputAbsoluteGroups == RDKit::StereoGroupAbsOptions::AlwaysInclude ||
       (outputAbsoluteGroups ==
            RDKit::StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist &&
        !newGroups.empty())) &&
      (!absGroupAtoms.empty() || !absGroupBonds.empty())) {
    std::sort(absGroupAtoms.begin(), absGroupAtoms.end(),
              [](const RDKit::Atom *a, const RDKit::Atom *b) {
                return a->getIdx() < b->getIdx();
              });
    std::sort(absGroupBonds.begin(), absGroupBonds.end(),
              [](const RDKit::Bond *a, const RDKit::Bond *b) {
                return a->getIdx() < b->getIdx();
              });

    newGroups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                           absGroupAtoms, absGroupBonds, 0);
  }

  bestNewMol->setStereoGroups(newGroups);

  mol = std::unique_ptr<RDKit::ROMol>(bestNewMol.release());

  return;
}
void canonicalizeStereoGroups(std::unique_ptr<ROMol> &mol,
                              StereoGroupAbsOptions outputAbsoluteGroups) {
  // this returns a mol that has a caononical rep for the enhanced stereo
  // groups it expands the given mol to all possible non-stereo-group mols,
  // then determines a single set of stereo groups that uniquely represent
  // that group.

  // if there are both OR and AND groups, the AND groups are done first
  // by haveing the working routine call itself for pre-prosing the AND
  // groups

  auto sgCount = mol->getStereoGroups().size();
  if (sgCount == 0 ||
      (sgCount == 1 && mol->getStereoGroups()[0].getGroupType() ==
                           StereoGroupType::STEREO_ABSOLUTE)) {
    if (outputAbsoluteGroups == StereoGroupAbsOptions::AlwaysInclude) {
      addSingleAbsGroup(*mol);
    } else {
      clearStereoGroups(*mol);
    }
    return;
  }

  //  see if it is a simple compound, which has only one stereo group (not
  //  abs)
  // and if has two or fewer atoms in the group, and all chiral atoms are in
  // that group.
  //
  // furthermore, if the findMeso routine matches the simple group, it
  // should be removed
  //

  if (sgCount == 1) {
    const auto &sg = mol->getStereoGroups()[0];
    const auto &sgats = sg.getAtoms();
    const auto &sgBonds = sg.getBonds();
    if (sgats.size() <= 2 && sgBonds.size() == 0) {
      bool isSimple = true;

      for (auto &atom : mol->atoms()) {
        if ((atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
             atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW) &&
            std::find(sgats.begin(), sgats.end(), atom) == sgats.end()) {
          isSimple = false;
          break;
        }
      }
      if (isSimple) {
        for (auto &bond : mol->bonds()) {
          if ((bond->getStereo() == Bond::BondStereo::STEREOATROPCW ||
               bond->getStereo() == Bond::BondStereo::STEREOATROPCCW) &&
              std::find(sgBonds.begin(), sgBonds.end(), bond) ==
                  sgBonds.end()) {
            isSimple = false;
            break;
          }
        }

        if (isSimple) {
          auto res = Chirality::findMesoCenters(*mol);
          if (res.size() == 1) {
            if (outputAbsoluteGroups ==
                RDKit::StereoGroupAbsOptions::AlwaysInclude) {
              addSingleAbsGroup(*mol);
            } else {
              clearStereoGroups(*mol);
            }
          }

          return;  // we will not process the simple ones.   If the meso atoms
                   // were found, the one group was removed
        }
      }
    }
  }

  bool foundOrGroup = false;
  for (auto &stg : mol->getStereoGroups()) {
    if (stg.getGroupType() == StereoGroupType::STEREO_OR) {
      foundOrGroup = true;
      break;
    }
  }

  if (mol->needsUpdatePropertyCache()) {
    mol->updatePropertyCache(true);
  }

  // get the non-stereo rankings - these do NOT change as we iterate over
  // the enhanced possibilties.  They also do not change if the re-entrant
  // call is made

  std::vector<unsigned int> ranks(mol->getNumAtoms());

  const bool breakTies = false;
  const bool includeChirality = false;
  const bool includeIsotopes = false;
  const bool includeAtomMaps = true;
  const bool useNonStereoRanks = false;
  const bool includeChiralPresence = true;
  const bool includeStereoGroups = false;

  Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality, includeIsotopes,
                      includeAtomMaps, includeChiralPresence,
                      includeStereoGroups, useNonStereoRanks);

  for (auto atom : mol->atoms()) {
    atom->setProp(common_properties::_CanonicalRankingNumber,
                  ranks[atom->getIdx()]);
  }

  auto savedStereoGroups = mol->getStereoGroups();
  try {
    if (!foundOrGroup) {
      RDKit::canonicalizeStereoGroups_internal(mol, StereoGroupType::STEREO_AND,
                                               outputAbsoluteGroups);
    } else {
      RDKit::canonicalizeStereoGroups_internal(mol, StereoGroupType::STEREO_OR,
                                               outputAbsoluteGroups);
    }

    // Fix up the mol - round trip through smiles

    mol->clearComputedProps();

    SmilesWriteParams wp;
    wp.canonical = false;
    auto finalSmiles = MolToCXSmiles(*mol, wp);
    SmilesParserParams ps;
    ps.sanitize = false;
    mol.reset(SmilesToMol(finalSmiles, ps));

    return;
  } catch (const RigorousEnhancedStereoException &e) {
    mol->setStereoGroups(savedStereoGroups);
    return;
  }
}

}  // namespace RDKit
