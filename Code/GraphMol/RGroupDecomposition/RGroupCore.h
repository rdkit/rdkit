//
//  Copyright (C) 2020 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_CORE
#define RGROUP_CORE

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include "../RDKitBase.h"
#include "RGroupUtils.h"
#include "GraphMol/Substruct/SubstructMatch.h"

// #define VERBOSE 1

namespace RDKit {

//! RCore is the core common to a series of molecules
struct RCore {
  boost::shared_ptr<RWMol> core;
  boost::shared_ptr<RWMol> labelledCore;

  // A list of user indices for when onlyMatchAtRGroups = True
  std::set<int> core_atoms_with_user_labels;
  // Number of user labelled rgroups in the core
  size_t numberUserRGroups = 0;
  RCore(){};
  RCore(const RWMol &c, bool onlyMatchAtRGroups = false) : core(new RWMol(c)) {
    // Remove this from constructor if the create new core path can be removed from RGroupDecomposition::add
    if (onlyMatchAtRGroups) {
      findIndicesWithRLabel();
    }
    countUserRGroups();
  }

  void countUserRGroups() {
    numberUserRGroups = 0;
    for (const auto atom : core->atoms()) {
      int label;
      if (atom->getPropIfPresent(RLABEL, label)) {
        if (label > 0) {
          ++numberUserRGroups;
        }
      }
    }
  }

  void findIndicesWithRLabel() {
    // First find all the core atoms that have user
    //  label and put their indices into core_atoms_with_user_labels
    for (const auto atom : core->atoms()) {
      if (atom->hasProp(RLABEL)) {
        core_atoms_with_user_labels.insert(atom->getIdx());
      }
    }
  }

  // Return a copy of core where dummy atoms are replaced by
  // the atomic number of the respective matching atom in mol
  ROMOL_SPTR replaceCoreDummiesWithMolMatches(
      bool &hasDummies, const ROMol &mol, const MatchVectType &match) const {
    auto coreReplacedDummies = boost::make_shared<RWMol>(*core);
    hasDummies = false;
    for (const auto &p : match) {
      auto atom = coreReplacedDummies->getAtomWithIdx(p.first);
      if (isAnyAtomWithMultipleNeighborsOrNotUserRLabel(*atom)) {
        hasDummies = true;
        auto molAtom = mol.getAtomWithIdx(p.second);
        replaceDummyAtom(*coreReplacedDummies, *atom, *molAtom);
      }
    }

    for (auto bond : coreReplacedDummies->bonds()) {
      if (bond->hasQuery()) {
        hasDummies = true;
        const auto molBond =
            mol.getBondBetweenAtoms(match[bond->getBeginAtomIdx()].second,
                                    match[bond->getEndAtomIdx()].second);
        Bond newBond(molBond->getBondType());
        coreReplacedDummies->replaceBond(bond->getIdx(), &newBond, true);
      }
    }

#ifdef VERBOSE
    std::cerr << "Original core smarts  " << MolToSmarts(*core) << std::endl;
    std::cerr << "Dummy replaced core smarts  "
              << MolToSmarts(*coreReplacedDummies) << std::endl;
#endif

    return coreReplacedDummies;
  }

  void replaceDummyAtom(RWMol &mol, Atom &atom, const Atom &other) const {
    PRECONDITION(atom.getAtomicNum() == 0, "Atom must be dummy");
    auto atomicNumber = other.getAtomicNum();
    if (atom.hasQuery()) {
      Atom newAtom(atomicNumber);
      newAtom.setIsAromatic(other.getIsAromatic());
      mol.replaceAtom(atom.getIdx(), &newAtom, false, true);
    } else {
      atom.setAtomicNum(other.getAtomicNum());
      atom.setIsAromatic(other.getIsAromatic());
    }
  }

  // Final core returned to user with dummy atoms and bonds set to those in the
  // match
  ROMOL_SPTR coreWithMatches(const ROMol &coreReplacedDummies) {
    auto finalCore = boost::make_shared<RWMol>(*labelledCore);
    for (size_t atomIdx = 0; atomIdx < coreReplacedDummies.getNumAtoms();
         ++atomIdx) {
      auto coreAtom = finalCore->getAtomWithIdx(atomIdx);
      auto templateAtom = coreReplacedDummies.getAtomWithIdx(atomIdx);
      auto unlabelledCoreAtom = core->getAtomWithIdx(atomIdx);
      if (coreAtom->getAtomicNum() == 0 && templateAtom->getAtomicNum() > 0 &&
          isAnyAtomWithMultipleNeighborsOrNotUserRLabel(*unlabelledCoreAtom)) {
        replaceDummyAtom(*finalCore, *coreAtom, *templateAtom);
        finalCore->getAtomWithIdx(atomIdx)->setNoImplicit(true);
      }
    }

    for (size_t bondIdx = 0; bondIdx < coreReplacedDummies.getNumBonds();
         ++bondIdx) {
      auto coreBond = finalCore->getBondWithIdx(bondIdx);
      if (coreBond->hasQuery()) {
        auto templateBond = coreReplacedDummies.getBondWithIdx(bondIdx);
        Bond newBond(templateBond->getBondType());
        finalCore->replaceBond(bondIdx, &newBond, true);
      }
    }

    finalCore->updatePropertyCache(false);
    return finalCore;
  }
};

}  // namespace RDKit
#endif
