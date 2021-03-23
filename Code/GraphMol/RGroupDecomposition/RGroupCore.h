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

  // Bitset: indices corresponding to atoms bearing user-defined labels are 1
  boost::dynamic_bitset<> core_atoms_with_user_labels;
  // Number of user labelled rgroups in the core
  size_t numberUserRGroups = 0;
  RCore(){};
  RCore(const RWMol &c) : core(new RWMol(c)) {
    findIndicesWithRLabel();
    countUserRGroups();
  }

  inline bool isCoreAtomUserLabelled(int idx) const {
    return core_atoms_with_user_labels.test(idx);
  }

  void countUserRGroups() {
    numberUserRGroups = core_atoms_with_user_labels.count();
  }

  void findIndicesWithRLabel() {
    // Find all the core atoms that have user
    // label and set their indices to 1 into core_atoms_with_user_labels
    core_atoms_with_user_labels.resize(core->getNumAtoms());
    for (const auto atom : core->atoms()) {
      int label;
      if (atom->getPropIfPresent(RLABEL, label) && label > 0) {
        core_atoms_with_user_labels.set(atom->getIdx());
      }
    }
  }

  // Return a copy of core where dummy atoms are replaced by
  // the respective matching atom in mol, while other atoms have
  // their aromatic flag and formal charge copied from from
  // the respective matching atom in mol
  ROMOL_SPTR replaceCoreAtomsWithMolMatches(bool &hasCoreDummies,
                                            const ROMol &mol,
                                            const MatchVectType &match) const {
    auto coreReplacedAtoms = boost::make_shared<RWMol>(*core);
    hasCoreDummies = false;
    for (const auto &p : match) {
      auto atom = coreReplacedAtoms->getAtomWithIdx(p.first);
      if (atom->getAtomicNum() == 0) {
        hasCoreDummies = true;
      }
      if (isAtomWithMultipleNeighborsOrNotUserRLabel(*atom)) {
        auto molAtom = mol.getAtomWithIdx(p.second);
        replaceCoreAtom(*coreReplacedAtoms, *atom, *molAtom);
      }
    }

    for (auto bond : coreReplacedAtoms->bonds()) {
      if (bond->hasQuery()) {
        hasCoreDummies = true;
        const auto molBond =
            mol.getBondBetweenAtoms(match[bond->getBeginAtomIdx()].second,
                                    match[bond->getEndAtomIdx()].second);
        Bond newBond(molBond->getBondType());
        newBond.setIsAromatic(molBond->getIsAromatic());
        coreReplacedAtoms->replaceBond(bond->getIdx(), &newBond, true);
      }
    }

#ifdef VERBOSE
    std::cerr << "Original core smarts  " << MolToSmarts(*core) << std::endl;
    std::cerr << "Dummy replaced core smarts  "
              << MolToSmarts(*coreReplacedAtoms) << std::endl;
#endif
    return coreReplacedAtoms;
  }

  void replaceCoreAtom(RWMol &mol, Atom &atom, const Atom &other) const {
    auto atomicNumber = other.getAtomicNum();
    auto targetAtom = &atom;
    bool wasDummy = (atom.getAtomicNum() == 0);
    if (wasDummy) {
      if (atom.hasQuery()) {
        Atom newAtom(atomicNumber);
        auto atomIdx = atom.getIdx();
        mol.replaceAtom(atomIdx, &newAtom, false, true);
        targetAtom = mol.getAtomWithIdx(atomIdx);
      } else {
        atom.setAtomicNum(atomicNumber);
      }
    }
    targetAtom->setIsAromatic(other.getIsAromatic());
    targetAtom->setFormalCharge(other.getFormalCharge());
    if (wasDummy) {
      targetAtom->setNoImplicit(true);
      unsigned int numHs = 0;
      const auto &otherMol = other.getOwningMol();
      for (const auto &nbri :
           boost::make_iterator_range(otherMol.getAtomNeighbors(&other))) {
        const auto nbrAtom = otherMol[nbri];
        if (nbrAtom->getAtomicNum() == 1) {
          ++numHs;
        }
      }
      targetAtom->setNumExplicitHs(numHs + other.getTotalNumHs());
      targetAtom->updatePropertyCache(false);
    }
  }

  // Final core returned to user with dummy atoms and bonds set to those in the
  // match
  RWMOL_SPTR coreWithMatches(const ROMol &coreReplacedAtoms) const {
    auto finalCore = boost::make_shared<RWMol>(*labelledCore);
    for (size_t atomIdx = 0; atomIdx < coreReplacedAtoms.getNumAtoms();
         ++atomIdx) {
      auto coreAtom = finalCore->getAtomWithIdx(atomIdx);
      auto templateAtom = coreReplacedAtoms.getAtomWithIdx(atomIdx);
      auto unlabelledCoreAtom = core->getAtomWithIdx(atomIdx);
      if (templateAtom->getAtomicNum() > 0 &&
          isAtomWithMultipleNeighborsOrNotUserRLabel(*unlabelledCoreAtom)) {
        replaceCoreAtom(*finalCore, *coreAtom, *templateAtom);
      }
    }

    for (size_t bondIdx = 0; bondIdx < coreReplacedAtoms.getNumBonds();
         ++bondIdx) {
      auto coreBond = finalCore->getBondWithIdx(bondIdx);
      if (coreBond->hasQuery()) {
        auto templateBond = coreReplacedAtoms.getBondWithIdx(bondIdx);
        Bond newBond(templateBond->getBondType());
        newBond.setIsAromatic(templateBond->getIsAromatic());
        finalCore->replaceBond(bondIdx, &newBond, true);
      }
    }

    finalCore->updatePropertyCache(false);
    return finalCore;
  }
};

}  // namespace RDKit
#endif
