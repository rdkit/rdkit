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

namespace RDKit {

//! RCore is the core common to a series of molecules
struct RCore {
  boost::shared_ptr<RWMol> core;
  boost::shared_ptr<RWMol> labelledCore;
  std::set<int> core_atoms_with_user_labels;
  RCore(){};
  RCore(const RWMol &c, bool onlyMatchAtRGroups = false) : core(new RWMol(c)) {
    if (onlyMatchAtRGroups) {
      findIndicesWithRLabel();
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
  std::unique_ptr<ROMol> replaceCoreDummiesWithMolMatchesOriginal(
      const ROMol &mol, const MatchVectType &match) const {
    std::unique_ptr<ROMol> coreReplacedDummies(new ROMol(*core));
    for (const auto &p : match) {
      auto a = coreReplacedDummies->getAtomWithIdx(p.first);
      if (!a->getAtomicNum() && !a->hasProp(RLABEL)) {
        a->setAtomicNum(mol.getAtomWithIdx(p.second)->getAtomicNum());
      }
    }
    return coreReplacedDummies;
  }

  // Return a copy of core where dummy atoms are replaced by
  // the atomic number of the respective matching atom in mol
  std::shared_ptr<ROMol> replaceCoreDummiesWithMolMatches(
      const ROMol &mol, const MatchVectType &match) const {
    auto coreReplacedDummies = std::make_shared<RWMol>(*core);
    for (const auto &p : match) {
      auto atom = coreReplacedDummies->getAtomWithIdx(p.first);
      if (isIndexAnyRLabelOrMultipleConnectedUserRlabel(*atom)) {
        auto molAtom = mol.getAtomWithIdx(p.second);
        if (atom->hasQuery()) {
          auto atomicNumber = molAtom->getAtomicNum();
          auto newAtom = new Atom(atomicNumber);
          coreReplacedDummies->replaceAtom(p.first, newAtom, false, true);
          delete newAtom;
          atom = coreReplacedDummies->getAtomWithIdx(p.first);
        } else {
          atom->setAtomicNum(mol.getAtomWithIdx(p.second)->getAtomicNum());
        }
        atom->setIsAromatic(molAtom->getIsAromatic());
      }
    }

    for (auto bond : coreReplacedDummies->bonds()) {
      if (bond->hasQuery()) {
        const auto molBond =
            mol.getBondBetweenAtoms(match[bond->getBeginAtomIdx()].second,
                                    match[bond->getEndAtomIdx()].second);
        auto newBond = new Bond(molBond->getBondType());
        coreReplacedDummies->replaceBond(bond->getIdx(), newBond, true);
        delete newBond;
      }
    }

    std::cerr << "Original core smarts  " << MolToSmarts(*core) << std::endl;
    std::cerr << "Dummy replaced core smarts  "
              << MolToSmarts(*coreReplacedDummies) << std::endl;
    return coreReplacedDummies;
  }
};

}  // namespace RDKit
#endif
