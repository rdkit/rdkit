//
//  Copyright (c) 2017-2023, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RGroupData.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <regex>

namespace RDKit {

void RGroupData::add(boost::shared_ptr<ROMol> newMol,
                     const std::vector<int> &rlabel_attachments) {
  // some fragments can be add multiple times if they are cyclic
  for (auto &mol : mols) {
    if (newMol.get() == mol.get()) {
      return;
    }
  }

  if (mols.size() > 0) {
    // don't add extraneous hydrogens
    if (isMolHydrogen(*newMol)) {
      return;
    }
    if (is_hydrogen) {
      // if we are adding a heavy attachment to hydrogens, discard the
      // hydrogen and start over
      combinedMol = nullptr;
      smilesVect.clear();
      attachments.clear();
      mols.clear();
    }
  }

  labelled = false;
  std::copy(rlabel_attachments.begin(), rlabel_attachments.end(),
            std::inserter(attachments, attachments.end()));

  mols.push_back(newMol);
  static const std::regex remove_isotopes_regex("\\[\\d*\\*\\]");
  // remove the isotope labels from the SMILES string to avoid
  // that identical R-group are perceived as different when
  // MCS alignment is not used (NoAlign flag)
  smilesVect.push_back(std::regex_replace(MolToSmiles(*newMol, true),
                                          remove_isotopes_regex, "*"));
  if (!combinedMol.get()) {
    combinedMol = boost::shared_ptr<RWMol>(new RWMol(*mols[0].get()));
  } else {
    ROMol *m = combineMols(*combinedMol.get(), *newMol.get());
    single_fragment = false;
    m->updateProps(*combinedMol.get());
    combinedMol.reset(new RWMol(*m));
    delete m;
  }
  smiles = getSmiles();
  combinedMol->setProp(common_properties::internalRgroupSmiles, smiles);
  computeIsHydrogen();
  is_linker = single_fragment && attachments.size() > 1;
}

std::map<int, int> RGroupData::getNumBondsToRlabels() const {
  std::map<int, int> rlabelsUsedCount;

  for (ROMol::AtomIterator atIt = combinedMol->beginAtoms();
       atIt != combinedMol->endAtoms(); ++atIt) {
    Atom *atom = *atIt;
    int rlabel;
    if (atom->getPropIfPresent<int>(RLABEL, rlabel)) {
      rlabelsUsedCount[rlabel] += 1;
    }
  }
  return rlabelsUsedCount;
}

std::string RGroupData::toString() const {
  auto attachmentString = std::accumulate(
      attachments.cbegin(), attachments.cend(), std::string(),
      [](std::string s, int a) {
        return s.empty() ? std::to_string(a)
                         : std::move(s) + ',' + std::to_string(a);
      });
  std::stringstream ss;
  ss << "RG " << attachmentString << " " << getSmiles();
  return ss.str();
}

void RGroupData::computeIsHydrogen() {  // is the rgroup all Hs
  for (const auto &mol : mols) {
    if (!isMolHydrogen(*mol)) {
      is_hydrogen = false;
      return;
    }
  }
  is_hydrogen = true;
}

bool RGroupData::isMolHydrogen(ROMol &mol) {
  for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    auto atom = *atIt;
    if (atom->getAtomicNum() > 1) {
      return false;
    } else if (atom->getAtomicNum() == 0 && !atom->hasProp(SIDECHAIN_RLABELS)) {
      return false;
    }
  }
  return true;
}

//! compute the canonical smiles for the attachments (bug: removes dupes since
//! we are using a set...)
std::string RGroupData::getSmiles() const {
  std::string s;
  for (const auto &it : smilesVect) {
    if (s.length()) {
      s += ".";
    }
    s += it;
  }
  return s;
}

}  // namespace RDKit
