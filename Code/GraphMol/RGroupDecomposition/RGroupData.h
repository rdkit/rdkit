//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_DATA
#define RGROUP_DATA

#include "../RDKitBase.h"
#include "RGroupUtils.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <DataStructs/ExplicitBitVect.h>
#include <boost/scoped_ptr.hpp>
#include <set>
#include <vector>
#include <regex>

namespace RDKit {

//! A single rgroup attached to a given core.
struct RGroupData {
  boost::shared_ptr<RWMol> combinedMol;
  std::vector<boost::shared_ptr<ROMol>> mols;  // All the mols in the rgroup
  std::vector<std::string> smilesVect;         // used for rgroup equivalence
  std::string
      smiles;  // smiles for all the mols in the rgroup (with attachments)
  std::set<int> attachments;  // core attachment points
  std::unique_ptr<ExplicitBitVect>
      fingerprint;  // fingerprint for score calculations
  std::vector<int> fingerprintOnBits;
  bool is_hydrogen = false;
  bool single_fragment = true;
  bool multiple_attachments = false;
  bool is_linker = false;
  bool labelled = false;

 private:
  RGroupData(const RGroupData &rhs);

 public:
  RGroupData() {}

  void add(boost::shared_ptr<ROMol> newMol,
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

  std::map<int, int> getNumBondsToRlabels() const {
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

  std::string toString() const {
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

 private:
  void computeIsHydrogen() {  // is the rgroup all Hs
    for (const auto &mol : mols) {
      if (!isMolHydrogen(*mol)) {
        is_hydrogen = false;
        return;
      }
    }
    is_hydrogen = true;
  }

  bool isMolHydrogen(ROMol &mol) {
    for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      auto atom = *atIt;
      if (atom->getAtomicNum() > 1) {
        return false;
      } else if (atom->getAtomicNum() == 0 &&
                 !atom->hasProp(SIDECHAIN_RLABELS)) {
        return false;
      }
    }
    return true;
  }

  //! compute the canonical smiles for the attachments (bug: removes dupes since
  //! we are using a set...)
  std::string getSmiles() const {
    std::string s;
    for (const auto &it : smilesVect) {
      if (s.length()) {
        s += ".";
      }
      s += it;
    }
    return s;
  }
};
}  // namespace RDKit

#endif
