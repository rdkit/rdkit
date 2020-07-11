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
#include <boost/scoped_ptr.hpp>
#include <set>
#include <vector>


namespace RDKit
{

//! A single rgroup attached to a given core.  
struct RGroupData {
  boost::shared_ptr<RWMol> combinedMol;
  std::vector<boost::shared_ptr<ROMol>> mols;  // All the mols in the rgroup
  std::set<std::string> smilesSet;             // used for rgroup equivalence
  std::string smiles;                          // smiles for all the mols in the rgroup (with attachments)
  std::set<int> attachments;                   // core attachment points
  bool is_hydrogen = false;
  bool single_fragment = true;
  bool multiple_attachments = false;
  bool is_linker = false;
  bool labelled = false;

 private:
  RGroupData(const RGroupData &rhs);

 public:
  RGroupData() : combinedMol(), mols(), smilesSet(), smiles(), attachments() {}

  void add(boost::shared_ptr<ROMol> newMol,
           const std::vector<int> &rlabel_attachments) {
    // some fragments can be add multiple times if they are cyclic
    for (auto &mol : mols) {
      if (newMol.get() == mol.get()) {
        return;
      }
    }

    labelled = false;
    std::copy(rlabel_attachments.begin(), rlabel_attachments.end(),
              std::inserter(attachments, attachments.end()));

    mols.push_back(newMol);
    std::string smi = MolToSmiles(*newMol, true);
    // REVIEW: we probably shouldn't be using a set here... the merging of
    // duplicates is likely not what we want
    smilesSet.insert(smi);
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

 private:
  void computeIsHydrogen() {  // is the rgroup all Hs
    for (const auto &mol : mols) {
      for (ROMol::AtomIterator atIt = mol->beginAtoms();
           atIt != mol->endAtoms(); ++atIt) {
        if ((*atIt)->getAtomicNum() > 1) {
	  is_hydrogen = false;
	  return;
        }
      }
    }
    is_hydrogen = true;
  }

  //! compute the canonical smiles for the attachments (bug: removes dupes since we are using a set...)
  std::string getSmiles() const {
    std::string s;
    for (const auto &it : smilesSet) {
      if (s.length()) {
        s += ".";
      }
      s += it;
    }
    return s;
  }
};
}

#endif
