//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Fragment.h"
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <boost/tokenizer.hpp>
#include <utility>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>

namespace RDKit {
namespace MolStandardize {

// constructor
FragmentRemover::FragmentRemover() {
  BOOST_LOG(rdInfoLog) << "Initializing FragmentRemover\n";
  FragmentCatalogParams fparams(defaultCleanupParameters.fragmentFile);
  //  unsigned int numfg = fparams->getNumFuncGroups();
  //  TEST_ASSERT(fparams->getNumFuncGroups() == 61);
  this->d_fcat = new FragmentCatalog(&fparams);
  this->LEAVE_LAST = true;
  this->SKIP_IF_ALL_MATCH = false;
}

// overloaded constructor
FragmentRemover::FragmentRemover(const std::string fragmentFile,
                                 bool leave_last, bool skip_if_all_match) {
  std::string fname = !fragmentFile.empty()
                          ? fragmentFile
                          : defaultCleanupParameters.fragmentFile;
  FragmentCatalogParams fparams(fname);
  this->d_fcat = new FragmentCatalog(&fparams);
  if (!this->d_fcat) {
    throw ValueErrorException(
        "could not open fragment catalog parameter file " + fname);
  }
  this->LEAVE_LAST = leave_last;
  this->SKIP_IF_ALL_MATCH = skip_if_all_match;
}

FragmentRemover::FragmentRemover(
    const std::vector<std::pair<std::string, std::string>> &data,
    bool leave_last, bool skip_if_all_match) {
  FragmentCatalogParams fparams(data);
  this->d_fcat = new FragmentCatalog(&fparams);
  if (!this->d_fcat) {
    throw ValueErrorException("could not process input data");
  }
  this->LEAVE_LAST = leave_last;
  this->SKIP_IF_ALL_MATCH = skip_if_all_match;
}

// overloaded constructor
FragmentRemover::FragmentRemover(std::istream &fragmentStream, bool leave_last,
                                 bool skip_if_all_match) {
  FragmentCatalogParams fparams(fragmentStream);
  this->d_fcat = new FragmentCatalog(&fparams);
  if (!this->d_fcat) {
    throw ValueErrorException("could not constract fragment catalog");
  }
  this->LEAVE_LAST = leave_last;
  this->SKIP_IF_ALL_MATCH = skip_if_all_match;
}

// Destructor
FragmentRemover::~FragmentRemover() { delete d_fcat; };

ROMol *FragmentRemover::remove(const ROMol &mol) {
  auto molcp = new RWMol(mol);
  removeInPlace(*molcp);
  return static_cast<ROMol *>(molcp);
}

void FragmentRemover::removeInPlace(RWMol &mol) {
  BOOST_LOG(rdInfoLog) << "Running FragmentRemover\n";
  PRECONDITION(this->d_fcat, "");
  const FragmentCatalogParams *fparams = this->d_fcat->getCatalogParams();

  PRECONDITION(fparams, "");

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  bool sanitizeFrags = false;
  // provides the list of atom numbers in each fragment
  std::vector<std::vector<int>> atomFragMapping;

  // track original fragment index with the fragment itself
  std::vector<std::pair<boost::shared_ptr<ROMol>, unsigned int>> frags;
  for (const auto &frag :
       MolOps::getMolFrags(mol, sanitizeFrags, nullptr, &atomFragMapping)) {
    frags.emplace_back(frag, frags.size());
  }

  for (auto &fgci : fgrps) {
    size_t oCount = frags.size();
    if (!oCount) {
      break;
    }
    auto tfrags = frags;
    // remove any fragments that this
    frags.erase(std::remove_if(
                    frags.begin(), frags.end(),
                    [&fgci](const std::pair<boost::shared_ptr<ROMol>,
                                            unsigned int> &frag) -> bool {
                      return fgci->getNumAtoms() == frag.first->getNumAtoms() &&
                             fgci->getNumBonds() == frag.first->getNumBonds() &&
                             SubstructMatch(*frag.first, *fgci).size() > 0;
                    }),
                frags.end());
    if (this->LEAVE_LAST && !this->SKIP_IF_ALL_MATCH && frags.empty()) {
      // All the remaining fragments match this pattern - leave them all
      frags = tfrags;
      break;
    }
    if (frags.size() != oCount) {
      BOOST_LOG(rdInfoLog) << "Removed fragment: "
                           << fgci->getProp<std::string>(
                                  common_properties::_Name)
                           << "\n";
    }
  }
  if (frags.empty()) {
    if (this->SKIP_IF_ALL_MATCH) {
      BOOST_LOG(rdInfoLog)
          << "All fragments matched; original molecule returned." << std::endl;
    } else {
      mol.beginBatchEdit();
      for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
        mol.removeAtom(i);
      }
      mol.commitBatchEdit();
    }
    return;
  }

  boost::dynamic_bitset<> atomsToRemove(mol.getNumAtoms());
  atomsToRemove.set();
  // loop over remaining fragments and track atoms we aren't keeping
  for (const auto &frag : frags) {
    unsigned int fragIdx = frag.second;
    for (auto atomIdx : atomFragMapping[fragIdx]) {
      atomsToRemove.set(atomIdx, false);
    }
  }
  // remove the atoms that need to go
  mol.beginBatchEdit();
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    if (atomsToRemove[i]) {
      mol.removeAtom(i);
    }
  }
  mol.commitBatchEdit();
}

bool isOrganic(const ROMol &mol, const std::vector<int> &indices) {
  // Returns true if fragment contains at least one carbon atom.
  for (auto idx : indices) {
    if (mol.getAtomWithIdx(idx)->getAtomicNum() == 6) {
      return true;
    }
  }
  return false;
}

LargestFragmentChooser::LargestFragmentChooser(
    const LargestFragmentChooser &other) {
  BOOST_LOG(rdInfoLog) << "Initializing LargestFragmentChooser\n";
  preferOrganic = other.preferOrganic;
  useAtomCount = other.useAtomCount;
  countHeavyAtomsOnly = other.countHeavyAtomsOnly;
}

ROMol *LargestFragmentChooser::choose(const ROMol &mol) const {
  auto res = new RWMol(mol);
  chooseInPlace(*res);
  // resanitize the molecule
  MolOps::sanitizeMol(*res);

  return static_cast<ROMol *>(res);
}

void LargestFragmentChooser::chooseInPlace(RWMol &mol) const {
  BOOST_LOG(rdInfoLog) << "Running LargestFragmentChooser\n";

  if (!mol.getNumAtoms()) {
    return;
  }

  std::vector<std::vector<int>> frags;
  MolOps::getMolFrags(mol, frags);
  if (frags.size() == 1) {
    // nothing to do
    return;
  }

  LargestFragmentChooser::Largest l;

  SmilesWriteParams ps;
  int bestFragment = -1;
  for (auto fidx = 0u; fidx < frags.size(); ++fidx) {
    const auto &frag = frags[fidx];
    std::string smiles = MolFragmentToSmiles(mol, ps, frag);
    BOOST_LOG(rdInfoLog) << "Fragment: " << smiles << "\n";
    bool organic = isOrganic(mol, frag);
    if (this->preferOrganic) {
      // Skip this fragment if not organic and we already have an organic
      // fragment as the largest so far
      if (bestFragment >= 0 && l.Organic && !organic) {
        continue;
      }
      // Reset largest if it wasn't organic and this fragment is organic
      // if largest and organic and not largest['organic']:
      if (bestFragment >= 0 && organic && !l.Organic) {
        bestFragment = -1;
      }
    }
    unsigned int numatoms = 0;
    if (this->useAtomCount) {
      for (const auto idx : frag) {
        ++numatoms;
        if (!this->countHeavyAtomsOnly) {
          numatoms += mol.getAtomWithIdx(idx)->getTotalNumHs();
        }
      }
      // Skip this fragment if fewer atoms than the largest
      if (bestFragment >= 0 && (numatoms < l.NumAtoms)) {
        continue;
      }
    }

    // Skip this fragment if equal number of atoms but weight is lower
    double weight = 0.0;
    for (auto idx : frag) {
      const auto atom = mol.getAtomWithIdx(idx);
      // it's not important to be perfect here
      weight += 100 * atom->getAtomicNum() + atom->getIsotope() -
                atom->getFormalCharge() * .1;
      if (!this->countHeavyAtomsOnly) {
        weight += atom->getTotalNumHs();
      }
    }

    if (bestFragment >= 0 && (!this->useAtomCount || numatoms == l.NumAtoms) &&
        (weight < l.Weight)) {
      continue;
    }

    // Skip this fragment if equal number of atoms and equal weight but smiles
    // comes last alphabetically
    if (bestFragment >= 0 && (!this->useAtomCount || numatoms == l.NumAtoms) &&
        (weight == l.Weight) && (smiles > l.Smiles)) {
      continue;
    }

    BOOST_LOG(rdInfoLog) << "New largest fragment: " << smiles << " ("
                         << numatoms << ")\n";
    // Otherwise this is the largest so far
    l.Smiles = smiles;
    bestFragment = fidx;
    l.NumAtoms = numatoms;
    l.Weight = weight;
    l.Organic = organic;
  }
  mol.beginBatchEdit();
  for (auto fi = 0; fi < static_cast<int>(frags.size()); ++fi) {
    if (fi == bestFragment) {
      continue;
    }
    for (auto i : frags[fi]) {
      mol.removeAtom(i);
    }
  }
  mol.commitBatchEdit();
}

LargestFragmentChooser::Largest::Largest() : Smiles(""), Fragment(nullptr) {}

LargestFragmentChooser::Largest::Largest(std::string &smiles,
                                         boost::shared_ptr<ROMol> fragment,
                                         unsigned int &numatoms, double &weight,
                                         bool &organic)
    : Smiles(smiles),
      Fragment(std::move(fragment)),
      NumAtoms(numatoms),
      Weight(weight),
      Organic(organic) {}

}  // namespace MolStandardize
}  // namespace RDKit
