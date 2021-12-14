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
#include <GraphMol/Descriptors/MolDescriptors.h>
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
      return new ROMol(mol);
    }
    return new ROMol();
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
  auto *removed = new RWMol(mol);
  removed->beginBatchEdit();
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    if (atomsToRemove[i]) {
      removed->removeAtom(i);
    }
  }
  removed->commitBatchEdit();
  return static_cast<ROMol *>(removed);
}

bool isOrganic(const ROMol &frag) {
  // Returns true if fragment contains at least one carbon atom.
  for (const auto at : frag.atoms()) {
    if (at->getAtomicNum() == 6) {
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

ROMol *LargestFragmentChooser::choose(const ROMol &mol) {
  BOOST_LOG(rdInfoLog) << "Running LargestFragmentChooser\n";

  if (!mol.getNumAtoms()) {
    return new ROMol(mol);
  }
  std::vector<boost::shared_ptr<ROMol>> frags = MolOps::getMolFrags(mol);
  LargestFragmentChooser::Largest l;

  for (const auto &frag : frags) {
    std::string smiles = MolToSmiles(*frag);
    BOOST_LOG(rdInfoLog) << "Fragment: " << smiles << "\n";
    bool organic = isOrganic(*frag);
    if (this->preferOrganic) {
      // Skip this fragment if not organic and we already have an organic
      // fragment as the largest so far
      if (l.Fragment != nullptr && l.Organic && !organic) {
        continue;
      }
      // Reset largest if it wasn't organic and this fragment is organic
      // if largest and organic and not largest['organic']:
      if (l.Fragment != nullptr && organic && !l.Organic) {
        l.Fragment = nullptr;
      }
    }
    unsigned int numatoms = 0;
    if (this->useAtomCount) {
      for (const auto at : frag->atoms()) {
        ++numatoms;
        if (!this->countHeavyAtomsOnly) {
          numatoms += at->getTotalNumHs();
        }
      }
      // Skip this fragment if fewer atoms than the largest
      if (l.Fragment != nullptr && (numatoms < l.NumAtoms)) {
        continue;
      }
    }

    // Skip this fragment if equal number of atoms but weight is lower
    double weight = Descriptors::calcExactMW(*frag);
    if (l.Fragment != nullptr &&
        (!this->useAtomCount || numatoms == l.NumAtoms) &&
        (weight < l.Weight)) {
      continue;
    }

    // Skip this fragment if equal number of atoms and equal weight but smiles
    // comes last alphabetically
    if (l.Fragment != nullptr &&
        (!this->useAtomCount || numatoms == l.NumAtoms) &&
        (weight == l.Weight) && (smiles > l.Smiles)) {
      continue;
    }

    BOOST_LOG(rdInfoLog) << "New largest fragment: " << smiles << " ("
                         << numatoms << ")\n";
    // Otherwise this is the largest so far
    l.Smiles = smiles;
    l.Fragment = frag;
    l.NumAtoms = numatoms;
    l.Weight = weight;
    l.Organic = organic;
  }

  return new ROMol(*(l.Fragment));
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
